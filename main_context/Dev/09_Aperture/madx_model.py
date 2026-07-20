"""
madx_model.py

Thin cpymad/MAD-X execution layer for the ISIS RCS optics GUI backend.

Responsibilities
----------------
- Start a MAD-X/cpymad instance.
- Load the base ISIS lattice files.
- Apply one generated machine-state override file.
- Use the requested sequence.
- Run TWISS.
- Return Twiss and summary information as pandas DataFrames.

Important
---------
This layer deliberately does not load 2023.strength.

The operational state should now be supplied by the generated machine-state
override file written by machine_state_writer.py. That override captures:
    - beam momentum / Brho
    - main magnet scaling
    - trim quadrupole settings
    - harmonic tune variables
    - corrector kicks
"""

import os
from datetime import datetime
from pathlib import Path

import pandas as pd

from tune_control import strength_to_current


try:
    from cpymad.madx import Madx
except ImportError as exc:
    Madx = None
    _CPYMAD_IMPORT_ERROR = exc
else:
    _CPYMAD_IMPORT_ERROR = None


class MadxModel:
    """
    Thin wrapper around a single cpymad MAD-X instance.

    Parameters
    ----------
    lattice_folder : str
        Folder containing the MAD-X lattice files.

    sequence_name : str
        Name of the MAD-X sequence to use.

    machine_state_file : str, optional
        Generated machine-state override file. This should normally be produced
        by write_machine_state_file(...).

    output_dir : str
        Directory for logs and optional TWISS output files.

    injected_beam_file, strength_file, elements_file, sequence_file : str
        Standard lattice file names inside lattice_folder.

    logfile : str
        Name of the cpymad log file written inside output_dir.

    echo, warn : bool
        MAD-X option flags.
    """

    def __init__(
        self,
        lattice_folder,
        sequence_name="synchrotron",
        machine_state_file=None,
        output_dir="./madx_runs",
        strength_file="ISIS.strength",
        elements_file="ISIS.elements",
        sequence_file="ISIS.sequence",
        aperture_file=None,
        logfile="cpymad_logfile.log",
        echo=True,
        warn=True,
    ):
        self.lattice_folder = lattice_folder
        self.sequence_name = sequence_name
        self.machine_state_file = machine_state_file
        self.output_dir = output_dir

        self.strength_file = strength_file
        self.elements_file = elements_file
        self.sequence_file = sequence_file
        self.aperture_file = aperture_file

        self.logfile = logfile
        self.echo = echo
        self.warn = warn

        self.madx = None
        self.logfile_path = None
        self.twiss_df = None
        self.summary_df = None
        self.aperture_df = None

        self.loaded_files = []
        self.metadata = {
            "created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "lattice_folder": self.lattice_folder,
            "sequence_name": self.sequence_name,
            "machine_state_file": self.machine_state_file,
            "output_dir": self.output_dir,
            "strength_file": self.strength_file,
            "elements_file": self.elements_file,
            "sequence_file": self.sequence_file,
            "aperture_file": self.aperture_file,
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _require_cpymad(self):
        if Madx is None:
            raise ImportError(
                "cpymad could not be imported. Install cpymad before using "
                "MadxModel."
            ) from _CPYMAD_IMPORT_ERROR

    def _full_lattice_path(self, filename):
        return os.path.join(self.lattice_folder, filename)

    def _check_file_exists(self, path, label=None):
        if not os.path.isfile(path):
            if label is None:
                label = "file"
            raise FileNotFoundError(f"Missing {label}: {path}")

    def _call_file(self, path, label=None):
        self._check_file_exists(path, label=label)
        self.madx.call(file=path)
        self.loaded_files.append(path)
        return path

    def _ensure_started(self):
        if self.madx is None:
            self.start()
        return self.madx

    @staticmethod
    def _format_madx_value(value):
        if value is None:
            return "0.0"
        return f"{float(value):.16g}"

    def _assign_global(self, name, value, assignment=":="):
        self.madx.input(
            f"{name} {assignment} {self._format_madx_value(value)};"
        )

    def _get_global_float(self, name):
        try:
            return float(self.madx.globals[name])
        except Exception:
            return float(getattr(self.madx.globals, name))

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self):
        """
        Start a fresh cpymad MAD-X instance.
        """

        self._require_cpymad()

        os.makedirs(self.output_dir, exist_ok=True)

        self.logfile_path = os.path.join(self.output_dir, self.logfile)
        logfile_handle = open(self.logfile_path, "w")

        self.madx = Madx(stdout=logfile_handle)
        self.madx.options.echo = bool(self.echo)
        self.madx.options.warn = bool(self.warn)

        self.metadata["logfile_path"] = self.logfile_path
        self.metadata["started"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        return self.madx

    def reset(self):
        """
        Reset the MAD-X instance by starting a fresh one.

        This is safer than trying to undo arbitrary MAD-X global state.
        """

        self.madx = None
        self.twiss_df = None
        self.summary_df = None
        self.aperture_df = None
        self.loaded_files = []

        return self.start()

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    def load_lattice(
        self,
        machine_state_file=None,
        aperture_file=None,
        use_sequence=True,
    ):
        """
        Load the ISIS lattice and optional machine-state override.

        Load order
        ----------
        1. ISIS.elements
        2. ISIS.sequence
        3. ISIS.strength
        4. aperture file, if supplied
        5. generated machine-state override file, if supplied

        Notes
        -----
        2023.strength is intentionally not loaded here.

        ISIS.injected_beam is intentionally not loaded here. The beam momentum
        should be supplied by the generated machine-state file, using the
        selected cycle-time beam state.

        If an aperture file is supplied, it is loaded after the sequence and
        base strength file, but before the machine-state override.
        """

        self._ensure_started()

        if machine_state_file is not None:
            self.machine_state_file = machine_state_file
            self.metadata["machine_state_file"] = machine_state_file

        if aperture_file is not None:
            self.aperture_file = aperture_file
            self.metadata["aperture_file"] = aperture_file

        elements_path = self._full_lattice_path(self.elements_file)
        sequence_path = self._full_lattice_path(self.sequence_file)
        strength_path = self._full_lattice_path(self.strength_file)

        aperture_path = None
        if self.aperture_file is not None:
            aperture_path = self._full_lattice_path(self.aperture_file)

        self._call_file(elements_path, label="elements file")
        self._call_file(sequence_path, label="sequence file")
        self._call_file(strength_path, label="base strength file")

        if aperture_path is not None:
            self._call_file(aperture_path, label="aperture file")

        if self.machine_state_file is not None:
            self.apply_machine_state_file(self.machine_state_file)

        if use_sequence:
            self.use_sequence()

        self.metadata["loaded_lattice"] = True
        self.metadata["loaded_files"] = list(self.loaded_files)

        return self.madx

    def apply_machine_state_file(self, machine_state_file):
        """
        Apply a generated machine-state override file.

        This is expected to be a MAD-X-readable file produced by
        machine_state_writer.py.
        """

        self._ensure_started()

        self._check_file_exists(machine_state_file, label="machine-state file")

        self.madx.call(file=machine_state_file)
        self.loaded_files.append(machine_state_file)

        self.machine_state_file = machine_state_file
        self.metadata["machine_state_file"] = machine_state_file

        return machine_state_file

    def apply_machine_state(self, machine_state, include_trim_quads=True):
        """
        Apply a MachineState object directly to the active MAD-X instance.
        """

        self._ensure_started()

        beam = machine_state.beam_summary_dict()

        self._assign_global("brho", beam["brho_Tm"])
        self.madx.input(
            "beam, particle=proton, "
            f"pc={self._format_madx_value(beam['momentum_GeV_c'])};"
        )

        for name, value in machine_state.main_magnet_scaling.items():
            self._assign_global(name, value)

        if include_trim_quads:
            self._assign_global("kqtd", machine_state.kqtd)
            self._assign_global("kqtf", machine_state.kqtf)

        for name, value in machine_state.harmonic_tunes.items():
            self._assign_global(name, value, assignment="=")

        for name, value in machine_state.hd_corrector_kicks_rad.items():
            self._assign_global(name, value)

        for name, value in machine_state.vd_corrector_kicks_rad.items():
            self._assign_global(name, value)

        self.metadata["machine_state_applied"] = True
        self.metadata["machine_state_cycle_time_ms"] = machine_state.cycle_time_ms

        return machine_state

    def match_tune_from_machine_state(
        self,
        machine_state,
        requested_dq1=None,
        requested_dq2=None,
        step=1.0e-4,
        calls=50000,
        tolerance=1.0e-6,
        chrom=True,
        apply_machine_state=True,
        include_trim_quads=True,
        run_twiss_after_match=True,
    ):
        """
        Match the tune in MAD-X using the requested tunes in MachineState.

        The matched kqtf/kqtd values and equivalent currents are written back
        into the supplied MachineState object.
        """

        self._ensure_started()

        if not self.metadata.get("loaded_lattice", False):
            raise RuntimeError(
                "Load the lattice before matching tunes: call load_lattice() first."
            )

        if machine_state.tune_method != "madx_match":
            raise ValueError(
                "match_tune_from_machine_state requires "
                "machine_state.tune_method == 'madx_match'."
            )

        if apply_machine_state:
            self.apply_machine_state(
                machine_state,
                include_trim_quads=include_trim_quads,
            )

        self.use_sequence()

        self.madx.command.match(chrom=chrom)
        self.madx.command.vary(name="kqtd", step=float(step))
        self.madx.command.vary(name="kqtf", step=float(step))

        if requested_dq1 is None or requested_dq2 is None:
            self.madx.command.global_(
                sequence=self.sequence_name,
                q1=float(machine_state.requested_qx),
                q2=float(machine_state.requested_qy),
            )
        else:
            self.madx.command.global_(
                sequence=self.sequence_name,
                q1=float(machine_state.requested_qx),
                q2=float(machine_state.requested_qy),
                dq1=float(requested_dq1),
                dq2=float(requested_dq2),
            )

        self.madx.command.jacobian(
            calls=int(calls),
            tolerance=float(tolerance),
        )
        self.madx.command.endmatch()

        if run_twiss_after_match:
            self.madx.twiss(sequence=self.sequence_name)
            self.summary_df = self.get_summary_df(refresh=True)

        kqtd = self._get_global_float("kqtd")
        kqtf = self._get_global_float("kqtf")

        beam_state = machine_state.beam_state
        brho_Tm = float(beam_state.brho_Tm)
        pn = float(beam_state.normalised_momentum)

        iqtd_A = strength_to_current(
            strength=kqtd,
            brho_Tm=brho_Tm,
            pn=pn,
        )
        iqtf_A = strength_to_current(
            strength=kqtf,
            brho_Tm=brho_Tm,
            pn=pn,
        )

        summary = {}
        try:
            summary = self.get_summary_dict()
        except RuntimeError:
            pass

        result = {
            "kqtd": float(kqtd),
            "kqtf": float(kqtf),
            "iqtd_A": float(iqtd_A),
            "iqtf_A": float(iqtf_A),
            "matched_qx": summary.get("q1", None),
            "matched_qy": summary.get("q2", None),
            "matched_dqx": summary.get("dq1", None),
            "matched_dqy": summary.get("dq2", None),
        }

        machine_state.kqtd = result["kqtd"]
        machine_state.kqtf = result["kqtf"]
        machine_state.iqtd_A = result["iqtd_A"]
        machine_state.iqtf_A = result["iqtf_A"]
        machine_state.metadata["madx_match"] = {
            "sequence_name": self.sequence_name,
            "requested_qx": machine_state.requested_qx,
            "requested_qy": machine_state.requested_qy,
            "requested_dq1": requested_dq1,
            "requested_dq2": requested_dq2,
            "matched_qx": result["matched_qx"],
            "matched_qy": result["matched_qy"],
            "matched_dqx": result["matched_dqx"],
            "matched_dqy": result["matched_dqy"],
        }

        self.metadata["madx_match"] = dict(machine_state.metadata["madx_match"])

        return result

    def match_orbit(
        self,
        constraints,
        vary_names,
        horizontal=True,
        step=1.0e-4,
        calls=50000,
        tolerance=1.0e-6,
        chrom=False,
        max_kick=None,
        run_twiss_after_match=True,
        twiss_columns=None,
    ):
        """
        Match one closed-orbit plane by varying selected MAD-X kick globals.

        This is intentionally low-level: callers are responsible for building
        measurement constraints and choosing the correct plane's vary list.
        """

        self._ensure_started()

        if not self.metadata.get("loaded_lattice", False):
            raise RuntimeError(
                "Load the lattice before matching orbit: call load_lattice() first."
            )

        constraints = list(constraints or [])
        vary_names = [str(name) for name in list(vary_names or [])]

        if not constraints:
            raise ValueError("At least one orbit constraint is required.")

        if not vary_names:
            raise ValueError("At least one MAD-X variable must be selected to vary.")

        plane_key = "x" if horizontal else "y"
        other_plane_key = "y" if horizontal else "x"
        filtered_constraints = []
        for constraint in constraints:
            filtered = dict(constraint)
            filtered.pop(other_plane_key, None)
            if plane_key in filtered:
                filtered_constraints.append(filtered)

        if not filtered_constraints:
            raise ValueError(f"No {plane_key!r} constraints available for orbit match.")

        self.madx.command.match(chrom=chrom)

        for name in vary_names:
            if max_kick is None:
                self.madx.command.vary(name=name, step=float(step))
            else:
                self.madx.command.vary(
                    name=name,
                    step=float(step),
                    lower=-float(max_kick),
                    upper=float(max_kick),
                )

        for constraint in filtered_constraints:
            self.madx.command.constraint(**constraint)

        self.madx.command.jacobian(
            calls=int(calls),
            tolerance=float(tolerance),
        )
        self.madx.command.endmatch()

        twiss_df = None
        summary = {}
        if run_twiss_after_match:
            twiss_df = self.run_twiss(columns=twiss_columns)
            summary = self.get_summary_dict()

        matched_kicks = {}
        for name in vary_names:
            try:
                matched_kicks[name] = self._get_global_float(name)
            except Exception:
                matched_kicks[name] = None

        self.metadata["orbit_match"] = {
            "sequence_name": self.sequence_name,
            "plane": "H" if horizontal else "V",
            "n_constraints": len(filtered_constraints),
            "n_vary": len(vary_names),
            "vary_names": list(vary_names),
            "max_kick": max_kick,
            "calls": int(calls),
            "tolerance": float(tolerance),
        }

        return {
            "plane": "H" if horizontal else "V",
            "constraints": filtered_constraints,
            "vary_names": list(vary_names),
            "matched_kicks_rad": matched_kicks,
            "summary": summary,
            "twiss_df": twiss_df,
            "metadata": dict(self.metadata["orbit_match"]),
        }

    def correct_orbit(
        self,
        plane,
        corrector_names,
        monitor_pattern=None,
        monitor_names=None,
        model_table="bare",
        output_prefix=None,
        mode="svd",
        cond=1,
        ncorr=0,
        error=1.0e-7,
        corzero=1,
        monerror=0,
        monscale=0,
        flag="ring",
        run_twiss_after_correct=True,
        corrected_table="corrected",
        twiss_columns=None,
    ):
        """
        Run MAD-X CORRECT for one plane using selected steering correctors.

        The caller must create the reference TWISS table named by model_table
        before calling this method. This wrapper only owns the MAD-X correction
        command and the usekick/usemonitor selection needed by the GUI layer.
        """

        self._ensure_started()

        if not self.metadata.get("loaded_lattice", False):
            raise RuntimeError(
                "Load the lattice before correcting orbit: call load_lattice() first."
            )

        plane_key = str(plane).upper()
        if plane_key in ("H", "X", "HORIZONTAL"):
            madx_plane = "x"
            default_monitor_pattern = ".*HM.*"
        elif plane_key in ("V", "Y", "VERTICAL"):
            madx_plane = "y"
            default_monitor_pattern = ".*VM.*"
        else:
            raise ValueError("plane must be H/X/horizontal or V/Y/vertical.")

        corrector_names = [str(name) for name in list(corrector_names or [])]
        if not corrector_names:
            raise ValueError("At least one corrector must be selected.")

        if output_prefix is None:
            output_prefix = Path(self.output_dir) / f"orbit_correct_{madx_plane}"
        output_prefix = Path(output_prefix)
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        clist_path = output_prefix.with_name(output_prefix.name + "_clist.tfs")
        mlist_path = output_prefix.with_name(output_prefix.name + "_mlist.tfs")
        # MAD-X lower-cases command input, including file paths.
        clist_path = Path(str(clist_path).lower())
        mlist_path = Path(str(mlist_path).lower())
        clist_path.parent.mkdir(parents=True, exist_ok=True)
        mlist_path.parent.mkdir(parents=True, exist_ok=True)

        self.madx.command.usekick(
            sequence=self.sequence_name,
            status="off",
            pattern=".*",
        )
        for name in corrector_names:
            element_name = name.replace("_kick", "").upper()
            self.madx.command.usekick(
                sequence=self.sequence_name,
                status="on",
                pattern=f".*{element_name}$",
            )

        self.madx.command.usemonitor(
            sequence=self.sequence_name,
            status="off",
            pattern=".*",
        )
        if monitor_names is not None:
            for monitor_name in monitor_names:
                self.madx.command.usemonitor(
                    sequence=self.sequence_name,
                    status="on",
                    pattern=f".*{str(monitor_name).upper()}$",
                )
        else:
            self.madx.command.usemonitor(
                sequence=self.sequence_name,
                status="on",
                pattern=monitor_pattern or default_monitor_pattern,
            )

        command_kwargs = {
            "model": str(model_table),
            "sequence": self.sequence_name,
            "plane": madx_plane,
            "flag": flag,
            "error": float(error),
            "mode": str(mode),
            "corzero": int(corzero),
            "monerror": int(monerror),
            "monscale": int(monscale),
            "clist": str(clist_path),
            "mlist": str(mlist_path),
        }
        if str(mode).lower() == "svd":
            command_kwargs["cond"] = cond
        elif str(mode).lower() == "micado":
            command_kwargs["ncorr"] = int(ncorr)

        self.madx.command.correct(**command_kwargs)

        twiss_df = None
        summary = {}
        if run_twiss_after_correct:
            self.madx.twiss(sequence=self.sequence_name, table=corrected_table)
            table = getattr(self.madx.table, corrected_table)
            twiss_df = table.dframe().reset_index(drop=True)
            twiss_df = self._normalise_dataframe_columns(twiss_df)
            if twiss_columns is not None:
                twiss_df = twiss_df[twiss_columns]
            self.summary_df = self.get_summary_df(refresh=True)
            summary = self.get_summary_dict()

        result = {
            "plane": "H" if madx_plane == "x" else "V",
            "madx_plane": madx_plane,
            "corrector_names": list(corrector_names),
            "monitor_pattern": monitor_pattern or default_monitor_pattern,
            "monitor_names": None if monitor_names is None else list(monitor_names),
            "model_table": str(model_table),
            "clist_path": str(clist_path),
            "mlist_path": str(mlist_path),
            "summary": summary,
            "twiss_df": twiss_df,
        }

        self.metadata["orbit_correct"] = dict(result)
        self.metadata["orbit_correct"].pop("twiss_df", None)

        return result

    def apply_error_table(self, error_table_path, table_name="error_table"):
        """
        Apply a MAD-X error table before running TWISS.

        Parameters
        ----------
        error_table_path : str
            Path to a MAD-X-readable TFS error table.

        table_name : str
            Name used internally by MAD-X for the loaded table.

        Notes
        -----
        This expects a standard MAD-X error table suitable for:

            readtable, file=..., table=...;
            seterr, table=...;

        Important
        ---------
        MAD-X clears applied errors when USE is called. Therefore the expected
        order is:

            load_lattice(use_sequence=True)
            apply_error_table(...)
            run_twiss(...)

        run_twiss() must not call USE after this method.
        """

        self._ensure_started()

        self._check_file_exists(error_table_path, label="MAD-X error table")

        self.madx.input(
            f'readtable, file="{error_table_path}", table={table_name};'
        )

        self.madx.input(
            f"seterr, table={table_name};"
        )

        self.loaded_files.append(error_table_path)
        self.metadata["error_table_path"] = error_table_path
        self.metadata["error_table_name"] = table_name
        self.metadata["errors_applied"] = True

        return error_table_path
        
    def use_sequence(self, sequence_name=None):
        """
        Select the active MAD-X sequence.
        """

        self._ensure_started()

        if sequence_name is not None:
            self.sequence_name = sequence_name
            self.metadata["sequence_name"] = sequence_name

        self.madx.use(sequence=self.sequence_name)

        return self.sequence_name

    # ------------------------------------------------------------------
    # TWISS
    # ------------------------------------------------------------------

    def run_twiss(
        self,
        sequence_name=None,
        file_out=None,
        save_twiss=False,
        columns=None,
    ):
        """
        Run MAD-X TWISS and return the Twiss table as a pandas DataFrame.

        Parameters
        ----------
        sequence_name : str, optional
            Sequence name to pass to MAD-X TWISS. If None, self.sequence_name is
            used.

            Important: this method does not call USE. If a different sequence is
            supplied here, it only changes the sequence argument used by TWISS.

        file_out : str, optional
            Output TFS filename. If None and save_twiss=True, a default file is
            written inside output_dir.

        save_twiss : bool
            If True, write the TWISS table to file.

        columns : list[str], optional
            Columns to request from MAD-X. If None, a standard optics/orbit set
            is requested.

        Important
        ---------
        This method deliberately does not call self.use_sequence(...).

        MAD-X clears applied errors when USE is called. Therefore the safe error
        table workflow is:

            model.load_lattice(use_sequence=True)
            model.apply_error_table(error_table_path)
            twiss_df = model.run_twiss()
        """

        self._ensure_started()

        if sequence_name is not None:
            self.sequence_name = sequence_name
            self.metadata["sequence_name"] = sequence_name

        if columns is None:
            columns = self.default_twiss_columns()

        column_string = ", ".join(columns)

        self.madx.input('set, format="12.12f";')
        self.madx.input("select, flag=twiss, clear;")
        self.madx.input(f"select, flag=twiss, column={column_string};")

        if save_twiss:
            if file_out is None:
                file_out = os.path.join(
                    self.output_dir,
                    f"{self.sequence_name}_twiss.tfs",
                )

            self.madx.twiss(sequence=self.sequence_name, file=file_out)
            self.metadata["twiss_file"] = file_out

        else:
            self.madx.twiss(sequence=self.sequence_name)

        self.twiss_df = self.madx.table.twiss.dframe()
        self.twiss_df = self._normalise_dataframe_columns(self.twiss_df)

        self.summary_df = self.get_summary_df(refresh=True)

        self.metadata["last_twiss"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.metadata["twiss_rows"] = len(self.twiss_df)

        return self.twiss_df

    def run_aperture(
        self,
        file_out=None,
        dqf=3.24,
        betaqfx=18.6669,
        interval=0.1,
        columns=None,
    ):
        """
        Run MAD-X APERTURE and return the written table as a DataFrame.

        The table is read back from the written TFS file because direct cpymad
        table decoding can fail for some aperture row names.
        """

        self._ensure_started()

        if columns is None:
            columns = self.default_aperture_columns()

        if file_out is None:
            file_out = os.path.join(
                self.output_dir,
                f"{self.sequence_name}_madx_aperture.tfs",
            )

        output_dir = os.path.dirname(str(file_out))
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        if self.twiss_df is None:
            self.madx.twiss(sequence=self.sequence_name)
            self.summary_df = self.get_summary_df(refresh=True)

        column_string = ", ".join(columns)

        self.madx.input("select, flag=aperture, clear;")
        self.madx.input(f"select, flag=aperture, column={column_string};")
        self.madx.input('set, format="12.12f";')
        self.madx.input(
            f'aperture, range=#s/#e, dqf={float(dqf)}, '
            f'betaqfx={float(betaqfx)}, interval={float(interval)}, '
            f'file="{file_out}";'
        )

        actual_file_out = file_out
        if not os.path.exists(actual_file_out):
            lower_file_out = str(file_out).lower()
            if os.path.exists(lower_file_out):
                actual_file_out = lower_file_out

        aperture_df = self._read_tfs_table(actual_file_out)
        aperture_df = self._normalise_dataframe_columns(aperture_df)

        self.aperture_df = aperture_df
        self.metadata["aperture_file_out"] = str(actual_file_out)
        self.metadata["aperture_rows"] = len(aperture_df)
        self.metadata["aperture_interval"] = float(interval)

        return aperture_df

    @staticmethod
    def default_aperture_columns():
        """
        Standard MAD-X APERTURE columns for aperture-margin workflows.
        """

        return [
            "name",
            "n1",
            "n1x_m",
            "n1y_m",
            "apertype",
            "aper_1",
            "aper_2",
            "aper_3",
            "aper_4",
            "rtol",
            "xtol",
            "ytol",
            "s",
            "betx",
            "bety",
            "dx",
            "dy",
            "x",
            "y",
            "px",
            "py",
            "on_ap",
            "on_elem",
            "spec",
            "x_pos_hit",
            "y_pos_hit",
        ]

    @staticmethod
    def _read_tfs_table(path):
        """
        Read a MAD-X TFS file using tfs if available, with a simple fallback.
        """

        try:
            import tfs

            return tfs.read(path)
        except Exception:
            pass

        with open(path, "r", errors="replace") as handle:
            lines = handle.readlines()

        header_line = None
        data_lines = []

        for line in lines:
            stripped = line.strip()
            if stripped.startswith("*"):
                header_line = stripped
            elif stripped and not stripped.startswith("@") and not stripped.startswith("$"):
                data_lines.append(stripped)

        if header_line is None:
            raise RuntimeError(f"Could not find TFS header line in {path}")

        columns = header_line.split()[1:]
        rows = []
        for line in data_lines:
            parts = line.split()
            if len(parts) == len(columns):
                rows.append(parts)

        df = pd.DataFrame(rows, columns=columns)
        for col in df.columns:
            if col.lower() not in ("name", "apertype", "on_ap", "on_elem", "spec"):
                df[col] = pd.to_numeric(df[col], errors="coerce")

        return df

    @staticmethod
    def default_twiss_columns():
        """
        Standard TWISS columns for optics GUI development.
        """

        return [
            "keyword",
            "name",
            "s",
            "l",
            "betx",
            "alfx",
            "mux",
            "bety",
            "alfy",
            "muy",
            "x",
            "px",
            "y",
            "py",
            "t",
            "pt",
            "dx",
            "dpx",
            "dy",
            "dpy",
            "wx",
            "phix",
            "dmux",
            "wy",
            "phiy",
            "dmuy",
            "ddx",
            "ddpx",
            "ddy",
            "ddpy",
            "r11",
            "r12",
            "r21",
            "r22",
            "energy",
            "angle",
            "k0l",
            "k0sl",
            "k1l",
            "k1sl",
            "k2l",
            "k2sl",
            "hkick",
            "vkick",
            "tilt",
        ]

    @staticmethod
    def _normalise_dataframe_columns(df):
        """
        Return a copy with lowercase column names.

        This keeps downstream plotting and processing consistent with the
        existing helper style.
        """

        df_out = df.copy()
        df_out.columns = [str(col).lower() for col in df_out.columns]
        return df_out

    # ------------------------------------------------------------------
    # Results
    # ------------------------------------------------------------------

    def get_twiss_df(self):
        """
        Return the most recently calculated TWISS DataFrame.
        """

        if self.twiss_df is None:
            raise RuntimeError("No TWISS table available. Run run_twiss() first.")

        return self.twiss_df

    def get_summary_dict(self):
        """
        Return MAD-X summary table as a dictionary.
        """

        self._ensure_started()

        if "summ" not in list(self.madx.table):
            raise RuntimeError("MAD-X summary table is not available yet.")

        summary = {}

        for key, value in self.madx.table.summ.items():
            try:
                summary[key.lower()] = value[0]
            except Exception:
                summary[key.lower()] = value

        return summary

    def get_summary_df(self, refresh=False):
        """
        Return MAD-X summary table as a one-row pandas DataFrame.
        """

        if self.summary_df is None or refresh:
            summary = self.get_summary_dict()
            self.summary_df = pd.DataFrame([summary])

        return self.summary_df

    def get_tunes(self):
        """
        Return Qx, Qy if available from the MAD-X summary table.
        """

        summary = self.get_summary_dict()

        qx = summary.get("q1", None)
        qy = summary.get("q2", None)

        return qx, qy

    def get_chromaticities(self):
        """
        Return dQx, dQy if available from the MAD-X summary table.
        """

        summary = self.get_summary_dict()

        dqx = summary.get("dq1", None)
        dqy = summary.get("dq2", None)

        return dqx, dqy

    def get_metadata(self):
        """
        Return run metadata.
        """

        metadata = dict(self.metadata)
        metadata["loaded_files"] = list(self.loaded_files)
        return metadata

    # ------------------------------------------------------------------
    # Convenience checks
    # ------------------------------------------------------------------

    def print_loaded_files(self):
        """
        Print the files loaded into MAD-X in order.
        """

        for path in self.loaded_files:
            print(path)

    def print_summary(self):
        """
        Print a compact MAD-X summary.
        """

        summary = self.get_summary_dict()

        qx = summary.get("q1", None)
        qy = summary.get("q2", None)
        dqx = summary.get("dq1", None)
        dqy = summary.get("dq2", None)

        print(f"Sequence: {self.sequence_name}")
        print(f"Qx = {qx}")
        print(f"Qy = {qy}")
        print(f"dQx = {dqx}")
        print(f"dQy = {dqy}")

    # ------------------------------------------------------------------
    # Simple one-shot workflow
    # ------------------------------------------------------------------

    def run(
        self,
        machine_state_file=None,
        error_table_path=None,
        error_table_name="error_table",
        save_twiss=False,
        twiss_file_out=None,
    ):
        """
        One-shot convenience method.

        Equivalent to:

            load_lattice(...)
            optionally apply_error_table(...)
            run_twiss(...)

        Notes
        -----
        The order is important:

            1. load_lattice(use_sequence=True)
            2. apply_error_table(...)
            3. run_twiss(...)

        run_twiss() deliberately does not call USE, so applied errors are
        preserved.
        """

        self.load_lattice(
            machine_state_file=machine_state_file,
            use_sequence=True,
        )

        if error_table_path is not None:
            self.apply_error_table(
                error_table_path=error_table_path,
                table_name=error_table_name,
            )

        return self.run_twiss(
            save_twiss=save_twiss,
            file_out=twiss_file_out,
        )
