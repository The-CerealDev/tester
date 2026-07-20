"""
snapshot.py

Whole-system snapshot orchestration for the ISIS RCS optics GUI backend.

This layer assembles existing lower-layer outputs for one or more machine
states. It does not own physics calculations, MAD-X execution details, or
plotting. The outputs are DataFrame-backed and ready for notebooks, plotting
helpers and the future GUI.
"""

from dataclasses import asdict, dataclass, field, replace
from datetime import datetime
from pathlib import Path
import copy
import re

import pandas as pd

from .aperture import (
    DEFAULT_SOURCE_APERTURE_CSV,
    evaluate_aperture_margins,
    normalise_madx_aperture_table,
    read_source_aperture_csv,
)
from .cycle_time import RCSRamp
from .envelope import EnvelopeInputs, evaluate_envelope_from_twiss
from .machine_state import MachineState
from .madx_model import MadxModel
from .orbit_branch import (
    OrbitBranch,
    OrbitBranchConfig,
    compare_orbit_results,
    extract_orbit_df,
    summarise_orbit_df,
)
from .orbit_correction import correct_bpm_measurements_with_madx_correct, orbit_residuals
from .tune import (
    build_tune_programme_table,
    build_working_point_table,
    extract_tune_summary,
    make_tune_diagram_inputs,
)
from .tune_matching import normalise_harmonic_inputs


TIMESTAMP_FORMAT = "%Y%m%d_%H%M%S"
WARNING_COLUMNS = ["timestamp", "snapshot_id", "severity", "source", "code", "message", "context"]
PLOT_MANIFEST_COLUMNS = [
    "timestamp",
    "snapshot_id",
    "aspect",
    "plot_name",
    "path",
    "dpi",
    "title",
    "has_title",
    "has_axis_labels",
    "has_legend",
]


@dataclass
class SnapshotWarning:
    """
    Lightweight warning record for GUI tables and saved metadata.
    """

    timestamp: str
    snapshot_id: str
    severity: str
    source: str
    code: str
    message: str
    context: dict = field(default_factory=dict)


@dataclass
class SnapshotCorrectorSettings:
    """
    GUI/archiver-facing corrector input for one snapshot.
    """

    hd_corrector_kicks_rad: dict = None
    vd_corrector_kicks_rad: dict = None
    hd_corrector_currents_A: dict = None
    vd_corrector_currents_A: dict = None
    prefer: str = "kicks"
    source: str = "manual"
    metadata: dict = field(default_factory=dict)


@dataclass
class SnapshotOrbitCorrectionConfig:
    """
    Read-only MAD-X orbit-correction request for one snapshot.
    """

    plane: str
    bpm_measurements: object
    correctors: object = None
    fit_knobs: list = None
    label: str = None
    reference_table: str = "bare"
    correction_mode: str = "svd"
    correction_cond: float = 1
    correction_ncorr: int = 0
    monitor_pattern: str = None
    max_fit_kick_rad: float = 0.005
    fit_step: float = 1.0e-4
    fit_calls: int = 50000
    fit_tolerance: float = 1.0e-6
    use_error_bounds: bool = False
    twiss_columns: list = None
    metadata: dict = field(default_factory=dict)


@dataclass
class SnapshotOrbitCorrectionResult:
    """
    Snapshot-level wrapper around orbit_correction.MadxCorrectResult.
    """

    label: str
    plane: str
    result: object


@dataclass
class SnapshotPlotSaveConfig:
    """
    Central plot-save settings for snapshot-generated figures.
    """

    enabled: bool = False
    dpi: int = 200
    root_dir: str = "plots"
    timestamp: str = None
    metadata: dict = field(default_factory=dict)


@dataclass
class SnapshotRunPaths:
    """
    Resolved timestamped output locations for a snapshot or series.
    """

    timestamp: str
    root_dir: str
    snapshot_dir: str = None
    series_dir: str = None
    madx_dir: str = None
    plots_dir: str = None


@dataclass
class SnapshotConfig:
    """
    Configuration for one complete machine snapshot at one cycle time.
    """

    cycle_time_ms: float
    label: str = None
    case: str = "nominal"
    snapshot_id: str = None

    lattice_folder: str = "../Lattice_Files/00_Simplified_Lattice"
    sequence_name: str = "synchrotron"
    output_dir: str = None
    aperture_file: str = "ISIS.aperture"

    main_magnet_mode: str = "rcs_bare"
    requested_qx: float = 4.31
    requested_qy: float = 3.83
    base_qx: float = 4.31
    base_qy: float = 3.83
    tune_method: str = "di_wright"
    harmonics: dict = field(default_factory=dict)

    error_table_paths: list = field(default_factory=list)
    error_table_name: str = "error_table"

    hd_corrector_kicks_rad: dict = None
    vd_corrector_kicks_rad: dict = None
    hd_corrector_currents_A: dict = None
    vd_corrector_currents_A: dict = None
    corrector_prefer: str = "kicks"
    corrector_settings: SnapshotCorrectorSettings = None

    run_envelope: bool = True
    envelope_inputs: EnvelopeInputs = field(default_factory=EnvelopeInputs)

    run_aperture: bool = True
    source_aperture_path: str = str(DEFAULT_SOURCE_APERTURE_CSV)
    aperture_interval: float = 0.1

    branch_configs: list = field(default_factory=list)
    orbit_correction_configs: list = field(default_factory=list)
    plot_save_config: SnapshotPlotSaveConfig = None
    run_timestamp: str = None
    twiss_columns: list = None
    metadata: dict = field(default_factory=dict)

    def resolved_label(self):
        if self.label is not None:
            return str(self.label)
        return f"{self.case}_{float(self.cycle_time_ms):g}ms"

    def resolved_snapshot_id(self, index=None):
        if self.snapshot_id is not None:
            return str(self.snapshot_id)
        if index is None:
            return self.resolved_label()
        return f"snapshot_{int(index):03d}"


@dataclass
class SnapshotResult:
    """
    Evaluated result for one complete machine snapshot.
    """

    config: SnapshotConfig
    snapshot_id: str
    label: str
    case: str
    machine_state: MachineState
    beam_summary: dict
    twiss_df: pd.DataFrame
    madx_summary_df: pd.DataFrame
    tune_summary_df: pd.DataFrame
    orbit_df: pd.DataFrame
    orbit_summary_df: pd.DataFrame
    envelope_result: object = None
    aperture_result: object = None
    source_aperture_result: object = None
    branch_results: list = field(default_factory=list)
    branch_comparison_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    corrector_settings_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    corrector_summary_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_results: list = field(default_factory=list)
    orbit_correction_summary_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_bpm_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_correctors_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_bpm_comparison_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_before_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_after_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    orbit_correction_warnings_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    warnings_df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=WARNING_COLUMNS))
    plot_manifest_df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=PLOT_MANIFEST_COLUMNS))
    run_paths: SnapshotRunPaths = None
    metadata: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def available_tables(self):
        names = [
            "twiss",
            "madx_summary",
            "tune_summary",
            "orbit",
            "orbit_summary",
            "corrector_settings",
            "corrector_summary",
            "warnings",
        ]
        if self.envelope_result is not None:
            names.extend(["envelope", "envelope_summary"])
        if self.aperture_result is not None:
            names.extend(["aperture", "aperture_aligned", "aperture_summary"])
        if self.source_aperture_result is not None:
            names.extend(["source_aperture", "source_aperture_aligned", "source_aperture_summary"])
        if self.branch_results:
            names.extend(["branches", "branch_comparison"])
        if self.orbit_correction_results:
            names.extend(
                [
                    "orbit_correction_summary",
                    "orbit_correction_bpm",
                    "orbit_correction_bpm_comparison",
                    "orbit_correction_correctors",
                    "orbit_correction_before",
                    "orbit_correction_after",
                    "orbit_correction_warnings",
                ]
            )
        if not self.plot_manifest_df.empty:
            names.append("plot_manifest")
        return names

    def table(self, name):
        key = str(name).strip().lower()
        tables = {
            "twiss": self.twiss_df,
            "madx_summary": self.madx_summary_df,
            "tune_summary": self.tune_summary_df,
            "orbit": self.orbit_df,
            "orbit_summary": self.orbit_summary_df,
            "branch_comparison": self.branch_comparison_df,
            "corrector_settings": self.corrector_settings_df,
            "corrector_summary": self.corrector_summary_df,
            "warnings": self.warnings_df,
            "orbit_correction_summary": self.orbit_correction_summary_df,
            "orbit_correction_bpm": self.orbit_correction_bpm_df,
            "orbit_correction_bpm_comparison": self.orbit_correction_bpm_comparison_df,
            "orbit_correction_correctors": self.orbit_correction_correctors_df,
            "orbit_correction_before": self.orbit_correction_before_df,
            "orbit_correction_after": self.orbit_correction_after_df,
            "orbit_correction_warnings": self.orbit_correction_warnings_df,
            "plot_manifest": self.plot_manifest_df,
        }
        if self.envelope_result is not None:
            tables["envelope"] = self.envelope_result.envelope_df
            tables["envelope_summary"] = self.envelope_result.summary_df
        if self.aperture_result is not None:
            tables["aperture"] = self.aperture_result.aperture_df
            tables["aperture_aligned"] = self.aperture_result.aligned_df
            tables["aperture_summary"] = self.aperture_result.summary_df
        if self.source_aperture_result is not None:
            tables["source_aperture"] = self.source_aperture_result.aperture_df
            tables["source_aperture_aligned"] = self.source_aperture_result.aligned_df
            tables["source_aperture_summary"] = self.source_aperture_result.summary_df
        if key == "branches":
            rows = []
            for result in self.branch_results:
                rows.append(
                    {
                        "snapshot_id": self.snapshot_id,
                        "branch": result.name,
                        "twiss_rows": len(result.twiss_df),
                        "orbit_rows": len(result.orbit_df),
                        **dict(result.metadata or {}),
                    }
                )
            return pd.DataFrame(rows)
        if key not in tables:
            raise KeyError(f"Unknown snapshot table {name!r}. Available tables: {self.available_tables()}")
        return tables[key].copy()

    def to_summary_dataframe(self):
        row = {
            "snapshot_id": self.snapshot_id,
            "label": self.label,
            "case": self.case,
            "cycle_time_ms": float(self.config.cycle_time_ms),
            "set_qx": float(self.config.requested_qx),
            "set_qy": float(self.config.requested_qy),
            "twiss_rows": len(self.twiss_df),
            "orbit_rows": len(self.orbit_df),
            "has_envelope": self.envelope_result is not None,
            "has_aperture": self.aperture_result is not None,
            "n_branches": len(self.branch_results),
            "n_orbit_corrections": len(self.orbit_correction_results),
            "n_warnings": len(self.warnings),
        }
        if not self.tune_summary_df.empty:
            tune_row = self.tune_summary_df.iloc[0].to_dict()
            row.update(
                {
                    "predicted_qx": tune_row.get("predicted_qx", tune_row.get("actual_qx")),
                    "predicted_qy": tune_row.get("predicted_qy", tune_row.get("actual_qy")),
                    "dqx": tune_row.get("dqx"),
                    "dqy": tune_row.get("dqy"),
                    "madx_dqx_dpt": tune_row.get("madx_dqx_dpt"),
                    "madx_dqy_dpt": tune_row.get("madx_dqy_dpt"),
                    "lorentz_beta": tune_row.get("lorentz_beta"),
                    "iqtf_A": tune_row.get("iqtf_A"),
                    "iqtd_A": tune_row.get("iqtd_A"),
                    "kqtf": tune_row.get("kqtf"),
                    "kqtd": tune_row.get("kqtd"),
                }
            )
        return pd.DataFrame([row])


@dataclass
class SnapshotSeriesConfig:
    """
    Configuration for an ordered set of snapshots.
    """

    snapshots: list
    label: str = "snapshot_series"
    output_dir: str = None
    run_timestamp: str = None
    xlims: tuple = (4.0, 4.5)
    ylims: tuple = (3.5, 4.0)
    orders: tuple = (1, 2, 3, 4)
    periodicity: int = 10
    metadata: dict = field(default_factory=dict)


@dataclass
class SnapshotSeriesResult:
    """
    Evaluated result for multiple snapshots.
    """

    config: SnapshotSeriesConfig
    snapshots: list
    summary_df: pd.DataFrame
    tune_programme_df: pd.DataFrame
    working_points_df: pd.DataFrame
    resonance_lines_df: pd.DataFrame
    resonance_proximity_df: pd.DataFrame
    warnings_df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=WARNING_COLUMNS))
    plot_manifest_df: pd.DataFrame = field(default_factory=lambda: pd.DataFrame(columns=PLOT_MANIFEST_COLUMNS))
    run_paths: SnapshotRunPaths = None
    metadata: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def available_tables(self):
        return [
            "summary",
            "tune_programme",
            "working_points",
            "resonance_lines",
            "resonance_proximity",
            "warnings",
            "plot_manifest",
        ]

    def table(self, name):
        key = str(name).strip().lower()
        tables = {
            "summary": self.summary_df,
            "tune_programme": self.tune_programme_df,
            "working_points": self.working_points_df,
            "resonance_lines": self.resonance_lines_df,
            "resonance_proximity": self.resonance_proximity_df,
            "warnings": self.warnings_df,
            "plot_manifest": self.plot_manifest_df,
        }
        if key not in tables:
            raise KeyError(f"Unknown snapshot series table {name!r}. Available tables: {self.available_tables()}")
        return tables[key].copy()

    def to_summary_dataframe(self):
        return self.summary_df.copy()


def copy_snapshot_config(config, **overrides):
    """
    Return a deep-copied SnapshotConfig with explicit overrides.
    """

    if not isinstance(config, SnapshotConfig):
        raise TypeError("config must be a SnapshotConfig.")
    copied = copy.deepcopy(config)
    for key, value in overrides.items():
        if not hasattr(copied, key):
            raise AttributeError(f"SnapshotConfig has no field {key!r}.")
        setattr(copied, key, value)
    return copied


def copy_snapshot_result_config(result, **overrides):
    """
    Return a copied config derived from an evaluated SnapshotResult.
    """

    if not isinstance(result, SnapshotResult):
        raise TypeError("result must be a SnapshotResult.")
    return copy_snapshot_config(result.config, **overrides)


def _now_timestamp():
    return datetime.now().strftime(TIMESTAMP_FORMAT)


def _safe_name(name):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_")
    return safe or "snapshot"


def _normalise_timestamp(timestamp=None):
    return str(timestamp or _now_timestamp())


def _require_output_dir(output_dir, owner="SnapshotConfig"):
    if output_dir is None or str(output_dir).strip() == "":
        raise ValueError(f"{owner}.output_dir must be explicitly supplied.")
    return Path(output_dir)


def _display_error_table_paths(paths):
    return "+".join(Path(path).name for path in (paths or [])) or None


def _warning_record(snapshot_id, source, message, severity="warning", code=None, context=None):
    return SnapshotWarning(
        timestamp=datetime.now().isoformat(timespec="seconds"),
        snapshot_id=str(snapshot_id),
        severity=str(severity),
        source=str(source),
        code=str(code or source),
        message=str(message),
        context={} if context is None else dict(context),
    )


def _warnings_to_dataframe(warnings):
    rows = []
    for warning in warnings or []:
        if isinstance(warning, SnapshotWarning):
            rows.append(asdict(warning))
        elif isinstance(warning, dict):
            row = {column: warning.get(column) for column in WARNING_COLUMNS}
            row["context"] = row["context"] or {}
            rows.append(row)
        else:
            rows.append(
                asdict(
                    _warning_record(
                        snapshot_id="",
                        source="snapshot",
                        message=str(warning),
                    )
                )
            )
    return pd.DataFrame(rows, columns=WARNING_COLUMNS)


def create_snapshot_run_paths(config, snapshot_id, run_dir=None, timestamp=None):
    timestamp = _normalise_timestamp(timestamp or config.run_timestamp)
    if run_dir is None:
        root = _require_output_dir(config.output_dir)
        snapshot_dir = root / f"{timestamp}_{_safe_name(snapshot_id)}"
    else:
        snapshot_dir = Path(run_dir)
        root = snapshot_dir.parent
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    madx_dir = snapshot_dir / "madx"
    plots_dir = snapshot_dir / "plots"
    madx_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    return SnapshotRunPaths(
        timestamp=timestamp,
        root_dir=str(root),
        snapshot_dir=str(snapshot_dir),
        madx_dir=str(madx_dir),
        plots_dir=str(plots_dir),
    )


def create_series_run_paths(config, timestamp=None):
    timestamp = _normalise_timestamp(timestamp or config.run_timestamp)
    output_dir = config.output_dir
    if output_dir is None and config.snapshots:
        output_dir = config.snapshots[0].output_dir
    root = _require_output_dir(output_dir, owner="SnapshotSeriesConfig")
    series_dir = root / f"{timestamp}_{_safe_name(config.label)}"
    series_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = series_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    return SnapshotRunPaths(
        timestamp=timestamp,
        root_dir=str(root),
        series_dir=str(series_dir),
        plots_dir=str(plots_dir),
    )


def _snapshot_output_dir(config, snapshot_id):
    return Path(create_snapshot_run_paths(config, snapshot_id).snapshot_dir)


def _validate_snapshot_config(config):
    if not isinstance(config, SnapshotConfig):
        raise TypeError("config must be a SnapshotConfig.")
    _require_output_dir(config.output_dir)
    if config.run_aperture and not config.run_envelope:
        raise ValueError("run_aperture=True requires run_envelope=True.")
    if config.corrector_prefer not in ("kicks", "currents"):
        raise ValueError("corrector_prefer must be 'kicks' or 'currents'.")
    if config.corrector_settings is not None and not isinstance(config.corrector_settings, SnapshotCorrectorSettings):
        raise TypeError("corrector_settings must be a SnapshotCorrectorSettings object.")
    if config.corrector_settings is not None and config.corrector_settings.prefer not in ("kicks", "currents"):
        raise ValueError("corrector_settings.prefer must be 'kicks' or 'currents'.")
    if config.plot_save_config is not None and not isinstance(config.plot_save_config, SnapshotPlotSaveConfig):
        raise TypeError("plot_save_config must be a SnapshotPlotSaveConfig object.")
    if config.lattice_folder is None or not Path(config.lattice_folder).exists():
        raise FileNotFoundError(f"Missing lattice folder: {config.lattice_folder}")
    for path in config.error_table_paths or []:
        if not Path(path).exists():
            raise FileNotFoundError(f"Missing error table: {path}")
    if config.run_aperture and config.source_aperture_path is not None and not Path(config.source_aperture_path).exists():
        raise FileNotFoundError(f"Missing source aperture CSV: {config.source_aperture_path}")
    for branch in config.branch_configs or []:
        if not isinstance(branch, (OrbitBranchConfig, dict)):
            raise TypeError("branch_configs entries must be OrbitBranchConfig objects or dictionaries.")
    for correction in config.orbit_correction_configs or []:
        if not isinstance(correction, (SnapshotOrbitCorrectionConfig, dict)):
            raise TypeError("orbit_correction_configs entries must be SnapshotOrbitCorrectionConfig objects or dictionaries.")
    return True


def corrector_settings_from_manual(
    hd_corrector_kicks_rad=None,
    vd_corrector_kicks_rad=None,
    hd_corrector_currents_A=None,
    vd_corrector_currents_A=None,
    prefer="kicks",
    source="manual",
    metadata=None,
):
    return SnapshotCorrectorSettings(
        hd_corrector_kicks_rad=hd_corrector_kicks_rad,
        vd_corrector_kicks_rad=vd_corrector_kicks_rad,
        hd_corrector_currents_A=hd_corrector_currents_A,
        vd_corrector_currents_A=vd_corrector_currents_A,
        prefer=prefer,
        source=source,
        metadata={} if metadata is None else dict(metadata),
    )


def corrector_settings_from_dataframe(dataframe, cycle_time_ms=None, prefer="currents", source="epics_archiver"):
    """
    Build SnapshotCorrectorSettings from a GUI/EPICS-style table.

    Expected columns are corrector plus either current_A or kick_rad. Optional
    plane/family columns may identify HD/VD; otherwise names starting with hd
    or vd are used.
    """

    df = pd.DataFrame(dataframe).copy()
    if cycle_time_ms is not None and "cycle_time_ms" in df.columns:
        times = pd.to_numeric(df["cycle_time_ms"], errors="coerce")
        df = df[times == float(cycle_time_ms)].copy()
    if df.empty:
        raise ValueError("No corrector rows matched the requested cycle time.")
    if "corrector" not in df.columns:
        for candidate in ("name", "device", "pv", "corrector_name"):
            if candidate in df.columns:
                df["corrector"] = df[candidate]
                break
    if "corrector" not in df.columns:
        raise ValueError("Corrector table must contain a corrector/name/device column.")

    current_column = next((c for c in ("current_A", "current", "value_A", "value") if c in df.columns), None)
    kick_column = next((c for c in ("kick_rad", "kick", "value_rad") if c in df.columns), None)
    if current_column is None and kick_column is None:
        raise ValueError("Corrector table must contain current_A/current or kick_rad/kick values.")

    hd_currents, vd_currents, hd_kicks, vd_kicks = {}, {}, {}, {}
    for _, row in df.iterrows():
        name = str(row["corrector"])
        plane = str(row.get("plane", row.get("family", ""))).upper()
        is_hd = plane.startswith("H") or name.lower().startswith("hd")
        is_vd = plane.startswith("V") or name.lower().startswith("vd")
        if not (is_hd or is_vd):
            raise ValueError(f"Cannot infer HD/VD plane for corrector {name!r}.")
        if current_column is not None:
            target = hd_currents if is_hd else vd_currents
            target[name] = float(row[current_column])
        if kick_column is not None:
            target = hd_kicks if is_hd else vd_kicks
            target[name] = float(row[kick_column])

    return SnapshotCorrectorSettings(
        hd_corrector_kicks_rad=hd_kicks or None,
        vd_corrector_kicks_rad=vd_kicks or None,
        hd_corrector_currents_A=hd_currents or None,
        vd_corrector_currents_A=vd_currents or None,
        prefer=prefer,
        source=source,
        metadata={"cycle_time_ms": cycle_time_ms},
    )


def _resolved_corrector_settings(config):
    if config.corrector_settings is not None:
        return config.corrector_settings
    return SnapshotCorrectorSettings(
        hd_corrector_kicks_rad=config.hd_corrector_kicks_rad,
        vd_corrector_kicks_rad=config.vd_corrector_kicks_rad,
        hd_corrector_currents_A=config.hd_corrector_currents_A,
        vd_corrector_currents_A=config.vd_corrector_currents_A,
        prefer=config.corrector_prefer,
        source="snapshot_config",
    )


def _corrector_tables(machine_state, source=None):
    rows = []
    for family, kicks, currents in (
        ("HD", machine_state.hd_corrector_kicks_rad, machine_state.hd_corrector_currents_A),
        ("VD", machine_state.vd_corrector_kicks_rad, machine_state.vd_corrector_currents_A),
    ):
        for name in kicks:
            rows.append(
                {
                    "family": family,
                    "corrector": name,
                    "kick_rad": float(kicks[name]),
                    "kick_mrad": 1.0e3 * float(kicks[name]),
                    "current_A": float(currents[name]),
                    "source": source,
                }
            )
    settings_df = pd.DataFrame(rows)
    if settings_df.empty:
        summary_df = pd.DataFrame()
    else:
        summary_df = (
            settings_df.groupby("family")
            .agg(
                n_correctors=("corrector", "count"),
                max_abs_current_A=("current_A", lambda s: float(s.abs().max())),
                max_abs_kick_mrad=("kick_mrad", lambda s: float(s.abs().max())),
            )
            .reset_index()
        )
    return settings_df, summary_df


def _build_machine_state(config, beam_state):
    harmonics = normalise_harmonic_inputs(config.harmonics)
    error_path = _display_error_table_paths(config.error_table_paths)
    corrector_settings = _resolved_corrector_settings(config)
    return MachineState.from_defaults(
        beam_state=beam_state,
        main_magnet_mode=config.main_magnet_mode,
        requested_qx=config.requested_qx,
        requested_qy=config.requested_qy,
        base_qx=config.base_qx,
        base_qy=config.base_qy,
        tune_method=config.tune_method,
        harmonic_tunes=harmonics,
        error_table_path=error_path,
        hd_corrector_kicks_rad=corrector_settings.hd_corrector_kicks_rad,
        vd_corrector_kicks_rad=corrector_settings.vd_corrector_kicks_rad,
        hd_corrector_currents_A=corrector_settings.hd_corrector_currents_A,
        vd_corrector_currents_A=corrector_settings.vd_corrector_currents_A,
        corrector_prefer=corrector_settings.prefer,
        metadata={
            **dict(config.metadata or {}),
            "snapshot_label": config.resolved_label(),
            "snapshot_case": config.case,
            "corrector_source": corrector_settings.source,
            "corrector_metadata": dict(corrector_settings.metadata or {}),
        },
    )


def _summary_df_from_tune(config, summary_df, machine_state):
    summary = extract_tune_summary(summary_df, source="madx_summary")
    beta_rel = float(machine_state.beam_state.beta)
    madx_dqx_dpt = summary["dqx"]
    madx_dqy_dpt = summary["dqy"]
    row = {
        "cycle_time_ms": float(config.cycle_time_ms),
        "set_qx": float(config.requested_qx),
        "set_qy": float(config.requested_qy),
        "predicted_qx": summary["qx"],
        "predicted_qy": summary["qy"],
        "dqx": None if madx_dqx_dpt is None else beta_rel * madx_dqx_dpt,
        "dqy": None if madx_dqy_dpt is None else beta_rel * madx_dqy_dpt,
        "madx_dqx_dpt": madx_dqx_dpt,
        "madx_dqy_dpt": madx_dqy_dpt,
        "lorentz_beta": beta_rel,
        "iqtf_A": machine_state.iqtf_A,
        "iqtd_A": machine_state.iqtd_A,
        "kqtf": machine_state.kqtf,
        "kqtd": machine_state.kqtd,
        "source": summary["source"],
    }
    return pd.DataFrame([row])


def _apply_error_tables(model, config):
    applied = []
    for index, path in enumerate(config.error_table_paths or []):
        table_name = config.error_table_name
        if len(config.error_table_paths) > 1:
            table_name = f"{config.error_table_name}_{index}"
        applied.append(model.apply_error_table(path, table_name=table_name))
    return applied


def _build_branch_results(config, machine_state, output_dir):
    results = []
    for index, branch in enumerate(config.branch_configs or []):
        if isinstance(branch, OrbitBranchConfig):
            branch_config = replace(
                branch,
                machine_state=branch.machine_state or machine_state,
                lattice_folder=branch.lattice_folder or config.lattice_folder,
                sequence_name=branch.sequence_name or config.sequence_name,
                output_dir=branch.output_dir or str(output_dir / f"branch_{index:02d}_{branch.name}"),
            )
        elif isinstance(branch, dict):
            branch_config = OrbitBranchConfig(
                machine_state=branch.get("machine_state", machine_state),
                lattice_folder=branch.get("lattice_folder", config.lattice_folder),
                sequence_name=branch.get("sequence_name", config.sequence_name),
                output_dir=branch.get("output_dir", str(output_dir / f"branch_{index:02d}_{branch.get('name', 'branch')}")),
                **{key: value for key, value in branch.items() if key not in {"machine_state", "lattice_folder", "sequence_name", "output_dir"}},
            )
        else:
            raise TypeError("branch_configs entries must be OrbitBranchConfig objects or dictionaries.")
        results.append(OrbitBranch(branch_config).run())
    return results


def _branch_comparison(reference_orbit, branch_results, snapshot_id):
    rows = []
    for result in branch_results:
        comparison = compare_orbit_results(reference_orbit, result.orbit_df)
        comparison.insert(0, "branch", result.name)
        comparison.insert(0, "snapshot_id", snapshot_id)
        rows.append(comparison)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def _normalise_orbit_correction_config(config):
    if isinstance(config, SnapshotOrbitCorrectionConfig):
        return config
    if isinstance(config, dict):
        return SnapshotOrbitCorrectionConfig(**config)
    raise TypeError("orbit correction entries must be SnapshotOrbitCorrectionConfig objects or dictionaries.")


def _flatten_summary(prefix, summary):
    row = {}
    for key, value in dict(summary or {}).items():
        if isinstance(value, dict):
            for nested_key, nested_value in value.items():
                row[f"{key}_{nested_key}"] = nested_value
        else:
            row[key] = value
    if prefix:
        return {f"{prefix}_{key}": value for key, value in row.items()}
    return row


def _bpm_comparison_table(result, snapshot_id, correction_label):
    before = orbit_residuals(result.measured_twiss_df, result.bpm_measurements, result.plane)
    after = orbit_residuals(result.corrected_twiss_df, result.bpm_measurements, result.plane)
    compare = before.rename(
        columns={
            "model_mm": "before_model_mm",
            "residual_mm": "before_residual_mm",
        }
    )
    after = after.rename(
        columns={
            "model_mm": "after_model_mm",
            "residual_mm": "after_residual_mm",
        }
    )
    keep_after = ["bpm", "plane", "after_model_mm", "after_residual_mm"]
    compare = compare.merge(after[keep_after], on=["bpm", "plane"], how="left")
    if "closed_orbit_mm_err" in result.bpm_measurements.columns:
        errors = result.bpm_measurements[["bpm", "plane", "closed_orbit_mm_err"]].copy()
        errors = errors.rename(columns={"closed_orbit_mm_err": "measurement_error_mm"})
        compare = compare.merge(errors, on=["bpm", "plane"], how="left")
    else:
        compare["measurement_error_mm"] = pd.NA
    compare["residual_change_mm"] = compare["after_residual_mm"] - compare["before_residual_mm"]
    compare.insert(0, "correction_label", correction_label)
    compare.insert(0, "snapshot_id", snapshot_id)
    ordered = [
        "snapshot_id",
        "correction_label",
        "plane",
        "bpm",
        "enabled",
        "s",
        "matched_name",
        "measurement_mm",
        "measurement_error_mm",
        "before_model_mm",
        "before_residual_mm",
        "after_model_mm",
        "after_residual_mm",
        "residual_change_mm",
    ]
    return compare[[column for column in ordered if column in compare.columns]]


def _build_orbit_correction_results(config, machine_state, run_paths, snapshot_id):
    results = []
    warning_records = []
    summary_rows = []
    bpm_rows = []
    bpm_comparison_rows = []
    corrector_rows = []
    before_rows = []
    after_rows = []

    for index, raw_correction in enumerate(config.orbit_correction_configs or []):
        correction = _normalise_orbit_correction_config(raw_correction)
        label = correction.label or f"orbit_correction_{index:02d}_{correction.plane}"
        try:
            result = correct_bpm_measurements_with_madx_correct(
                plane=correction.plane,
                lattice_folder=config.lattice_folder,
                bpm_measurements=correction.bpm_measurements,
                machine_state=machine_state,
                correctors=correction.correctors,
                fit_knobs=correction.fit_knobs,
                sequence_name=config.sequence_name,
                output_dir=Path(run_paths.snapshot_dir) / "orbit_correction" / _safe_name(label),
                reference_table=correction.reference_table,
                correction_mode=correction.correction_mode,
                correction_cond=correction.correction_cond,
                correction_ncorr=correction.correction_ncorr,
                monitor_pattern=correction.monitor_pattern,
                max_fit_kick_rad=correction.max_fit_kick_rad,
                fit_step=correction.fit_step,
                fit_calls=correction.fit_calls,
                fit_tolerance=correction.fit_tolerance,
                use_error_bounds=correction.use_error_bounds,
                twiss_columns=correction.twiss_columns,
                metadata={
                    **dict(correction.metadata or {}),
                    "snapshot_id": snapshot_id,
                    "correction_label": label,
                },
            )
        except Exception as exc:
            warning_records.append(
                _warning_record(
                    snapshot_id=snapshot_id,
                    source="orbit_correction",
                    code="orbit_correction_failed",
                    severity="error",
                    message=str(exc),
                    context={"label": label, "plane": correction.plane},
                )
            )
            continue

        results.append(SnapshotOrbitCorrectionResult(label=label, plane=result.plane, result=result))

        summary_row = {
            "snapshot_id": snapshot_id,
            "correction_label": label,
            **_flatten_summary(None, result.summary),
        }
        summary_rows.append(summary_row)

        bpm = result.bpm_measurements.copy()
        bpm.insert(0, "correction_label", label)
        bpm.insert(0, "snapshot_id", snapshot_id)
        bpm_rows.append(bpm)
        bpm_comparison_rows.append(_bpm_comparison_table(result, snapshot_id, label))

        correctors = result.correctors.copy()
        correctors.insert(0, "correction_label", label)
        correctors.insert(0, "snapshot_id", snapshot_id)
        corrector_rows.append(correctors)

        before_rows.append(
            {
                "snapshot_id": snapshot_id,
                "correction_label": label,
                "plane": result.plane,
                **result.monitor_summary_before,
            }
        )
        after_rows.append(
            {
                "snapshot_id": snapshot_id,
                "correction_label": label,
                "plane": result.plane,
                **result.monitor_summary_after,
            }
        )
        for warning in result.warnings or []:
            warning_records.append(
                _warning_record(
                    snapshot_id=snapshot_id,
                    source="orbit_correction",
                    code="orbit_correction_warning",
                    message=warning,
                    context={"label": label, "plane": result.plane},
                )
            )

    return {
        "results": results,
        "summary_df": pd.DataFrame(summary_rows),
        "bpm_df": pd.concat(bpm_rows, ignore_index=True) if bpm_rows else pd.DataFrame(),
        "bpm_comparison_df": pd.concat(bpm_comparison_rows, ignore_index=True) if bpm_comparison_rows else pd.DataFrame(),
        "correctors_df": pd.concat(corrector_rows, ignore_index=True) if corrector_rows else pd.DataFrame(),
        "before_df": pd.DataFrame(before_rows),
        "after_df": pd.DataFrame(after_rows),
        "warnings_df": _warnings_to_dataframe(warning_records),
        "warnings": warning_records,
    }


def _figure_from_plot_object(plot_object):
    if hasattr(plot_object, "figure"):
        return plot_object.figure
    if isinstance(plot_object, (list, tuple)) and plot_object:
        return _figure_from_plot_object(plot_object[0])
    if hasattr(plot_object, "savefig"):
        return plot_object
    raise TypeError("plot_object must be a Matplotlib figure, axis or axes collection.")


def _plot_title(plot_object):
    if hasattr(plot_object, "get_title"):
        return plot_object.get_title()
    if isinstance(plot_object, (list, tuple)) and plot_object:
        titles = [_plot_title(item) for item in plot_object if hasattr(item, "get_title")]
        return "; ".join(title for title in titles if title)
    if hasattr(plot_object, "axes") and plot_object.axes:
        titles = [axis.get_title() for axis in plot_object.axes if axis.get_title()]
        return "; ".join(titles)
    return ""


def _plot_completeness(plot_object):
    figure = _figure_from_plot_object(plot_object)
    axes = list(getattr(figure, "axes", []))
    if not axes:
        return {"has_title": False, "has_axis_labels": False, "has_legend": False}
    titles = [axis.get_title() for axis in axes]
    labels = [(axis.get_xlabel(), axis.get_ylabel()) for axis in axes]
    legends = [axis.get_legend() for axis in axes]
    return {
        "has_title": any(bool(title) for title in titles),
        "has_axis_labels": all(bool(xlabel) and bool(ylabel) for xlabel, ylabel in labels),
        "has_legend": any(legend is not None for legend in legends),
    }


def save_snapshot_plot(
    snapshot,
    plot_object,
    plot_name,
    aspect,
    dpi=None,
    timestamp=None,
):
    """
    Save one Matplotlib figure under the snapshot's timestamped plot tree.
    """

    if not isinstance(snapshot, SnapshotResult):
        raise TypeError("snapshot must be a SnapshotResult.")
    if snapshot.run_paths is None:
        raise ValueError("snapshot.run_paths is required to save plots.")

    save_config = snapshot.config.plot_save_config or SnapshotPlotSaveConfig()
    timestamp = _normalise_timestamp(timestamp or save_config.timestamp or snapshot.run_paths.timestamp)
    dpi = int(dpi or save_config.dpi or 200)
    aspect = _safe_name(aspect)
    plot_name = _safe_name(plot_name)
    root_dir = Path(snapshot.run_paths.snapshot_dir) / str(save_config.root_dir or "plots")
    output_dir = root_dir / aspect
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{timestamp}_{_safe_name(snapshot.snapshot_id)}_{plot_name}.png"

    figure = _figure_from_plot_object(plot_object)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    completeness = _plot_completeness(plot_object)
    row = {
        "timestamp": timestamp,
        "snapshot_id": snapshot.snapshot_id,
        "aspect": aspect,
        "plot_name": plot_name,
        "path": str(path),
        "dpi": dpi,
        "title": _plot_title(plot_object),
        **completeness,
    }
    snapshot.plot_manifest_df = pd.concat(
        [snapshot.plot_manifest_df, pd.DataFrame([row], columns=PLOT_MANIFEST_COLUMNS)],
        ignore_index=True,
    )
    return path


def save_snapshot_plots(snapshot, plots, dpi=None, timestamp=None):
    """
    Save a sequence of plot definitions and return the updated manifest table.
    """

    for item in plots:
        save_snapshot_plot(
            snapshot=snapshot,
            plot_object=item["plot_object"],
            plot_name=item["plot_name"],
            aspect=item.get("aspect", "misc"),
            dpi=dpi,
            timestamp=timestamp,
        )
    return snapshot.table("plot_manifest")


def build_machine_snapshot(config, snapshot_index=None, run_dir=None, run_timestamp=None):
    """
    Build one complete machine snapshot.
    """

    _validate_snapshot_config(config)

    snapshot_id = config.resolved_snapshot_id(snapshot_index)
    label = config.resolved_label()
    run_paths = create_snapshot_run_paths(
        config,
        snapshot_id=snapshot_id,
        run_dir=run_dir,
        timestamp=run_timestamp,
    )
    output_dir = Path(run_paths.snapshot_dir)

    beam_state = RCSRamp().state_at(config.cycle_time_ms)
    machine_state = _build_machine_state(config, beam_state)
    corrector_settings = _resolved_corrector_settings(config)
    corrector_settings_df, corrector_summary_df = _corrector_tables(
        machine_state,
        source=corrector_settings.source,
    )

    model = MadxModel(
        lattice_folder=config.lattice_folder,
        sequence_name=config.sequence_name,
        aperture_file=config.aperture_file if config.run_aperture else None,
        output_dir=run_paths.madx_dir,
    )
    model.load_lattice(use_sequence=False)
    model.apply_machine_state(machine_state)
    model.use_sequence()
    applied_errors = _apply_error_tables(model, config)

    twiss_df = model.run_twiss(columns=config.twiss_columns)
    madx_summary_df = model.get_summary_df(refresh=True)
    tune_summary_df = _summary_df_from_tune(config, madx_summary_df, machine_state)
    orbit_df = extract_orbit_df(twiss_df)
    orbit_summary_df = pd.DataFrame([summarise_orbit_df(orbit_df)])

    warning_records = []
    envelope_result = None
    if config.run_envelope:
        envelope_result = evaluate_envelope_from_twiss(
            twiss_df,
            inputs=config.envelope_inputs,
            beam_state=machine_state.beam_state,
            source="snapshot",
            metadata={"snapshot_id": snapshot_id, "label": label},
        )
        warning_records.extend(
            _warning_record(snapshot_id, "envelope", warning, code="envelope_warning")
            for warning in (envelope_result.warnings or [])
        )

    aperture_result = None
    source_aperture_result = None
    if config.run_aperture:
        raw_aperture_df = model.run_aperture(interval=config.aperture_interval)
        normalised_aperture_df = normalise_madx_aperture_table(raw_aperture_df)
        aperture_result = evaluate_aperture_margins(
            normalised_aperture_df,
            envelope_result,
            label=f"{label} MAD-X aperture",
            source="madx_aperture",
        )
        warning_records.extend(
            _warning_record(snapshot_id, "aperture", warning, code="aperture_warning")
            for warning in (aperture_result.warnings or [])
        )

        source_aperture_df = read_source_aperture_csv(config.source_aperture_path)
        source_aperture_result = evaluate_aperture_margins(
            source_aperture_df,
            envelope_result,
            label=f"{label} source aperture",
            source="source_aperture",
        )
        warning_records.extend(
            _warning_record(snapshot_id, "source_aperture", warning, code="source_aperture_warning")
            for warning in (source_aperture_result.warnings or [])
        )

    branch_results = _build_branch_results(config, machine_state, output_dir)
    branch_comparison_df = _branch_comparison(orbit_df, branch_results, snapshot_id)

    correction_bundle = _build_orbit_correction_results(config, machine_state, run_paths, snapshot_id)
    warning_records.extend(correction_bundle["warnings"])
    warnings_df = _warnings_to_dataframe(warning_records)

    metadata = {
        "snapshot_id": snapshot_id,
        "label": label,
        "case": config.case,
        "output_dir": str(output_dir),
        "run_paths": asdict(run_paths),
        "applied_error_tables": list(applied_errors),
        "error_table_display": _display_error_table_paths(config.error_table_paths),
        "machine_state": machine_state.summary_dict(),
        "madx": model.get_metadata(),
        **dict(config.metadata or {}),
    }

    return SnapshotResult(
        config=copy.deepcopy(config),
        snapshot_id=snapshot_id,
        label=label,
        case=config.case,
        machine_state=machine_state,
        beam_summary=machine_state.beam_summary_dict(),
        twiss_df=twiss_df,
        madx_summary_df=madx_summary_df,
        tune_summary_df=tune_summary_df,
        orbit_df=orbit_df,
        orbit_summary_df=orbit_summary_df,
        envelope_result=envelope_result,
        aperture_result=aperture_result,
        source_aperture_result=source_aperture_result,
        branch_results=branch_results,
        branch_comparison_df=branch_comparison_df,
        corrector_settings_df=corrector_settings_df,
        corrector_summary_df=corrector_summary_df,
        orbit_correction_results=correction_bundle["results"],
        orbit_correction_summary_df=correction_bundle["summary_df"],
        orbit_correction_bpm_df=correction_bundle["bpm_df"],
        orbit_correction_bpm_comparison_df=correction_bundle["bpm_comparison_df"],
        orbit_correction_correctors_df=correction_bundle["correctors_df"],
        orbit_correction_before_df=correction_bundle["before_df"],
        orbit_correction_after_df=correction_bundle["after_df"],
        orbit_correction_warnings_df=correction_bundle["warnings_df"],
        warnings_df=warnings_df,
        run_paths=run_paths,
        metadata=metadata,
        warnings=warning_records,
    )


def _series_summary_rows(snapshot_results):
    rows = []
    for result in snapshot_results:
        row = result.to_summary_dataframe().iloc[0].to_dict()
        row["warnings"] = list(result.warnings)
        rows.append(row)
    return pd.DataFrame(rows)


def _series_warnings(snapshot_results):
    frames = [result.table("warnings") for result in snapshot_results if "warnings" in result.available_tables()]
    if not frames:
        return pd.DataFrame(columns=WARNING_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def _series_plot_manifest(snapshot_results):
    frames = [result.table("plot_manifest") for result in snapshot_results if not result.table("plot_manifest").empty]
    if not frames:
        return pd.DataFrame(columns=PLOT_MANIFEST_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def _series_tune_programme(summary_df):
    rows = []
    for _, row in summary_df.iterrows():
        rows.append(
            {
                "snapshot_id": row["snapshot_id"],
                "label": row["label"],
                "case": row["case"],
                "cycle_time_ms": row["cycle_time_ms"],
                "set_qx": row["set_qx"],
                "set_qy": row["set_qy"],
                "predicted_qx": row.get("predicted_qx"),
                "predicted_qy": row.get("predicted_qy"),
                "dqx": row.get("dqx"),
                "dqy": row.get("dqy"),
                "madx_dqx_dpt": row.get("madx_dqx_dpt"),
                "madx_dqy_dpt": row.get("madx_dqy_dpt"),
                "lorentz_beta": row.get("lorentz_beta"),
                "iqtf_A": row.get("iqtf_A"),
                "iqtd_A": row.get("iqtd_A"),
                "kqtf": row.get("kqtf"),
                "kqtd": row.get("kqtd"),
                "source": "snapshot_series",
            }
        )
    return build_tune_programme_table(pd.DataFrame(rows))


def _attach_series_identity(table, summary_df, include_index=False):
    out = table.copy()
    identities = summary_df[["snapshot_id", "label", "case"]].reset_index(drop=True)
    if len(out) != len(identities):
        return out
    for column in reversed(list(identities.columns)):
        if column in out.columns:
            out[column] = identities[column].to_numpy()
        else:
            out.insert(0, column, identities[column].to_numpy())
    if include_index and "snapshot_index" not in out.columns:
        out.insert(0, "snapshot_index", range(len(out)))
    return out


def _drop_all_missing_columns(table):
    out = table.copy()
    empty_columns = [column for column in out.columns if out[column].isna().all()]
    if empty_columns:
        out = out.drop(columns=empty_columns)
    return out


def build_snapshot_series(config):
    """
    Build an ordered set of machine snapshots and resonance-diagram inputs.
    """

    if not isinstance(config, SnapshotSeriesConfig):
        raise TypeError("config must be a SnapshotSeriesConfig.")
    if not config.snapshots:
        raise ValueError("SnapshotSeriesConfig.snapshots must not be empty.")
    run_paths = create_series_run_paths(config)
    snapshots = [
        build_machine_snapshot(
            copy_snapshot_config(
                snapshot_config,
                output_dir=str(Path(run_paths.series_dir) / "snapshots"),
            ),
            snapshot_index=index,
            run_dir=Path(run_paths.series_dir) / "snapshots" / _safe_name(snapshot_config.resolved_snapshot_id(index)),
            run_timestamp=run_paths.timestamp,
        )
        for index, snapshot_config in enumerate(config.snapshots)
    ]
    summary_df = _series_summary_rows(snapshots)
    tune_programme_df = _series_tune_programme(summary_df)
    diagram = make_tune_diagram_inputs(
        tune_programme_df,
        xlims=config.xlims,
        ylims=config.ylims,
        orders=config.orders,
        periodicity=config.periodicity,
    )
    tune_programme = _drop_all_missing_columns(
        _attach_series_identity(diagram["programme"], summary_df, include_index=True)
    )
    working_points = _drop_all_missing_columns(
        _attach_series_identity(diagram["working_points"], summary_df, include_index=True)
    )
    proximity = diagram["resonance_proximity"].copy()
    if not proximity.empty and "point_index" in proximity.columns:
        identity = working_points[
            ["snapshot_index", "snapshot_id", "label", "case"]
        ].reset_index(drop=True)
        proximity = proximity.merge(
            identity,
            left_on="point_index",
            right_index=True,
            how="left",
        )
    warnings = []
    for result in snapshots:
        warnings.extend(result.warnings)
    warnings_df = _series_warnings(snapshots)
    plot_manifest_df = _series_plot_manifest(snapshots)
    return SnapshotSeriesResult(
        config=copy.deepcopy(config),
        snapshots=snapshots,
        summary_df=summary_df,
        tune_programme_df=tune_programme,
        working_points_df=working_points,
        resonance_lines_df=diagram["resonance_lines"],
        resonance_proximity_df=proximity,
        warnings_df=warnings_df,
        plot_manifest_df=plot_manifest_df,
        run_paths=run_paths,
        metadata={
            "label": config.label,
            "n_snapshots": len(snapshots),
            "run_paths": asdict(run_paths),
            "diagram": dict(diagram["metadata"]),
            **dict(config.metadata or {}),
        },
        warnings=warnings,
    )


def build_full_cycle_snapshot_series(
    cycle_times_ms,
    qx_values,
    qy_values,
    base_config=None,
    label="full_cycle",
    point_overrides=None,
    **series_kwargs,
):
    """
    Convenience wrapper for a time/tune programme snapshot series.
    """

    cycle_times_ms = list(cycle_times_ms)
    qx_values = list(qx_values)
    qy_values = list(qy_values)
    if not (len(cycle_times_ms) == len(qx_values) == len(qy_values)):
        raise ValueError("cycle_times_ms, qx_values and qy_values must have the same length.")
    point_overrides = [{} for _ in cycle_times_ms] if point_overrides is None else list(point_overrides)
    if len(point_overrides) != len(cycle_times_ms):
        raise ValueError("point_overrides must have the same length as cycle_times_ms.")

    if base_config is None:
        base_config = SnapshotConfig(cycle_time_ms=cycle_times_ms[0])
    if not isinstance(base_config, SnapshotConfig):
        raise TypeError("base_config must be a SnapshotConfig.")

    snapshots = []
    for index, (cycle_time, qx, qy) in enumerate(zip(cycle_times_ms, qx_values, qy_values)):
        overrides = {
            "cycle_time_ms": float(cycle_time),
            "requested_qx": float(qx),
            "requested_qy": float(qy),
            "label": f"{label}_{index:02d}",
            "snapshot_id": f"{label}_{index:03d}",
            "case": label,
            **dict(point_overrides[index] or {}),
        }
        snapshots.append(
            copy_snapshot_config(
                base_config,
                **overrides,
            )
        )

    return build_snapshot_series(
        SnapshotSeriesConfig(
            snapshots=snapshots,
            label=label,
            **series_kwargs,
        )
    )


def snapshot_configs_from_timepoint_table(table, base_config, corrector_table=None):
    """
    Build complete SnapshotConfig objects from a GUI/archiver timepoint table.
    """

    if not isinstance(base_config, SnapshotConfig):
        raise TypeError("base_config must be a SnapshotConfig.")
    df = pd.DataFrame(table).copy()
    if "cycle_time_ms" not in df.columns:
        raise ValueError("timepoint table must contain cycle_time_ms.")
    configs = []
    for index, row in df.reset_index(drop=True).iterrows():
        overrides = {
            "cycle_time_ms": float(row["cycle_time_ms"]),
            "label": row.get("label", f"timepoint_{index:02d}"),
            "snapshot_id": row.get("snapshot_id", f"timepoint_{index:03d}"),
        }
        for source, target in (
            ("requested_qx", "requested_qx"),
            ("set_qx", "requested_qx"),
            ("requested_qy", "requested_qy"),
            ("set_qy", "requested_qy"),
            ("base_qx", "base_qx"),
            ("base_qy", "base_qy"),
            ("tune_method", "tune_method"),
            ("main_magnet_mode", "main_magnet_mode"),
        ):
            if source in row.index and pd.notna(row[source]):
                overrides[target] = row[source]
        if "harmonics" in row.index and isinstance(row["harmonics"], dict):
            overrides["harmonics"] = row["harmonics"]
        if "error_table_paths" in row.index and row["error_table_paths"] is not None:
            value = row["error_table_paths"]
            overrides["error_table_paths"] = value if isinstance(value, list) else [value]
        if corrector_table is not None:
            overrides["corrector_settings"] = corrector_settings_from_dataframe(
                corrector_table,
                cycle_time_ms=float(row["cycle_time_ms"]),
            )
        configs.append(copy_snapshot_config(base_config, **overrides))
    return configs


def config_to_dataframe(configs):
    """
    Return a compact table of snapshot configs for review in notebooks/GUI.
    """

    rows = []
    for index, config in enumerate(configs):
        if not isinstance(config, SnapshotConfig):
            raise TypeError("configs must contain SnapshotConfig objects.")
        row = asdict(config)
        row["index"] = index
        row["label"] = config.resolved_label()
        row["snapshot_id"] = config.resolved_snapshot_id(index)
        rows.append(row)
    return pd.DataFrame(rows)
