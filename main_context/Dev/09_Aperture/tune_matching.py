"""
tune_matching.py

GUI-facing tune matching workflows for the ISIS RCS optics backend.

This layer keeps MAD-X execution inside MadxModel and wraps the existing
MachineState tune controls into small workflows that expose the values the GUI
needs to display: requested tunes, actual MAD-X tunes, matched tunes, trim quad
currents and MAD-X K values.
"""

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import re

import pandas as pd

from cycle_time import RCSRamp
from machine_state import MachineState
from machine_state_defaults import DEFAULT_HARMONICS
from madx_model import MadxModel


CANONICAL_HARMONIC_KEYS = tuple(DEFAULT_HARMONICS.keys())
_HARMONIC_LOOKUP = {
    re.sub(r"[^A-Z0-9]", "", key.upper()): key
    for key in CANONICAL_HARMONIC_KEYS
}


@dataclass
class TuneWorkflowResult:
    """
    Result container for one tune workflow.
    """

    workflow: str
    machine_state: MachineState
    requested_qx: float
    requested_qy: float
    di_wright_qx: float = None
    di_wright_qy: float = None
    actual_qx: float = None
    actual_qy: float = None
    matched_qx: float = None
    matched_qy: float = None
    iqtf_A: float = None
    iqtd_A: float = None
    kqtf: float = None
    kqtd: float = None
    harmonics: OrderedDict = None
    summary_df: pd.DataFrame = None
    twiss_df: pd.DataFrame = None
    metadata: dict = None
    warnings: list = None

    def trim_quad_display_dict(self):
        """
        Return trim quad currents and K values for compact GUI display.
        """

        return {
            "iqtf_A": self.iqtf_A,
            "iqtd_A": self.iqtd_A,
            "kqtf": self.kqtf,
            "kqtd": self.kqtd,
        }

    def to_summary_dataframe(self):
        """
        Return a one-row summary DataFrame for notebook or GUI display.
        """

        return pd.DataFrame(
            [
                {
                    "workflow": self.workflow,
                    "requested_qx": self.requested_qx,
                    "requested_qy": self.requested_qy,
                    "di_wright_qx": self.di_wright_qx,
                    "di_wright_qy": self.di_wright_qy,
                    "actual_qx": self.actual_qx,
                    "actual_qy": self.actual_qy,
                    "matched_qx": self.matched_qx,
                    "matched_qy": self.matched_qy,
                    "iqtf_A": self.iqtf_A,
                    "iqtd_A": self.iqtd_A,
                    "kqtf": self.kqtf,
                    "kqtd": self.kqtd,
                }
            ]
        )


def normalise_harmonic_inputs(harmonics=None, **kwargs):
    """
    Return a complete OrderedDict of canonical harmonic tune settings.

    Accepted input keys are case-insensitive and may include separators, so
    F8Sin, F8SIN, f8sin and F8_SIN all map to F8SIN.
    """

    merged = {}
    if harmonics is not None:
        merged.update(dict(harmonics))
    merged.update(kwargs)

    normalised = OrderedDict(DEFAULT_HARMONICS)

    for key, value in merged.items():
        lookup_key = re.sub(r"[^A-Z0-9]", "", str(key).upper())
        if lookup_key not in _HARMONIC_LOOKUP:
            valid = ", ".join(CANONICAL_HARMONIC_KEYS)
            raise KeyError(f"Unknown harmonic tune variable {key!r}. Valid keys: {valid}")

        normalised[_HARMONIC_LOOKUP[lookup_key]] = float(value)

    return normalised


def detect_explicit_harmonic_support(lattice_folder, strength_file="ISIS.strength"):
    """
    Check whether a lattice strength file exposes the harmonic tune accessors.
    """

    strength_path = Path(lattice_folder) / strength_file
    if not strength_path.is_file():
        raise FileNotFoundError(f"Missing strength file: {strength_path}")

    text = strength_path.read_text()
    present = OrderedDict((key, key in text) for key in CANONICAL_HARMONIC_KEYS)
    missing = [key for key, found in present.items() if not found]

    return {
        "supported": not missing,
        "present": present,
        "missing": missing,
        "strength_file": str(strength_path),
    }


def make_tune_machine_state(
    beam_state,
    requested_qx,
    requested_qy,
    tune_method,
    harmonics=None,
    main_magnet_mode="rcs_bare",
    metadata=None,
    calculate_tune=True,
):
    """
    Build a MachineState with canonical harmonic settings.
    """

    return MachineState.from_defaults(
        beam_state=beam_state,
        main_magnet_mode=main_magnet_mode,
        requested_qx=requested_qx,
        requested_qy=requested_qy,
        tune_method=tune_method,
        harmonic_tunes=normalise_harmonic_inputs(harmonics),
        metadata={} if metadata is None else dict(metadata),
        calculate_tune=calculate_tune,
    )


def _default_beam_state(cycle_time_ms):
    return RCSRamp().state_at(cycle_time_ms)


def _summary_value(summary, key):
    value = summary.get(key, None)
    if value is None:
        return None
    return float(value)


def _run_twiss_for_state(
    machine_state,
    lattice_folder,
    output_dir,
    sequence_name,
    twiss_columns=None,
):
    model = MadxModel(
        lattice_folder=lattice_folder,
        sequence_name=sequence_name,
        output_dir=output_dir,
    )
    model.load_lattice(use_sequence=False)
    model.apply_machine_state(machine_state)
    model.use_sequence()
    twiss_df = model.run_twiss(columns=twiss_columns)
    summary = model.get_summary_dict()
    return model, twiss_df, summary


def match_tune_with_madx(
    requested_qx,
    requested_qy,
    beam_state=None,
    cycle_time_ms=0.0,
    harmonics=None,
    lattice_folder="../Lattice_Files/00_Simplified_Lattice",
    output_dir="./tune_matching_tests/madx_match",
    sequence_name="synchrotron",
    main_magnet_mode="rcs_bare",
    requested_dq1=None,
    requested_dq2=None,
    step=1.0e-4,
    calls=50000,
    tolerance=1.0e-6,
    chrom=True,
    twiss_columns=None,
):
    """
    Use MAD-X matching to find trim quad settings for requested Qx/Qy.
    """

    if beam_state is None:
        beam_state = _default_beam_state(cycle_time_ms)

    support = detect_explicit_harmonic_support(lattice_folder)
    warnings = []
    if not support["supported"]:
        warnings.append(f"Lattice is missing harmonic tune accessors: {support['missing']}")

    machine_state = make_tune_machine_state(
        beam_state=beam_state,
        requested_qx=requested_qx,
        requested_qy=requested_qy,
        tune_method="madx_match",
        harmonics=harmonics,
        main_magnet_mode=main_magnet_mode,
        metadata={"workflow": "match_tune_with_madx"},
        calculate_tune=False,
    )

    model = MadxModel(
        lattice_folder=lattice_folder,
        sequence_name=sequence_name,
        output_dir=output_dir,
    )
    model.load_lattice(use_sequence=False)
    match_result = model.match_tune_from_machine_state(
        machine_state,
        requested_dq1=requested_dq1,
        requested_dq2=requested_dq2,
        step=step,
        calls=calls,
        tolerance=tolerance,
        chrom=chrom,
        apply_machine_state=True,
        include_trim_quads=True,
        run_twiss_after_match=True,
    )
    twiss_df = model.run_twiss(columns=twiss_columns)
    summary_df = model.get_summary_df(refresh=True)

    result = TuneWorkflowResult(
        workflow="match_tune_with_madx",
        machine_state=machine_state,
        requested_qx=float(requested_qx),
        requested_qy=float(requested_qy),
        matched_qx=match_result["matched_qx"],
        matched_qy=match_result["matched_qy"],
        iqtf_A=machine_state.iqtf_A,
        iqtd_A=machine_state.iqtd_A,
        kqtf=machine_state.kqtf,
        kqtd=machine_state.kqtd,
        harmonics=OrderedDict(machine_state.harmonic_tunes),
        summary_df=summary_df,
        twiss_df=twiss_df,
        metadata=model.get_metadata(),
        warnings=warnings,
    )
    result.summary_df = result.to_summary_dataframe()
    return result


def evaluate_di_wright_then_match(
    di_wright_qx,
    di_wright_qy,
    target_qx=None,
    target_qy=None,
    beam_state=None,
    cycle_time_ms=0.0,
    harmonics=None,
    lattice_folder="../Lattice_Files/00_Simplified_Lattice",
    output_dir="./tune_matching_tests/di_wright_then_match",
    sequence_name="synchrotron",
    main_magnet_mode="rcs_bare",
    requested_dq1=None,
    requested_dq2=None,
    step=1.0e-4,
    calls=50000,
    tolerance=1.0e-6,
    chrom=True,
    twiss_columns=None,
):
    """
    Evaluate Di Wright tune settings in MAD-X, then match to target Qx/Qy.

    If target_qx/target_qy are not supplied, the Di Wright input tunes are also
    used as the MAD-X matching target.
    """

    if beam_state is None:
        beam_state = _default_beam_state(cycle_time_ms)

    if target_qx is None:
        target_qx = di_wright_qx
    if target_qy is None:
        target_qy = di_wright_qy

    harmonics = normalise_harmonic_inputs(harmonics)

    support = detect_explicit_harmonic_support(lattice_folder)
    warnings = []
    if not support["supported"]:
        warnings.append(f"Lattice is missing harmonic tune accessors: {support['missing']}")

    di_state = make_tune_machine_state(
        beam_state=beam_state,
        requested_qx=di_wright_qx,
        requested_qy=di_wright_qy,
        tune_method="di_wright",
        harmonics=harmonics,
        main_magnet_mode=main_magnet_mode,
        metadata={"workflow": "evaluate_di_wright_actual"},
    )

    actual_model, actual_twiss_df, actual_summary = _run_twiss_for_state(
        machine_state=di_state,
        lattice_folder=lattice_folder,
        output_dir=str(Path(output_dir) / "actual"),
        sequence_name=sequence_name,
        twiss_columns=twiss_columns,
    )

    match_result = match_tune_with_madx(
        requested_qx=target_qx,
        requested_qy=target_qy,
        beam_state=beam_state,
        harmonics=harmonics,
        lattice_folder=lattice_folder,
        output_dir=str(Path(output_dir) / "matched"),
        sequence_name=sequence_name,
        main_magnet_mode=main_magnet_mode,
        requested_dq1=requested_dq1,
        requested_dq2=requested_dq2,
        step=step,
        calls=calls,
        tolerance=tolerance,
        chrom=chrom,
        twiss_columns=twiss_columns,
    )

    metadata = {
        "actual_model": actual_model.get_metadata(),
        "matched_model": match_result.metadata,
    }

    result = TuneWorkflowResult(
        workflow="evaluate_di_wright_then_match",
        machine_state=match_result.machine_state,
        requested_qx=float(target_qx),
        requested_qy=float(target_qy),
        di_wright_qx=float(di_wright_qx),
        di_wright_qy=float(di_wright_qy),
        actual_qx=_summary_value(actual_summary, "q1"),
        actual_qy=_summary_value(actual_summary, "q2"),
        matched_qx=match_result.matched_qx,
        matched_qy=match_result.matched_qy,
        iqtf_A=match_result.iqtf_A,
        iqtd_A=match_result.iqtd_A,
        kqtf=match_result.kqtf,
        kqtd=match_result.kqtd,
        harmonics=OrderedDict(harmonics),
        summary_df=None,
        twiss_df=match_result.twiss_df,
        metadata=metadata,
        warnings=warnings + match_result.warnings,
    )
    result.summary_df = result.to_summary_dataframe()
    result.metadata["actual_twiss_rows"] = len(actual_twiss_df)

    return result


def evaluate_di_wright_tune_programme(
    cycle_times_ms,
    qx_values,
    qy_values,
    harmonics=None,
    lattice_folder="../Lattice_Files/00_Simplified_Lattice",
    output_dir="./tune_matching_tests/di_wright_programme",
    sequence_name="synchrotron",
    main_magnet_mode="rcs_bare",
    twiss_columns=None,
):
    """
    Evaluate a Di Wright tune programme in MAD-X without tune matching.

    One MachineState is created for each cycle point. The returned table keeps
    the set tunes, actual MAD-X tunes, Di Wright trim currents and Di Wright K
    values together for GUI display and resonance plotting.
    """

    cycle_times_ms = list(cycle_times_ms)
    qx_values = list(qx_values)
    qy_values = list(qy_values)

    if not (len(cycle_times_ms) == len(qx_values) == len(qy_values)):
        raise ValueError("cycle_times_ms, qx_values and qy_values must have the same length.")

    harmonics = normalise_harmonic_inputs(harmonics)
    support = detect_explicit_harmonic_support(lattice_folder)
    warnings = []
    if not support["supported"]:
        warnings.append(f"Lattice is missing harmonic tune accessors: {support['missing']}")

    ramp = RCSRamp()
    rows = []

    for index, (cycle_time_ms, qx, qy) in enumerate(
        zip(cycle_times_ms, qx_values, qy_values)
    ):
        beam_state = ramp.state_at(cycle_time_ms)
        machine_state = make_tune_machine_state(
            beam_state=beam_state,
            requested_qx=qx,
            requested_qy=qy,
            tune_method="di_wright",
            harmonics=harmonics,
            main_magnet_mode=main_magnet_mode,
            metadata={
                "workflow": "evaluate_di_wright_tune_programme",
                "programme_index": index,
            },
        )

        model, twiss_df, summary = _run_twiss_for_state(
            machine_state=machine_state,
            lattice_folder=lattice_folder,
            output_dir=str(Path(output_dir) / f"point_{index:02d}"),
            sequence_name=sequence_name,
            twiss_columns=twiss_columns,
        )

        rows.append(
            {
                "index": index,
                "cycle_time_ms": float(cycle_time_ms),
                "set_qx": float(qx),
                "set_qy": float(qy),
                "actual_qx": _summary_value(summary, "q1"),
                "actual_qy": _summary_value(summary, "q2"),
                "actual_dqx": _summary_value(summary, "dq1"),
                "actual_dqy": _summary_value(summary, "dq2"),
                "iqtf_A": machine_state.iqtf_A,
                "iqtd_A": machine_state.iqtd_A,
                "kqtf": machine_state.kqtf,
                "kqtd": machine_state.kqtd,
                "brho_Tm": float(beam_state.brho_Tm),
                "normalised_momentum": float(beam_state.normalised_momentum),
                "twiss_rows": len(twiss_df),
                "output_dir": model.output_dir,
                "warnings": list(warnings),
            }
        )

    return pd.DataFrame(rows)


def evaluate_di_wright_tune_point(
    cycle_time_ms,
    qx,
    qy,
    harmonics=None,
    lattice_folder="../Lattice_Files/00_Simplified_Lattice",
    output_dir="./tune_matching_tests/di_wright_point",
    sequence_name="synchrotron",
    main_magnet_mode="rcs_bare",
    twiss_columns=None,
):
    """
    Evaluate one Di Wright tune point in MAD-X without tune matching.

    This is the single-point companion to evaluate_di_wright_tune_programme and
    returns TWISS data for plotting beta-function changes.
    """

    beam_state = RCSRamp().state_at(cycle_time_ms)
    harmonics = normalise_harmonic_inputs(harmonics)
    warnings = []
    support = detect_explicit_harmonic_support(lattice_folder)
    if not support["supported"]:
        warnings.append(f"Lattice is missing harmonic tune accessors: {support['missing']}")

    machine_state = make_tune_machine_state(
        beam_state=beam_state,
        requested_qx=qx,
        requested_qy=qy,
        tune_method="di_wright",
        harmonics=harmonics,
        main_magnet_mode=main_magnet_mode,
        metadata={"workflow": "evaluate_di_wright_tune_point"},
    )

    model, twiss_df, summary = _run_twiss_for_state(
        machine_state=machine_state,
        lattice_folder=lattice_folder,
        output_dir=output_dir,
        sequence_name=sequence_name,
        twiss_columns=twiss_columns,
    )

    result = TuneWorkflowResult(
        workflow="evaluate_di_wright_tune_point",
        machine_state=machine_state,
        requested_qx=float(qx),
        requested_qy=float(qy),
        di_wright_qx=float(qx),
        di_wright_qy=float(qy),
        actual_qx=_summary_value(summary, "q1"),
        actual_qy=_summary_value(summary, "q2"),
        iqtf_A=machine_state.iqtf_A,
        iqtd_A=machine_state.iqtd_A,
        kqtf=machine_state.kqtf,
        kqtd=machine_state.kqtd,
        harmonics=OrderedDict(harmonics),
        twiss_df=twiss_df,
        metadata=model.get_metadata(),
        warnings=warnings,
    )
    result.summary_df = result.to_summary_dataframe()
    return result
