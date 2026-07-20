"""
orbit_correction.py

Read-only COCU-style closed-orbit correction workflows.

This layer owns measurement normalisation, BPM/corrector selection, correction
suggestion tables and residual summaries. MAD-X execution remains inside
MadxModel, and operational corrector current conversion remains in correctors.py.
"""

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np
import pandas as pd

from .correctors import (
    kick_rad_to_current,
)
from .machine_state_defaults import (
    CORRECTOR_SUPERPERIODS,
    HD_CORRECTOR_NAMES,
    VD_CORRECTOR_NAMES,
)
from .madx_model import MadxModel


BPM_COLUMNS = (
    "bpm",
    "plane",
    "closed_orbit_mm",
    "closed_orbit_mm_err",
    "s",
    "enabled",
)

CORRECTOR_COLUMNS = (
    "corrector",
    "plane",
    "superperiod",
    "enabled",
    "initial_kick_rad",
    "matched_kick_rad",
    "delta_kick_rad",
    "delta_kick_mrad",
    "delta_current_A",
)


@dataclass
class MadxCorrectResult:
    """
    Two-stage MAD-X orbit correction result.

    A small set of BPM measurements is first used to fit a viable measured
    orbit with non-steering knobs. MAD-X CORRECT is then used to propose
    steering-dipole corrections against the reference orbit.
    """

    plane: str
    reference_twiss_df: pd.DataFrame
    measured_twiss_df: pd.DataFrame
    corrected_twiss_df: pd.DataFrame
    bpm_measurements: pd.DataFrame
    seed_constraints: list
    seed_fit_knobs: list
    seed_fit_kicks_rad: dict
    correctors: pd.DataFrame
    correction_table_df: pd.DataFrame
    monitor_summary_before: dict
    monitor_summary_after: dict
    summary: dict
    warnings: list
    metadata: dict


class MeasurementSource:
    """
    Read-only source interface for future machine/EPICS measurements.
    """

    def read(self):
        raise NotImplementedError("MeasurementSource subclasses must implement read().")


class DataFrameMeasurementSource(MeasurementSource):
    """
    Read-only measurement source backed by an in-memory DataFrame.
    """

    def __init__(self, dataframe):
        self.dataframe = normalise_bpm_measurements(dataframe)

    def read(self):
        return self.dataframe.copy()


def normalise_plane(plane):
    plane = str(plane).upper()
    if plane in ("H", "X", "HORIZONTAL"):
        return "H"
    if plane in ("V", "Y", "VERTICAL"):
        return "V"
    raise ValueError("plane must be one of H, V, horizontal or vertical.")


def plane_coordinate(plane):
    return "x" if normalise_plane(plane) == "H" else "y"


def normalise_bpm_measurements(data, enabled_default=True):
    """
    Return a BPM measurement DataFrame with canonical columns.
    """

    df = pd.DataFrame(data).copy()

    required = ["bpm", "plane", "closed_orbit_mm"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(f"BPM measurements are missing columns: {missing}")

    if "closed_orbit_mm_err" not in df.columns:
        df["closed_orbit_mm_err"] = np.nan
    if "s" not in df.columns:
        df["s"] = np.nan
    if "enabled" not in df.columns:
        df["enabled"] = bool(enabled_default)

    df["bpm"] = df["bpm"].astype(str)
    df["plane"] = df["plane"].map(normalise_plane)
    df["closed_orbit_mm"] = pd.to_numeric(df["closed_orbit_mm"], errors="coerce")
    df["closed_orbit_mm_err"] = pd.to_numeric(df["closed_orbit_mm_err"], errors="coerce")
    df["s"] = pd.to_numeric(df["s"], errors="coerce")
    df["enabled"] = df["enabled"].astype(bool)

    if not np.isfinite(df["closed_orbit_mm"].to_numpy()).all():
        raise ValueError("closed_orbit_mm contains non-finite values.")

    return df.loc[:, BPM_COLUMNS].copy()


def bpm_measurements_from_twiss(
    twiss_df,
    plane,
    monitor_pattern=None,
    error_mm=0.1,
    enabled_default=True,
):
    """
    Build target BPM measurements from a TWISS table at monitor elements.
    """

    if not isinstance(twiss_df, pd.DataFrame):
        raise TypeError("twiss_df must be a pandas DataFrame.")

    plane = normalise_plane(plane)
    coordinate = plane_coordinate(plane)

    if monitor_pattern is None:
        monitor_pattern = "hm" if plane == "H" else "vm"

    required = ["name", "s", coordinate]
    missing = [column for column in required if column not in twiss_df.columns]
    if missing:
        raise ValueError(f"TWISS DataFrame is missing columns: {missing}")

    monitors = twiss_df[
        twiss_df["name"].astype(str).str.contains(monitor_pattern, case=False, na=False)
    ].copy()

    rows = []
    for _, row in monitors.iterrows():
        rows.append(
            {
                "bpm": str(row["name"]).split(":")[0],
                "plane": plane,
                "closed_orbit_mm": float(row[coordinate]) * 1.0e3,
                "closed_orbit_mm_err": float(error_mm),
                "s": float(row["s"]),
                "enabled": bool(enabled_default),
            }
        )

    return normalise_bpm_measurements(rows)


def set_bpm_enabled(bpm_measurements, bpm_names, enabled):
    """
    Return a copy with selected BPM rows enabled or disabled.
    """

    df = normalise_bpm_measurements(bpm_measurements)
    names = {str(name).lower() for name in bpm_names}
    mask = df["bpm"].str.lower().isin(names)
    df.loc[mask, "enabled"] = bool(enabled)
    return df


def default_corrector_selection(plane=None, enabled_default=True):
    """
    Return the operational HD/VD corrector selection table.
    """

    rows = []
    requested_plane = None if plane is None else normalise_plane(plane)

    for corrector_plane, names in (("H", HD_CORRECTOR_NAMES), ("V", VD_CORRECTOR_NAMES)):
        if requested_plane is not None and corrector_plane != requested_plane:
            continue
        for name in names:
            rows.append(
                {
                    "corrector": name,
                    "plane": corrector_plane,
                    "superperiod": _superperiod_from_corrector_name(name),
                    "enabled": bool(enabled_default),
                    "initial_kick_rad": 0.0,
                    "matched_kick_rad": 0.0,
                    "delta_kick_rad": 0.0,
                    "delta_kick_mrad": 0.0,
                    "delta_current_A": np.nan,
                }
            )

    return pd.DataFrame(rows, columns=CORRECTOR_COLUMNS)


def normalise_corrector_selection(data=None, plane=None, enabled_default=True):
    """
    Return a corrector selection DataFrame with canonical columns.
    """

    if data is None:
        return default_corrector_selection(plane=plane, enabled_default=enabled_default)

    df = pd.DataFrame(data).copy()
    if "corrector" not in df.columns:
        raise ValueError("Corrector selection is missing column: corrector")

    if "plane" not in df.columns:
        df["plane"] = df["corrector"].map(_plane_from_corrector_name)
    if "superperiod" not in df.columns:
        df["superperiod"] = df["corrector"].map(_superperiod_from_corrector_name)
    if "enabled" not in df.columns:
        df["enabled"] = bool(enabled_default)

    for column in CORRECTOR_COLUMNS:
        if column not in df.columns:
            df[column] = np.nan if column.endswith("_A") else 0.0

    df["corrector"] = df["corrector"].astype(str)
    df["plane"] = df["plane"].map(normalise_plane)
    df["superperiod"] = pd.to_numeric(df["superperiod"], errors="coerce").astype(int)
    df["enabled"] = df["enabled"].astype(bool)

    if plane is not None:
        df = df[df["plane"] == normalise_plane(plane)].copy()

    return df.loc[:, CORRECTOR_COLUMNS].reset_index(drop=True)


def set_corrector_enabled(correctors, corrector_names, enabled):
    """
    Return a copy with selected corrector rows enabled or disabled.
    """

    df = normalise_corrector_selection(correctors)
    names = {str(name).lower() for name in corrector_names}
    mask = df["corrector"].str.lower().isin(names)
    df.loc[mask, "enabled"] = bool(enabled)
    return df


def enabled_bpm_measurements(bpm_measurements, plane):
    df = normalise_bpm_measurements(bpm_measurements)
    plane = normalise_plane(plane)
    return df[(df["plane"] == plane) & (df["enabled"])].reset_index(drop=True)


def enabled_corrector_names(correctors, plane):
    df = normalise_corrector_selection(correctors, plane=plane)
    return df[df["enabled"]]["corrector"].astype(str).tolist()


def resolve_bpm_to_twiss_name(bpm_label, twiss_df, s_value=None):
    """
    Resolve a BPM label to a MAD-X element name from a TWISS DataFrame.
    """

    if twiss_df is None or "name" not in twiss_df.columns:
        return None

    label = str(bpm_label).lower()
    lookup = twiss_df.copy()
    lookup["name_str"] = lookup["name"].astype(str)
    mask = lookup["name_str"].str.lower().str.contains(re.escape(label), regex=True)
    candidates = lookup.loc[mask].copy()

    if candidates.empty:
        return None

    if "s" in candidates.columns and s_value is not None and np.isfinite(s_value):
        candidates["ds"] = np.abs(candidates["s"].astype(float) - float(s_value))
        selected = candidates.sort_values("ds").iloc[0]["name_str"]
    else:
        selected = candidates.iloc[0]["name_str"]

    return str(selected).split(":")[0]


def constraints_from_bpm_measurements(
    bpm_measurements,
    plane,
    twiss_df=None,
    min_error_mm=0.1,
    use_error_bounds=True,
):
    """
    Build MAD-X orbit constraints from enabled BPM measurements.
    """

    df = enabled_bpm_measurements(bpm_measurements, plane)
    coordinate = plane_coordinate(plane)
    constraints = []
    warnings = []

    constraint_type = None
    if use_error_bounds:
        try:
            from cpymad.types import Constraint
        except ImportError:
            Constraint = None
            warnings.append(
                "cpymad.types.Constraint is unavailable; using exact orbit constraints."
            )
        constraint_type = Constraint

    for _, row in df.iterrows():
        bpm = row["bpm"]
        resolved = resolve_bpm_to_twiss_name(bpm, twiss_df, row["s"])
        if resolved is None:
            warnings.append(f"No TWISS element matched BPM {bpm!r}; skipping it.")
            continue

        value_m = float(row["closed_orbit_mm"]) * 1.0e-3
        error_mm = row["closed_orbit_mm_err"]
        if not np.isfinite(error_mm):
            error_mm = min_error_mm
        error_m = max(float(error_mm), float(min_error_mm)) * 1.0e-3

        if constraint_type is not None:
            target = constraint_type(min=value_m - error_m, max=value_m + error_m)
        else:
            target = value_m

        constraints.append({"range": resolved, coordinate: target})

    if not constraints:
        raise ValueError("No enabled BPM measurements could be converted to constraints.")

    return constraints, warnings


def orbit_residuals(twiss_df, bpm_measurements, plane):
    """
    Compare model TWISS orbit against enabled BPM measurements.
    """

    df = enabled_bpm_measurements(bpm_measurements, plane)
    coordinate = plane_coordinate(plane)
    rows = []

    for _, row in df.iterrows():
        resolved = resolve_bpm_to_twiss_name(row["bpm"], twiss_df, row["s"])
        model_mm = np.nan
        residual_mm = np.nan

        if resolved is not None:
            mask = twiss_df["name"].astype(str).str.lower().str.split(":").str[0] == resolved.lower()
            matches = twiss_df.loc[mask]
            if matches.empty:
                matches = twiss_df[
                    twiss_df["name"].astype(str).str.lower().str.contains(
                        re.escape(str(row["bpm"]).lower()),
                        regex=True,
                    )
                ]
            if not matches.empty and coordinate in matches.columns:
                model_mm = float(matches.iloc[0][coordinate]) * 1.0e3
                residual_mm = model_mm - float(row["closed_orbit_mm"])

        rows.append(
            {
                "bpm": row["bpm"],
                "plane": normalise_plane(plane),
                "enabled": bool(row["enabled"]),
                "matched_name": resolved,
                "s": row["s"],
                "measurement_mm": float(row["closed_orbit_mm"]),
                "model_mm": model_mm,
                "residual_mm": residual_mm,
            }
        )

    return pd.DataFrame(rows)


def default_orbit_fit_knobs(plane, source="dipole"):
    """
    Return non-steering lattice kick knobs for fitting a viable measured orbit.

    These are deliberately separate from the operational HD/VD steering
    correctors used by MAD-X CORRECT.
    """

    plane = normalise_plane(plane)
    source = str(source).lower()

    if source != "dipole":
        raise ValueError("Only source='dipole' is currently supported.")

    suffix = "HKICK" if plane == "H" else "VKICK"
    return [f"R{superperiod}DIP_{suffix}" for superperiod in range(10)]


def correct_bpm_measurements_with_madx_correct(
    plane,
    lattice_folder,
    bpm_measurements,
    machine_state,
    correctors=None,
    fit_knobs=None,
    sequence_name="synchrotron",
    output_dir="./orbit_correction_runs/madx_correct",
    reference_table="bare",
    correction_mode="svd",
    correction_cond=1,
    correction_ncorr=0,
    monitor_pattern=None,
    max_fit_kick_rad=0.005,
    fit_step=1.0e-4,
    fit_calls=50000,
    fit_tolerance=1.0e-6,
    use_error_bounds=False,
    twiss_columns=None,
    metadata=None,
):
    """
    Fit a plausible measured orbit, then run MAD-X CORRECT with steering dipoles.

    This reproduces the operational split more closely than using steering
    correctors for both stages:
        1. fit BPM measurement data with non-steering lattice/error knobs;
        2. use MAD-X CORRECT to find selected HD/VD steering corrections;
        3. convert those MAD-X kicks to currents with ISIS calibrations.
    """

    if machine_state is None:
        raise ValueError("machine_state is required for this workflow.")

    plane = normalise_plane(plane)
    bpm_measurements = normalise_bpm_measurements(bpm_measurements)
    correctors = normalise_corrector_selection(correctors, plane=plane)
    vary_correctors = enabled_corrector_names(correctors, plane)
    if not vary_correctors:
        raise ValueError(f"No enabled {plane} correctors selected.")

    if fit_knobs is None:
        fit_knobs = default_orbit_fit_knobs(plane)
    fit_knobs = [str(name) for name in fit_knobs]

    output_dir = Path(output_dir)
    model = MadxModel(
        lattice_folder=lattice_folder,
        sequence_name=sequence_name,
        output_dir=output_dir / plane.lower() / "madx",
    )
    model.load_lattice(use_sequence=False)
    model.apply_machine_state(machine_state)
    model.use_sequence()

    model.madx.twiss(sequence=sequence_name, table=reference_table)
    reference_twiss = getattr(model.madx.table, reference_table).dframe().reset_index(drop=True)
    reference_twiss = model._normalise_dataframe_columns(reference_twiss)
    if twiss_columns is not None:
        reference_twiss = reference_twiss[twiss_columns]

    seed_constraints, warnings = constraints_from_bpm_measurements(
        bpm_measurements=bpm_measurements,
        plane=plane,
        twiss_df=reference_twiss,
        use_error_bounds=use_error_bounds,
    )

    seed_match = model.match_orbit(
        constraints=seed_constraints,
        vary_names=fit_knobs,
        horizontal=(plane == "H"),
        step=fit_step,
        calls=fit_calls,
        tolerance=fit_tolerance,
        max_kick=max_fit_kick_rad,
        run_twiss_after_match=True,
        twiss_columns=twiss_columns,
    )
    measured_twiss = seed_match["twiss_df"]

    correct_result = model.correct_orbit(
        plane=plane,
        corrector_names=vary_correctors,
        monitor_pattern=monitor_pattern,
        model_table=reference_table,
        output_prefix=output_dir / plane.lower() / f"madx_correct_{plane.lower()}",
        mode=correction_mode,
        cond=correction_cond,
        ncorr=correction_ncorr,
        run_twiss_after_correct=True,
        twiss_columns=twiss_columns,
    )
    corrected_twiss = correct_result["twiss_df"]

    correction_table = read_madx_correct_table(correct_result["clist_path"])
    corrector_result = corrector_suggestions_from_madx_correct_table(
        correction_table=correction_table,
        correctors=correctors,
        plane=plane,
        beam_state=machine_state.beam_state,
    )

    before_summary = monitor_orbit_summary(measured_twiss, plane, monitor_pattern=monitor_pattern)
    after_summary = monitor_orbit_summary(corrected_twiss, plane, monitor_pattern=monitor_pattern)

    summary = {
        "plane": plane,
        "read_only": True,
        "workflow": "bpm_match_then_madx_correct",
        "n_seed_bpm": int(len(enabled_bpm_measurements(bpm_measurements, plane))),
        "n_seed_fit_knobs": int(len(fit_knobs)),
        "n_enabled_correctors": int(len(vary_correctors)),
        "before": before_summary,
        "after": after_summary,
        "rms_change_mm": _rms_change(before_summary, after_summary),
    }

    result_metadata = dict(metadata or {})
    result_metadata.update(model.get_metadata())
    result_metadata.update(
        {
            "seed_fit_knobs": list(fit_knobs),
            "corrector_names": list(vary_correctors),
            "correction_mode": correction_mode,
            "correction_cond": correction_cond,
            "correction_ncorr": correction_ncorr,
            "clist_path": correct_result["clist_path"],
            "mlist_path": correct_result["mlist_path"],
            "machine_write": False,
        }
    )

    return MadxCorrectResult(
        plane=plane,
        reference_twiss_df=reference_twiss,
        measured_twiss_df=measured_twiss,
        corrected_twiss_df=corrected_twiss,
        bpm_measurements=bpm_measurements,
        seed_constraints=seed_constraints,
        seed_fit_knobs=fit_knobs,
        seed_fit_kicks_rad=seed_match["matched_kicks_rad"],
        correctors=corrector_result,
        correction_table_df=correction_table,
        monitor_summary_before=before_summary,
        monitor_summary_after=after_summary,
        summary=summary,
        warnings=warnings,
        metadata=result_metadata,
    )


def read_madx_correct_table(path):
    """
    Read a MAD-X CORRECT clist/mlist TFS-like table into a DataFrame.
    """

    header = None
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("*"):
                header = stripped.split()[1:]
                continue
            if stripped.startswith("$") or stripped.startswith("@"):
                continue
            if stripped.startswith('"') and header is not None:
                parts = stripped.replace('"', "").split()
                if len(parts) >= len(header):
                    rows.append(dict(zip(header, parts[: len(header)])))

    df = pd.DataFrame(rows)
    for column in df.columns:
        if column != "NAME":
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df


def corrector_suggestions_from_madx_correct_table(
    correction_table,
    correctors,
    plane,
    beam_state,
):
    """
    Map MAD-X CORRECT PX/PY.CORRECTION rows onto operational correctors.
    """

    plane = normalise_plane(plane)
    correction_column = "PX.CORRECTION" if plane == "H" else "PY.CORRECTION"
    df = normalise_corrector_selection(correctors, plane=plane)
    table = pd.DataFrame(correction_table).copy()

    correction_lookup = {}
    if not table.empty and "NAME" in table.columns and correction_column in table.columns:
        for _, row in table.iterrows():
            name = str(row["NAME"]).lower()
            value = row[correction_column]
            for corrector_name in df["corrector"]:
                element_name = str(corrector_name).replace("_kick", "").lower()
                if name.endswith(element_name):
                    correction_lookup[corrector_name] = float(value)

    rows = []
    for _, row in df.iterrows():
        corrector_name = row["corrector"]
        delta = correction_lookup.get(corrector_name, 0.0)
        delta_current = np.nan
        if beam_state is not None:
            try:
                delta_current = kick_rad_to_current(corrector_name, delta, beam_state)
            except Exception:
                delta_current = np.nan

        out = dict(row)
        out["initial_kick_rad"] = 0.0
        out["matched_kick_rad"] = float(delta)
        out["delta_kick_rad"] = float(delta)
        out["delta_kick_mrad"] = 1.0e3 * float(delta)
        out["delta_current_A"] = delta_current
        rows.append(out)

    return pd.DataFrame(rows, columns=CORRECTOR_COLUMNS)


def monitor_orbit_summary(twiss_df, plane, monitor_pattern=None):
    """
    Summarise the closed orbit at monitor elements for one plane.
    """

    plane = normalise_plane(plane)
    coordinate = plane_coordinate(plane)
    if monitor_pattern is None:
        monitor_pattern = "hm" if plane == "H" else "vm"

    monitors = twiss_df[
        twiss_df["name"].astype(str).str.contains(monitor_pattern, case=False, na=False)
    ].copy()
    if monitors.empty:
        return {
            "n_monitor": 0,
            "rms_orbit_mm": np.nan,
            "max_abs_orbit_mm": np.nan,
            "mean_orbit_mm": np.nan,
        }

    values = 1.0e3 * monitors[coordinate].astype(float).to_numpy()
    return {
        "n_monitor": int(len(values)),
        "rms_orbit_mm": float(np.sqrt(np.mean(values * values))),
        "max_abs_orbit_mm": float(np.max(np.abs(values))),
        "mean_orbit_mm": float(np.mean(values)),
    }


def plot_orbit_with_bpm(
    twiss_df,
    bpm_measurements,
    plane,
    ax=None,
    label="Model orbit",
    title=None,
    enabled_only=False,
    orbit_kwargs=None,
):
    """
    Plot one model orbit with selected BPM measurements overlaid.
    """

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    plane = normalise_plane(plane)
    coordinate = plane_coordinate(plane)
    orbit_kwargs = {} if orbit_kwargs is None else dict(orbit_kwargs)

    required = ["s", coordinate]
    missing = [column for column in required if column not in twiss_df.columns]
    if missing:
        raise ValueError(f"TWISS DataFrame is missing columns: {missing}")

    bpm_df = normalise_bpm_measurements(bpm_measurements)
    bpm_df = bpm_df[bpm_df["plane"] == plane].copy()
    if enabled_only:
        bpm_df = bpm_df[bpm_df["enabled"]].copy()

    ax.plot(
        twiss_df["s"],
        1.0e3 * twiss_df[coordinate],
        label=label,
        **orbit_kwargs,
    )

    for enabled, marker, color, alpha, bpm_label in (
        (True, "o", "black", 1.0, "Enabled BPM"),
        (False, "x", "0.5", 0.8, "Disabled BPM"),
    ):
        subset = bpm_df[bpm_df["enabled"] == enabled]
        if subset.empty:
            continue
        ax.errorbar(
            subset["s"],
            subset["closed_orbit_mm"],
            yerr=subset["closed_orbit_mm_err"],
            fmt=marker,
            color=color,
            alpha=alpha,
            capsize=3,
            linestyle="none",
            label=bpm_label,
        )

    ax.axhline(0.0, color="0.7", linewidth=0.8)
    ax.set_xlabel("s [m]")
    ax.set_ylabel(f"{coordinate} [mm]")
    ax.grid(True, which="both", linestyle=":", linewidth=0.6)
    ax.set_title(title or f"{plane} closed orbit")
    ax.legend(loc="best")
    return ax


def plot_corrector_suggestions(correctors, ax=None, value="delta_current_A", title=None):
    """
    Plot enabled corrector suggestions as currents or kicks.
    """

    import matplotlib.pyplot as plt

    if value not in ("delta_current_A", "delta_kick_mrad", "delta_kick_rad"):
        raise ValueError("value must be delta_current_A, delta_kick_mrad or delta_kick_rad.")

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    df = normalise_corrector_selection(correctors)
    enabled = df[df["enabled"]].copy()

    colors = np.where(enabled[value].astype(float) >= 0.0, "tab:red", "tab:blue")
    positions = np.arange(len(enabled))
    ax.bar(positions, enabled[value].astype(float), color=colors, alpha=0.85)
    ax.axhline(0.0, color="0.2", linewidth=0.8)
    ax.set_xticks(positions)
    ax.set_xticklabels(enabled["corrector"], rotation=45, ha="right")
    ax.set_ylabel(_suggestion_axis_label(value))
    ax.set_title(title or "Corrector suggestions")
    ax.grid(True, axis="y", linestyle=":", linewidth=0.6)
    return ax


def _superperiod_from_corrector_name(corrector_name):
    match = re.match(r"r([0-9])", str(corrector_name).lower())
    if not match:
        raise ValueError(f"Cannot parse superperiod from corrector name: {corrector_name}")
    superperiod = int(match.group(1))
    if superperiod not in CORRECTOR_SUPERPERIODS:
        raise ValueError(f"No operational steering corrector in superperiod {superperiod}.")
    return superperiod


def _plane_from_corrector_name(corrector_name):
    lower = str(corrector_name).lower()
    if "hd" in lower:
        return "H"
    if "vd" in lower:
        return "V"
    raise ValueError(f"Cannot parse plane from corrector name: {corrector_name}")


def _rms_change(before_summary, after_summary):
    before = before_summary.get("rms_residual_mm", before_summary.get("rms_orbit_mm", np.nan))
    after = after_summary.get("rms_residual_mm", after_summary.get("rms_orbit_mm", np.nan))
    if not np.isfinite(before) or not np.isfinite(after):
        return np.nan
    return float(after - before)


def _suggestion_axis_label(value):
    labels = {
        "delta_current_A": "delta current [A]",
        "delta_kick_mrad": "delta kick [mrad]",
        "delta_kick_rad": "delta kick [rad]",
    }
    return labels[value]
