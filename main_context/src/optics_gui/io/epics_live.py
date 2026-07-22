"""
epics_live.py

Live/archiver-facing helpers for the ISIS RCS optics GUI backend.

This module owns the "go get a real value from EPICS or the archiver"
boundary. It must not perform MAD-X execution or physics calculations --
it produces clean tables/objects for the existing io/measurements.py
adapters (and orbit_correction.py's BPM helpers) to consume unchanged.

EPICS/archiver access itself is injected as callables (fetch_value,
list_available_times) rather than called directly from here, because the
connection method (Channel Access vs. archiver HTTP API) is not yet
confirmed. This keeps every pure lookup rule in this file unit-testable
without a live network connection, and lets the real client be plugged in
later without changing this logic.
"""

import pandas as pd


def _format_ms(value):
    """
    Format a cycle-time-in-ms value to match the archiver's PV-name suffix
    convention, e.g. 0.0 -> "0", 2.5 -> "2.5", -0.4 -> "-.4", 0.5 -> ".5".
    """

    value = float(value)
    if value == int(value):
        return str(int(value))

    text = f"{value:g}"
    if text.startswith("0."):
        text = text[1:]
    elif text.startswith("-0."):
        text = "-" + text[2:]
    return text


def nearest_cycle_time(available_times_ms, requested_cycle_time_ms):
    """
    Return the latest available cycle time at or before the requested one.

    A magnet holds its last commanded current until something changes it,
    so "the most recent sample at or before the requested time" is the
    physically correct choice, not just a convenient fallback.

    Returns None if no available time is at or before the requested time
    (nothing to fall back to) -- callers must not silently substitute a
    later sample in that case.
    """

    requested = float(requested_cycle_time_ms)
    candidates = sorted(float(t) for t in available_times_ms if float(t) <= requested)
    if not candidates:
        return None
    return candidates[-1]


# ----------------------------------------------------------------------
# Correctors -- HD/VD steering, confirmed against the lattice and the
# archiver PV dump. Fully unblocked: no open questions remain here.
# ----------------------------------------------------------------------

from ..machine_state_defaults import CORRECTOR_SUPERPERIODS  # noqa: E402
from .measurements import corrector_settings_from_table  # noqa: E402

DEFAULT_HD_IOC = "DWHST_TEST"
DEFAULT_VD_IOC = "DWVST_TEST"


def corrector_device_name(superperiod, family):
    """
    Return the package-style corrector device name, e.g. "r0hd1".
    """

    family = str(family).upper()
    if family not in ("HD", "VD"):
        raise ValueError("family must be 'HD' or 'VD'.")
    return f"r{int(superperiod)}{family.lower()}1"


def corrector_pv_name(superperiod, family, cycle_time_ms, hd_ioc=DEFAULT_HD_IOC, vd_ioc=DEFAULT_VD_IOC):
    """
    Return the archiver PV name for one corrector at one cycle time.
    """

    family = str(family).upper()
    ioc = hd_ioc if family == "HD" else vd_ioc
    device = f"R{int(superperiod)}{family}1"
    return f"{ioc}::{device}:CURRENT:{_format_ms(cycle_time_ms)}MS"


def get_corrector_settings(
    cycle_time_ms,
    fetch_value,
    list_available_times,
    as_of=None,
    prefer="currents",
    superperiods=CORRECTOR_SUPERPERIODS,
    hd_ioc=DEFAULT_HD_IOC,
    vd_ioc=DEFAULT_VD_IOC,
):
    """
    Build a SnapshotCorrectorSettings object from live/archived corrector currents.

    fetch_value(pv_name, as_of=None) -> float
        Reads one PV's value: live if as_of is None, else at/just before that time.
    list_available_times(device, family) -> iterable of float
        Lists the cycle-time-ms suffixes actually archived for one magnet.

    Both are injected rather than called directly here, because the EPICS/
    archiver connection method is not yet confirmed -- swap in the real
    client without changing this function or its tests.

    Returns (settings, missing): settings is a SnapshotCorrectorSettings, or
    None if every corrector was missing (nothing to build a settings object
    from -- the caller decides how to surface that, rather than getting an
    unrelated-looking error from the table adapter underneath). missing is
    a list of (superperiod, family) pairs that had no available sample at
    or before cycle_time_ms, so an incomplete-but-nonempty snapshot is
    still visible rather than failing silently.
    """

    rows = []
    missing = []

    for superperiod in superperiods:
        for family in ("HD", "VD"):
            device = corrector_device_name(superperiod, family)
            available = list_available_times(device, family)
            resolved_time = nearest_cycle_time(available, cycle_time_ms)
            if resolved_time is None:
                missing.append((superperiod, family))
                continue

            pv_name = corrector_pv_name(
                superperiod, family, resolved_time, hd_ioc=hd_ioc, vd_ioc=vd_ioc
            )
            current_A = fetch_value(pv_name, as_of=as_of)

            rows.append(
                {
                    "cycle_time_ms": float(cycle_time_ms),
                    "device": device,
                    "plane": "H" if family == "HD" else "V",
                    "current_A": float(current_A),
                    "resolved_cycle_time_ms": resolved_time,
                    "pv_name": pv_name,
                }
            )

    if not rows:
        return None, missing

    settings = corrector_settings_from_table(
        rows,
        cycle_time_ms=cycle_time_ms,
        prefer=prefer,
        source="epics_archiver",
    )
    return settings, missing


# ----------------------------------------------------------------------
# BPM / measured orbit -- geometry is confirmed against the lattice, but
# the real position-readback PVs are NOT confirmed yet. The CHANGE:ORBIT_*
# PVs found in the archiver dump sit on the same test IOC as the correctors
# and look like distortion setpoints, not measurements -- don't wire
# pv_name_for_bpm to those until staff confirm the real readback PVs.
# ----------------------------------------------------------------------

from ..orbit_correction import bpm_measurements_from_twiss  # noqa: E402
from .measurements import normalise_bpm_table  # noqa: E402


def bpm_geometry_table(twiss_df, planes=("H", "V")):
    """
    Build the static bpm/plane/s geometry table from a TWISS DataFrame.

    This needs no EPICS/archiver access at all -- BPM positions are fixed
    lattice geometry, not something to query live. Run this once per
    lattice (or cache the result) rather than on every request; the
    closed_orbit_mm column from bpm_measurements_from_twiss is the model's
    prediction and is discarded here, not the real measurement.
    """

    frames = [bpm_measurements_from_twiss(twiss_df, plane=plane) for plane in planes]
    geometry = pd.concat(frames, ignore_index=True)
    return geometry.loc[:, ["bpm", "plane", "s"]].copy()


def get_bpm_measurements(geometry_table, fetch_value, pv_name_for_bpm, as_of=None, enabled_default=True):
    """
    Build a normalised BPM measurement table from live/archived closed-orbit readings.

    geometry_table: DataFrame with bpm/plane/s columns, e.g. from
        bpm_geometry_table(...). Static lattice geometry, never needs EPICS.
    fetch_value(pv_name, as_of=None) -> float
        Reads one BPM's closed_orbit_mm reading.
    pv_name_for_bpm(bpm_label, plane) -> str
        Maps a bpm label (e.g. "sp0_r0hm1") to its archiver PV name. This is
        the piece that is still blocked -- see the module note above.

    Returns a DataFrame in the canonical BPM shape (bpm, plane,
    closed_orbit_mm, closed_orbit_mm_err, s, enabled), ready for
    normalise_bpm_table / orbit correction / measured-orbit display.
    """

    rows = []
    for _, geo_row in geometry_table.iterrows():
        pv_name = pv_name_for_bpm(geo_row["bpm"], geo_row["plane"])
        closed_orbit_mm = fetch_value(pv_name, as_of=as_of)
        rows.append(
            {
                "bpm": geo_row["bpm"],
                "plane": geo_row["plane"],
                "closed_orbit_mm": float(closed_orbit_mm),
                "s": geo_row["s"],
                "enabled": enabled_default,
            }
        )

    return normalise_bpm_table(rows, enabled_default=enabled_default)


# ----------------------------------------------------------------------
# Trim quads / requested tune -- confirmed there is no "target tune" PV
# on the real machine (searched the full archiver PV list for TUNE, QX,
# QY, QH, QV, WORKING_POINT, SETPOINT, PROGRAM, RESONANCE, TARGET: zero
# matches). Tune is only ever set indirectly via trim-quad currents, so
# A3's set_qx/set_qy has to be derived from real QTD/QTF readings by
# reversing the model's own tune-control equations.
#
# NOT CONFIRMED WITH STAFF: there are 10 independent QTD/QTF pairs (one
# per superperiod), but the model's tune-control equations expect a
# single global iqtf_A/iqtd_A pair. get_requested_tune() averages across
# superperiods to bridge that gap -- that averaging is a physics
# assumption, not a confirmed rule. Treat the derived tune as a rough
# estimate until staff confirm how the real per-superperiod currents
# should combine.
# ----------------------------------------------------------------------

from ..machine_state_defaults import DEFAULT_BASE_QX, DEFAULT_BASE_QY  # noqa: E402
from ..tune_control import trim_quad_current_to_tune_di  # noqa: E402

DEFAULT_QT_IOC = "DWQ_TEST"
TRIM_QUAD_SUPERPERIODS = range(10)


def trim_quad_device_name(superperiod, family):
    """
    Return the package-style trim-quad device name, e.g. "r0qtd".
    """

    family = str(family).upper()
    if family not in ("QTD", "QTF"):
        raise ValueError("family must be 'QTD' or 'QTF'.")
    return f"r{int(superperiod)}{family.lower()}"


def trim_quad_pv_name(superperiod, family, cycle_time_ms, ioc=DEFAULT_QT_IOC):
    """
    Return the archiver PV name for one trim quad at one cycle time.
    """

    family = str(family).upper()
    device = f"R{int(superperiod)}{family}"
    return f"{ioc}::{device}:CURRENT:{_format_ms(cycle_time_ms)}MS"


def get_trim_quad_currents(
    cycle_time_ms,
    fetch_value,
    list_available_times,
    as_of=None,
    superperiods=TRIM_QUAD_SUPERPERIODS,
    ioc=DEFAULT_QT_IOC,
):
    """
    Fetch real per-superperiod QTD/QTF currents. Same injected-callable
    pattern as get_corrector_settings, for the same reason.

    Returns (currents, missing): currents is a dict keyed by
    (superperiod, family) -> current_A; missing is a list of
    (superperiod, family) pairs with no available sample at or before
    cycle_time_ms.
    """

    currents = {}
    missing = []

    for superperiod in superperiods:
        for family in ("QTD", "QTF"):
            device = trim_quad_device_name(superperiod, family)
            available = list_available_times(device, family)
            resolved_time = nearest_cycle_time(available, cycle_time_ms)
            if resolved_time is None:
                missing.append((superperiod, family))
                continue

            pv_name = trim_quad_pv_name(superperiod, family, resolved_time, ioc=ioc)
            currents[(superperiod, family)] = fetch_value(pv_name, as_of=as_of)

    return currents, missing


def get_requested_tune(
    cycle_time_ms,
    fetch_value,
    list_available_times,
    beam_state,
    as_of=None,
    superperiods=TRIM_QUAD_SUPERPERIODS,
    base_qx=DEFAULT_BASE_QX,
    base_qy=DEFAULT_BASE_QY,
):
    """
    Derive (set_qx, set_qy) for one A3 timepoint row from real trim-quad
    currents, by averaging the per-superperiod QTD/QTF readings into a
    single pair and reversing the model's own Di Wright equations
    (trim_quad_current_to_tune_di). See the module note above -- the
    averaging step is not confirmed with staff.

    Returns (row, missing): row is a dict {cycle_time_ms, set_qx, set_qy,
    iqtf_A, iqtd_A} ready to feed into snapshot_configs_from_table, or
    None if no trim quads had any data. missing is passed through from
    get_trim_quad_currents so an incomplete average is still visible.
    """

    currents, missing = get_trim_quad_currents(
        cycle_time_ms, fetch_value, list_available_times, as_of=as_of, superperiods=superperiods
    )

    qtd_values = [value for (superperiod, family), value in currents.items() if family == "QTD"]
    qtf_values = [value for (superperiod, family), value in currents.items() if family == "QTF"]

    if not qtd_values or not qtf_values:
        return None, missing

    iqtd_A = sum(qtd_values) / len(qtd_values)
    iqtf_A = sum(qtf_values) / len(qtf_values)

    qx, qy = trim_quad_current_to_tune_di(
        iqtf_A=iqtf_A,
        iqtd_A=iqtd_A,
        base_qx=base_qx,
        base_qy=base_qy,
        pn=float(beam_state.normalised_momentum),
    )

    row = {
        "cycle_time_ms": float(cycle_time_ms),
        "set_qx": qx,
        "set_qy": qy,
        "iqtf_A": iqtf_A,
        "iqtd_A": iqtd_A,
    }
    return row, missing
