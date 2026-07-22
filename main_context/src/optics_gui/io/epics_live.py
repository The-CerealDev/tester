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
#
# FOR FUTURE REFERENCE (not wired in -- deprioritized, BPM isn't heavily
# used in the control system yet; revisit when that changes): a live
# /glob search found a much stronger real-BPM candidate --
# RNG:DIAG:POS:R{superperiod}{HM|VM}{n}:{POSITION|IN/OUT|UP/DOWN}, 38
# devices, naming matches the lattice's real BPM labels (sp0_r0hm1 ->
# R0HM1) far better than CHANGE:ORBIT_*, and includes the raw-electrode-
# pair-plus-derived-position structure real BPM hardware has. Confirmed
# live via /getPVStatus that R0HM1:POSITION is real and "Being archived."
# Confirmed live via /data that it holds a genuinely varying value (not a
# flat setpoint like every other test-IOC PV this session) -- a 10-year
# window returned 12,762,602 samples, value moving from 0.0 to a nonzero
# reading. The catch: /data against this PV is very slow regardless of
# window size (a 5-minute window failed at 45s; only a 10-year window
# with a 240s timeout succeeded) -- whoever picks this back up needs to
# either confirm that's an intrinsic property of this densely-sampled PV
# class (needs a long fetch_value timeout, not the ~20s default) or retest
# a moderate window in isolation, since this session had already put
# heavy load on the archiver (real 429s/500s seen elsewhere) by the time
# this was tried, which may have been a confounding factor.
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
# Trim quads -- real per-superperiod QTD/QTF currents, one IOC per
# superperiod pair (same DWQ_TEST convention as the correctors). Still
# useful as raw diagnostic data even though get_requested_tune() below no
# longer derives tune from these (see the DWTRIM note further down).
# ----------------------------------------------------------------------

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


# ----------------------------------------------------------------------
# Requested tune / harmonics -- confirmed live on the real archiver, via
# a real DWTRIM IOC that earlier PV-name searches (TUNE, QX, QY, QH, QV,
# WORKING_POINT, SETPOINT, PROGRAM, RESONANCE, TARGET against the small
# TEST_PV.txt export) had missed entirely because that export didn't
# include DWTRIM at all:
#
#   DWTRIM::H_Q:AT_TIME:<ms>MS   -- real tune setpoint, horizontal.
#       DWTRIM::H_Q:AT_TIME:0MS fetched live = 4.331, exactly
#       DEFAULT_BASE_QX -- confirms this is the real set_qx PV, not a
#       derived/estimated value.
#   DWTRIM::V_Q:AT_TIME:<ms>MS   -- real tune setpoint, vertical.
#       DWTRIM::V_Q:AT_TIME:0MS fetched live = 3.731, exactly
#       DEFAULT_BASE_QY.
#   DWTRIM::{D7,D8,F8}{SIN,COS}:AT_TIME:<ms>MS -- harmonic-correction
#       amplitudes. DEFAULT_HARMONICS also expects F9SIN/F9COS; those two
#       were NOT found on this archiver -- exhaustively checked (every
#       case/separator variant, *F9* anywhere on the whole archiver not
#       just DWTRIM, DWTRIM::* re-confirmed at its real uncapped total of
#       213 PVs) so this isn't a naming miss on our end. tune_plots.py
#       shows why this is worth asking staff about rather than shrugging
#       off: D7/D8 modulate qtd_k and F8/F9 modulate qtf_k as a matched
#       pair of harmonics per family (7th+8th for D, 8th+9th for F) -- F9
#       is structurally as load-bearing as F8, which does exist live, so
#       this looks like this test rig just never had F9 wired up rather
#       than the real machine lacking 9th-harmonic correction entirely.
#       get_harmonic_tunes() reports both as missing rather than inventing
#       a value either way.
#
# So set_qx/set_qy is read directly per cycle time now, not derived by
# averaging the 10 superperiods' QTD/QTF currents and reversing the
# model's tune-control equations -- H_Q/V_Q are each already the single
# real setpoint for that instant, so there is nothing to average.
# ----------------------------------------------------------------------

from collections import OrderedDict  # noqa: E402

DEFAULT_DWTRIM_IOC = "DWTRIM"
HARMONIC_KEYS = ("D7SIN", "D7COS", "D8SIN", "D8COS", "F8SIN", "F8COS", "F9SIN", "F9COS")


def dwtrim_pv_name(signal, cycle_time_ms, ioc=DEFAULT_DWTRIM_IOC):
    """
    Return the archiver PV name for one DWTRIM ":AT_TIME:" signal at one
    cycle time, e.g. dwtrim_pv_name("H_Q", 0) -> "DWTRIM::H_Q:AT_TIME:0MS".
    """

    return f"{ioc}::{signal}:AT_TIME:{_format_ms(cycle_time_ms)}MS"


def get_requested_tune(
    cycle_time_ms,
    fetch_value,
    list_available_times_dwtrim,
    as_of=None,
    ioc=DEFAULT_DWTRIM_IOC,
):
    """
    Read (set_qx, set_qy) for one A3 timepoint row directly from the real
    DWTRIM::H_Q / DWTRIM::V_Q tune-setpoint PVs -- see the module note
    above for the live confirmation. Each is resolved independently to
    its own nearest-available cycle time at/before cycle_time_ms (no
    cross-superperiod averaging, since these are already single PVs).

    list_available_times_dwtrim(signal) -> iterable of float
        Lists the cycle-time-ms suffixes actually archived for one DWTRIM
        AT_TIME signal, e.g. epics_archiver_client.archiver_list_available_times_dwtrim.

    Returns (row, missing): row is a dict {cycle_time_ms, set_qx, set_qy}
    ready to feed into snapshot_configs_from_table, or None if either
    H_Q or V_Q had no available sample at or before cycle_time_ms.
    missing lists which of "H_Q"/"V_Q" (if any) were unavailable.
    """

    missing = []
    resolved = {}

    for signal in ("H_Q", "V_Q"):
        available = list_available_times_dwtrim(signal)
        resolved_time = nearest_cycle_time(available, cycle_time_ms)
        if resolved_time is None:
            missing.append(signal)
            continue
        pv_name = dwtrim_pv_name(signal, resolved_time, ioc=ioc)
        resolved[signal] = float(fetch_value(pv_name, as_of=as_of))

    if "H_Q" not in resolved or "V_Q" not in resolved:
        return None, missing

    row = {
        "cycle_time_ms": float(cycle_time_ms),
        "set_qx": resolved["H_Q"],
        "set_qy": resolved["V_Q"],
    }
    return row, missing


def get_harmonic_tunes(
    cycle_time_ms,
    fetch_value,
    list_available_times_dwtrim,
    as_of=None,
    ioc=DEFAULT_DWTRIM_IOC,
    keys=HARMONIC_KEYS,
):
    """
    Fetch the harmonic-correction amplitudes that MachineState.harmonic_tunes
    (DEFAULT_HARMONICS: D7SIN, D7COS, D8SIN, D8COS, F8SIN, F8COS, F9SIN,
    F9COS) expects, from the real DWTRIM IOC, one cycle time at a time.

    Returns (values, missing): values is an OrderedDict of whichever keys
    had a real sample at/before cycle_time_ms (in DEFAULT_HARMONICS order,
    ready to merge into MachineState.harmonic_tunes / a config's
    "harmonics" column); missing lists the rest -- confirmed live that
    F9SIN/F9COS are not present on this archiver, so callers should expect
    those two to come back missing rather than treat it as an error.
    """

    values = OrderedDict()
    missing = []

    for key in keys:
        available = list_available_times_dwtrim(key)
        resolved_time = nearest_cycle_time(available, cycle_time_ms)
        if resolved_time is None:
            missing.append(key)
            continue
        pv_name = dwtrim_pv_name(key, resolved_time, ioc=ioc)
        values[key] = float(fetch_value(pv_name, as_of=as_of))

    return values, missing


def get_timepoint_row(
    cycle_time_ms,
    fetch_value,
    list_available_times_dwtrim,
    as_of=None,
    ioc=DEFAULT_DWTRIM_IOC,
    snapshot_id=None,
    harmonic_keys=HARMONIC_KEYS,
):
    """
    Build one ready-to-use A3 timepoint_table row -- {cycle_time_ms, set_qx,
    set_qy, snapshot_id, harmonics}, the exact shape the student guide's A3
    example feeds straight into snapshot_configs_from_table. Combines
    get_requested_tune() and get_harmonic_tunes() so a caller building a
    multi-row programme table doesn't have to know both need calling and
    merging by hand -- see Dev/12_IO/archiver_live_app.py or the pipeline
    test this was verified against for the multi-row usage pattern.

    Returns (row, missing): row is None if H_Q/V_Q were unavailable (same
    condition as get_requested_tune); missing is {"tune": [...],
    "harmonics": [...]}, each passed straight through from the two
    underlying calls.
    """

    tune_row, tune_missing = get_requested_tune(
        cycle_time_ms, fetch_value, list_available_times_dwtrim, as_of=as_of, ioc=ioc
    )
    harmonics, harmonics_missing = get_harmonic_tunes(
        cycle_time_ms, fetch_value, list_available_times_dwtrim, as_of=as_of, ioc=ioc, keys=harmonic_keys
    )
    missing = {"tune": tune_missing, "harmonics": harmonics_missing}

    if tune_row is None:
        return None, missing

    row = dict(tune_row)
    if snapshot_id is not None:
        row["snapshot_id"] = snapshot_id
    row["harmonics"] = dict(harmonics)
    return row, missing
