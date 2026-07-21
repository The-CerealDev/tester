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
