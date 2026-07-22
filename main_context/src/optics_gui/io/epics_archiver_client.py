"""
epics_archiver_client.py

Real HTTP client for the ISIS archiver, against the real endpoints staff
gave directly:

    /glob?pv=<pattern>                          -- list PV names (wildcard search)
    /getPVStatus?pv=<pv>                        -- check whether a PV exists/is archived
    /data?pv=<pv>&from=<ISO8601>&to=<ISO8601>   -- fetch archived samples for one PV

Base host confirmed: http://athena.isis.rl.ac.uk:9506

Every response shape below was observed live against the real archiver
from this codebase (not guessed, and not just checked against a static
export):

- /glob returns a flat JSON array of plain PV name strings, e.g.
  ["AC:BEAM:ENERGY", "DWHST_TEST::R0HD1:CURRENT:0MS", ...]. A bare
  pv="*" is capped by the server at 500 results and does not include the
  correctors/trim-quads at all -- use a real prefix/suffix pattern
  (e.g. "DWHST_TEST::*", "*HD1:CURRENT*") to get useful results, the
  same way the worked examples used a real prefix, not a bare "*".

- /getPVStatus for one specific PV returns a JSON array containing one
  object: [{"appliance": ..., "connectionState": "true"/"false",
  "pvName": ..., "pvNameOnly": ..., "samplingPeriod": ...,
  "lastEvent": ..., ...}] -- the same per-PV record shape as a getAllPVs
  dump. getPVStatus?pv=* was tested and hangs/times out -- it is not
  built for bulk queries, /glob is.

- /data returns one JSON object with three parallel columnar dicts,
  each keyed by stringified row index: {"secs": {"0": ..., "1": ...},
  "nanos": {...}, "val": {...}}. The latest sample is the entry with the
  largest (secs, nanos) pair, not necessarily the highest-numbered index
  (verified in order for the one PV tested, but sorted explicitly here
  rather than assumed).

There is also a real DWTRIM IOC (confirmed live via /glob, separate from
the per-superperiod DWQ_TEST trim-quad IOC) exposing single ":AT_TIME:
<ms>MS" signals instead of ":CURRENT:<ms>MS" -- notably DWTRIM::H_Q and
DWTRIM::V_Q (the real tune setpoints: DWTRIM::H_Q:AT_TIME:0MS returned
4.331 live, matching DEFAULT_BASE_QX exactly; DWTRIM::V_Q:AT_TIME:0MS
returned 3.731, matching DEFAULT_BASE_QY), and the harmonic-correction
amplitudes DWTRIM::D7SIN/D7COS/D8SIN/D8COS/F8SIN/F8COS (F9SIN/F9COS,
also expected by DEFAULT_HARMONICS, were not found on this archiver).
"""

from datetime import datetime, timedelta, timezone
import re

import requests

DEFAULT_BASE_URL = "http://athena.isis.rl.ac.uk:9506"

_TIME_SUFFIX_RE = re.compile(r":CURRENT:([-\d.]+)MS$", re.IGNORECASE)
_AT_TIME_SUFFIX_RE = re.compile(r":AT_TIME:([-\d.]+)MS$", re.IGNORECASE)


def _iso(moment):
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


def archiver_list_pvs(pattern, base_url=DEFAULT_BASE_URL, timeout=20):
    """
    List PVs matching a glob pattern, e.g. "DWHST_TEST::R0HD1:CURRENT:*".

    Real endpoint and response shape, observed live: a flat JSON array
    of plain name strings. Avoid a bare "*" -- the server caps that at
    500 results and, tested live, that capped set does not include the
    correctors/trim-quads at all. Use a real prefix/suffix pattern.
    """

    response = requests.get(f"{base_url}/glob", params={"pv": pattern}, timeout=timeout)
    response.raise_for_status()
    entries = response.json()

    names = []
    for entry in entries:
        if isinstance(entry, str):
            names.append(entry)
        else:
            names.append(entry.get("pvName") or entry.get("pvNameOnly") or "")
    return [name for name in names if name]


def archiver_list_available_times(
    device,
    family,
    hd_ioc="DWHST_TEST",
    vd_ioc="DWVST_TEST",
    qtd_qtf_ioc="DWQ_TEST",
    base_url=DEFAULT_BASE_URL,
    timeout=20,
):
    """
    List the cycle-time-ms suffixes actually archived for one magnet,
    real list_available_times callable for get_corrector_settings /
    get_trim_quad_currents.
    """

    family_upper = str(family).upper()
    if family_upper == "HD":
        ioc = hd_ioc
    elif family_upper == "VD":
        ioc = vd_ioc
    else:
        ioc = qtd_qtf_ioc  # QTD / QTF

    pattern = f"{ioc}::{str(device).upper()}:CURRENT:*"
    names = archiver_list_pvs(pattern, base_url=base_url, timeout=timeout)

    times = []
    for name in names:
        match = _TIME_SUFFIX_RE.search(name)
        if match:
            times.append(float(match.group(1)))
    return times


def archiver_list_available_times_dwtrim(signal, ioc="DWTRIM", base_url=DEFAULT_BASE_URL, timeout=20):
    """
    List the cycle-time-ms suffixes actually archived for one DWTRIM
    ":AT_TIME:<ms>MS" signal, e.g. "H_Q", "V_Q", "D7SIN", "F8COS".

    Same nearest-cycle-time-fallback role as archiver_list_available_times,
    but DWTRIM signals are single global PVs (not one per superperiod) and
    use ":AT_TIME:" rather than ":CURRENT:" as the suffix marker -- real
    endpoint/shape confirmed live, see the module docstring.
    """

    pattern = f"{ioc}::{signal}:AT_TIME:*"
    names = archiver_list_pvs(pattern, base_url=base_url, timeout=timeout)

    times = []
    for name in names:
        match = _AT_TIME_SUFFIX_RE.search(name)
        if match:
            times.append(float(match.group(1)))
    return times


def archiver_get_pv_status(pv_name, base_url=DEFAULT_BASE_URL, timeout=10):
    """
    Check whether a PV exists/is archived, via /getPVStatus.

    Real response shape, observed live: a JSON array containing one
    object (the same per-PV record shape as a getAllPVs/glob-style
    dump) -- e.g. [{"pvName": ..., "connectionState": "true", ...}].
    Do not pass pv="*" here -- tested live, that hangs/times out;
    /glob is the bulk-listing endpoint, this one is single-PV only.
    """

    response = requests.get(f"{base_url}/getPVStatus", params={"pv": pv_name}, timeout=timeout)
    response.raise_for_status()
    return response.json()


def _extract_latest_value(payload):
    """
    Parse the real /data response shape, observed live against the ISIS
    archiver: one JSON object with three parallel columnar dicts, each
    keyed by stringified row index -- {"secs": {...}, "nanos": {...},
    "val": {...}}. Returns the val at the latest (secs, nanos) pair, or
    None if the payload doesn't match this shape at all (e.g. no samples
    for that PV/time range).
    """

    if not isinstance(payload, dict):
        return None

    secs = payload.get("secs")
    vals = payload.get("val")
    if not isinstance(secs, dict) or not isinstance(vals, dict) or not secs:
        return None

    nanos = payload.get("nanos") or {}
    latest_index = max(secs, key=lambda i: (secs[i], nanos.get(i, 0)))
    return float(vals[latest_index])


def archiver_fetch_value(pv_name, as_of=None, base_url=DEFAULT_BASE_URL, lookback_days=400, timeout=20):
    """
    Fetch one PV's most recent archived value at or before `as_of`
    (default: now), via /data. Real fetch_value callable for
    get_corrector_settings / get_bpm_measurements / get_trim_quad_currents.

    Response-shape parsing is confirmed against a real live response
    (see _extract_latest_value). Still raises a clear error naming the
    raw payload if a given PV/time-range genuinely has no samples,
    rather than returning something misleading.
    """

    to_time = as_of if as_of is not None else datetime.now(timezone.utc)
    from_time = to_time - timedelta(days=lookback_days)

    response = requests.get(
        f"{base_url}/data",
        params={"pv": pv_name, "from": _iso(from_time), "to": _iso(to_time)},
        timeout=timeout,
    )
    response.raise_for_status()
    payload = response.json()

    value = _extract_latest_value(payload)
    if value is None:
        raise ValueError(
            f"No archived samples found for {pv_name!r} in the {lookback_days}-day window "
            f"before {to_time.isoformat()}. Raw payload: {payload!r}"
        )
    return value
