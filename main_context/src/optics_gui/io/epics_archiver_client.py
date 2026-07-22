"""
epics_archiver_client.py

Real HTTP client for the ISIS archiver, against the real endpoints staff
gave directly (not the generic EPICS Archiver Appliance API this file
used to guess at):

    /glob?pv=<pattern>                          -- list PV names (wildcard search)
    /getPVStatus?pv=<pv>                        -- check whether a PV exists/is archived
    /data?pv=<pv>&from=<ISO8601>&to=<ISO8601>   -- fetch archived samples for one PV

Base host confirmed: http://athena.isis.rl.ac.uk:9506

What's actually verified here vs. still unconfirmed:

- /glob's path AND response shape are real, not guessed: TEST_PV.txt (a
  real PV export sitting in this same folder) is in exactly the
  array-of-{pvName, appliance, connectionState, samplingPeriod, lastEvent,
  ...} shape this file parses. archiver_list_pvs() /
  archiver_list_available_times() are tested against that real file.

- /data and /getPVStatus's response SHAPES are NOT yet confirmed -- nobody
  has pasted back a real response. archiver_fetch_value() and
  archiver_get_pv_status() hit the correct real URL and params, but the
  code that reads the JSON back out is a best-effort guess across a few
  plausible shapes, and raises a clear error naming the raw payload if it
  can't find a value it recognises, rather than silently returning a wrong
  number. Replace _extract_latest_value once a real /data response has
  been seen.
"""

from datetime import datetime, timedelta, timezone
import re

import requests

DEFAULT_BASE_URL = "http://athena.isis.rl.ac.uk:9506"

_TIME_SUFFIX_RE = re.compile(r":CURRENT:([-\d.]+)MS$", re.IGNORECASE)


def _iso(moment):
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


def archiver_list_pvs(pattern, base_url=DEFAULT_BASE_URL, timeout=10):
    """
    List PVs matching a glob pattern, e.g. "DWHST_TEST::R0HD1:CURRENT:*".
    Real endpoint and response shape, verified against TEST_PV.txt.
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
    timeout=10,
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


def archiver_get_pv_status(pv_name, base_url=DEFAULT_BASE_URL, timeout=10):
    """
    Check whether a PV exists/is archived, via /getPVStatus.

    NOT YET VERIFIED against a real response (see module docstring) --
    returns the parsed JSON (or raw text if it isn't JSON) as-is, so the
    real shape can be inspected the first time this actually runs against
    the network, instead of pretending to know a schema nobody has seen.
    """

    response = requests.get(f"{base_url}/getPVStatus", params={"pv": pv_name}, timeout=timeout)
    response.raise_for_status()
    try:
        return response.json()
    except ValueError:
        return response.text


def _extract_latest_value(payload):
    """
    Best-effort extraction across a few plausible /data response shapes.
    Not confirmed -- replace once a real response has been observed.
    """

    if isinstance(payload, list) and payload:
        first = payload[0]
        if isinstance(first, dict) and isinstance(first.get("data"), list) and first["data"]:
            last = first["data"][-1]
            if isinstance(last, dict):
                for key in ("val", "value", "y"):
                    if key in last:
                        return float(last[key])
        if isinstance(first, dict):
            for key in ("val", "value"):
                if key in first:
                    return float(first[key])
    if isinstance(payload, dict):
        for key in ("val", "value"):
            if key in payload and isinstance(payload[key], (int, float)):
                return float(payload[key])
        if isinstance(payload.get("data"), list) and payload["data"]:
            last = payload["data"][-1]
            if isinstance(last, dict):
                for key in ("val", "value", "y"):
                    if key in last:
                        return float(last[key])
    return None


def archiver_fetch_value(pv_name, as_of=None, base_url=DEFAULT_BASE_URL, lookback_days=400, timeout=10):
    """
    Fetch one PV's most recent archived value at or before `as_of`
    (default: now), via /data. Real fetch_value callable for
    get_corrector_settings / get_bpm_measurements / get_trim_quad_currents.

    Response-shape parsing is NOT yet confirmed -- see module docstring.
    Raises a clear error naming the raw payload rather than guessing a
    wrong number if the shape doesn't match what _extract_latest_value
    expects.
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
            f"Got a /data response for {pv_name!r} but couldn't find a recognised value "
            f"field in it -- the real response shape still needs confirming. "
            f"Raw payload: {payload!r}"
        )
    return value
