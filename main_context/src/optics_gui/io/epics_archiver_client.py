"""
epics_archiver_client.py

Real HTTP client for the EPICS Archiver Appliance, matching the two
standard endpoints:

    mgmt/bpl/getAllPVs        -- list PV names (glob pattern search)
    retrieval/data/getData.json -- fetch archived/live samples for one PV

The getAllPVs response shape here (pvName, pvNameOnly, appliance,
connectionState, samplingPeriod, lastEvent, ...) matches the archiver dump
already confirmed real against this project's lattice, so this client is
written against the real API, not a guess. The one missing piece is the
archiver's own host/port, which only staff can supply -- every function
here takes that as a required argument and raises a clear error if it is
left out, rather than silently pointing at a made-up address.

These functions are meant to be passed as the fetch_value /
list_available_times callables that epics_live.py's get_corrector_settings
and get_bpm_measurements already accept -- nothing in epics_live.py needs
to change to use them.
"""

from datetime import datetime, timedelta, timezone
import re

import requests

_TIME_SUFFIX_RE = re.compile(r":CURRENT:([-\d.]+)MS$", re.IGNORECASE)


def _iso(moment):
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


def archiver_fetch_value(pv_name, as_of=None, retrieval_base_url=None, lookback_days=400, timeout=10):
    """
    Fetch one PV's most recent archived value at or before `as_of`
    (default: now). This is the real fetch_value callable for
    get_corrector_settings / get_bpm_measurements.

    retrieval_base_url must point at the archiver's retrieval service,
    e.g. "http://<archiver-host>:17668/retrieval". Ask staff for the real
    host/port -- do not guess one.
    """

    if not retrieval_base_url:
        raise ValueError(
            "retrieval_base_url is required (e.g. 'http://<archiver-host>:17668/retrieval'). "
            "Get the real host/port from staff -- this is not something to guess."
        )

    to_time = as_of if as_of is not None else datetime.now(timezone.utc)
    from_time = to_time - timedelta(days=lookback_days)

    response = requests.get(
        f"{retrieval_base_url}/data/getData.json",
        params={"pv": pv_name, "from": _iso(from_time), "to": _iso(to_time)},
        timeout=timeout,
    )
    response.raise_for_status()
    payload = response.json()

    if not payload or not payload[0].get("data"):
        raise ValueError(
            f"No archived samples for {pv_name!r} in the {lookback_days}-day window "
            f"before {to_time.isoformat()}."
        )

    latest_sample = payload[0]["data"][-1]
    return float(latest_sample["val"])


def archiver_list_available_times(
    device,
    family,
    hd_ioc="DWHST_TEST",
    vd_ioc="DWVST_TEST",
    mgmt_base_url=None,
    timeout=10,
):
    """
    List the cycle-time-ms suffixes actually archived for one magnet, by
    querying getAllPVs with a glob pattern and parsing the :CURRENT:<ms>MS
    suffix out of each matching PV name. This is the real
    list_available_times callable for get_corrector_settings.

    mgmt_base_url must point at the archiver's management service, e.g.
    "http://<archiver-host>:17665/mgmt/bpl". Ask staff for the real
    host/port -- do not guess one.
    """

    if not mgmt_base_url:
        raise ValueError(
            "mgmt_base_url is required (e.g. 'http://<archiver-host>:17665/mgmt/bpl'). "
            "Get the real host/port from staff -- this is not something to guess."
        )

    ioc = hd_ioc if str(family).upper() == "HD" else vd_ioc
    pattern = f"{ioc}::{str(device).upper()}:CURRENT:*"

    response = requests.get(f"{mgmt_base_url}/getAllPVs", params={"pv": pattern}, timeout=timeout)
    response.raise_for_status()
    entries = response.json()

    times = []
    for entry in entries:
        name = entry if isinstance(entry, str) else entry.get("pvName", "")
        match = _TIME_SUFFIX_RE.search(name)
        if match:
            times.append(float(match.group(1)))
    return times
