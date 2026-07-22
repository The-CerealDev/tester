"""
archiver_live_app.py

Tiny Streamlit app wiring the real EPICS archiver client
(optics_gui.io.epics_archiver_client) to the real corrector/BPM/tune
builders (optics_gui.io.epics_live). This is the actual "live" tool --
unlike a hosted artifact, it runs locally and can reach the ISIS network
directly.

Run it once you're on the network:

    streamlit run Dev/12_IO/archiver_live_app.py
"""

import sys
from datetime import date, datetime, time, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import streamlit as st

from optics_gui.io import (
    archiver_fetch_value,
    archiver_list_available_times,
    archiver_list_available_times_dwtrim,
    bpm_geometry_table,
    get_bpm_measurements,
    get_corrector_settings,
    get_harmonic_tunes,
    get_requested_tune,
)
from optics_gui.io.epics_archiver_client import DEFAULT_BASE_URL
from optics_gui.machine_state import MachineState
from optics_gui.cycle_time import RCSRamp
from optics_gui.madx_model import MadxModel

st.set_page_config(page_title="ISIS Archiver Bridge", page_icon="\U0001F50C", layout="wide")

st.title("ISIS Archiver → optics_gui bridge")
st.caption(
    "Runs the real functions built this session: epics_archiver_client.py + epics_live.py, "
    "against the real confirmed endpoints (/glob, /getPVStatus, /data). Local tool, not a "
    "hosted page, so it can actually reach the archiver."
)

with st.sidebar:
    st.header("Archiver connection")
    base_url = st.text_input("Archiver base URL", value=DEFAULT_BASE_URL)
    cycle_time_ms = st.number_input("cycle_time_ms", value=0.0, step=0.5)
    st.caption("cycle_time_ms picks where in the ~10ms accelerator ramp -- every team uses this axis.")

    st.divider()
    st.header("Point in history")
    history_mode = st.radio("When", ["Live (now)", "As of a specific day"], index=0)
    selected_day = None
    if history_mode == "As of a specific day":
        selected_day = st.date_input("Day", value=date.today())
    st.caption(
        "Separate axis from cycle_time_ms: this picks which calendar day's archive to read "
        "from. If that exact day has no data, it automatically searches further back in time "
        "until it finds a real sample -- never invents or interpolates a value."
    )


def _end_of_day_utc(day):
    return datetime.combine(day, time(23, 59, 59, 999999), tzinfo=timezone.utc)


def make_fetch_value():
    if history_mode == "As of a specific day" and selected_day is not None:
        pinned_as_of = _end_of_day_utc(selected_day)

        def fetch_value(pv_name, as_of=None):  # as_of param intentionally ignored: pinned_as_of wins
            return archiver_fetch_value(
                pv_name, as_of=pinned_as_of, base_url=base_url, lookback_days=1, expand_search=True
            )

        return fetch_value

    def fetch_value(pv_name, as_of=None):
        return archiver_fetch_value(pv_name, as_of=as_of, base_url=base_url)

    return fetch_value


def make_list_available_times():
    def list_available_times(device, family):
        return archiver_list_available_times(device, family, base_url=base_url)

    return list_available_times


def make_list_available_times_dwtrim():
    def list_available_times_dwtrim(signal):
        return archiver_list_available_times_dwtrim(signal, base_url=base_url)

    return list_available_times_dwtrim


st.divider()
st.subheader("Corrector settings")
st.caption("Calls get_corrector_settings(...) for all 14 confirmed HD/VD correctors.")

if st.button("Fetch corrector settings", type="primary"):
    try:
        settings, missing = get_corrector_settings(
            cycle_time_ms=cycle_time_ms,
            fetch_value=make_fetch_value(),
            list_available_times=make_list_available_times(),
        )
    except Exception as exc:  # noqa: BLE001 -- surface the real error to the page, don't swallow it
        st.error(f"Fetch failed: {exc}")
    else:
        if settings is None:
            st.warning("No correctors had an available sample at or before this cycle time.")
        else:
            col1, col2 = st.columns(2)
            with col1:
                st.write("**HD currents (A)**")
                st.json(settings.hd_corrector_currents_A)
            with col2:
                st.write("**VD currents (A)**")
                st.json(settings.vd_corrector_currents_A)
        if missing:
            st.warning(f"No data for {len(missing)} corrector(s): {missing}")

st.divider()
st.subheader("Requested tune (A3)")
st.caption(
    "Calls get_requested_tune(...), reading DWTRIM::H_Q / DWTRIM::V_Q directly for this "
    "cycle time -- confirmed live to be the real tune setpoints (H_Q:AT_TIME:0MS == 4.331 == "
    "DEFAULT_BASE_QX, V_Q:AT_TIME:0MS == 3.731 == DEFAULT_BASE_QY). No averaging across the "
    "10 superperiods' trim-quad currents -- each cycle time is read individually."
)

if st.button("Fetch requested tune", type="primary"):
    try:
        row, missing = get_requested_tune(
            cycle_time_ms=cycle_time_ms,
            fetch_value=make_fetch_value(),
            list_available_times_dwtrim=make_list_available_times_dwtrim(),
        )
    except Exception as exc:  # noqa: BLE001
        st.error(f"Fetch failed: {exc}")
    else:
        if row is None:
            st.warning(f"No data at this cycle time for: {missing}")
        else:
            col1, col2 = st.columns(2)
            col1.metric("set_qx (H_Q)", f"{row['set_qx']:.4f}")
            col2.metric("set_qy (V_Q)", f"{row['set_qy']:.4f}")
            if missing:
                st.warning(f"Missing: {missing}")

st.divider()
st.subheader("Harmonic amplitudes")
st.caption(
    "Calls get_harmonic_tunes(...), reading the real DWTRIM::{D7,D8,F8,F9}{SIN,COS} PVs that "
    "MachineState.harmonic_tunes (DEFAULT_HARMONICS) expects. F9SIN/F9COS are confirmed NOT "
    "present on this archiver -- they'll show up under 'missing', not a fabricated 0.0."
)

if st.button("Fetch harmonic amplitudes", type="primary"):
    try:
        values, missing = get_harmonic_tunes(
            cycle_time_ms=cycle_time_ms,
            fetch_value=make_fetch_value(),
            list_available_times_dwtrim=make_list_available_times_dwtrim(),
        )
    except Exception as exc:  # noqa: BLE001
        st.error(f"Fetch failed: {exc}")
    else:
        st.json(dict(values))
        if missing:
            st.warning(f"No data for: {missing}")

st.divider()
st.subheader("BPM measurements")
st.caption(
    "Geometry (bpm/plane/s) comes from a real MAD-X TWISS run -- no EPICS needed for that part. "
    "The PV name pattern below is NOT confirmed by staff yet; edit it once you have the real readback PV names."
)

lattice_folder = st.text_input(
    "Lattice folder", value=str(REPO_ROOT / "Dev" / "Lattice_Files" / "00_Simplified_Lattice")
)
bpm_pv_pattern = st.text_input(
    "BPM PV name pattern (use {bpm} and {plane})",
    value="UNCONFIRMED::{bpm}:{plane}",
)

if st.button("Fetch BPM measurements", type="primary"):
    with st.spinner("Running MAD-X for BPM geometry..."):
        beam_state = RCSRamp().state_at(cycle_time_ms)
        machine_state = MachineState.from_defaults(beam_state=beam_state)
        model = MadxModel(
            lattice_folder=lattice_folder,
            sequence_name="synchrotron",
            aperture_file=None,
            output_dir=str(REPO_ROOT / "Dev" / "12_IO" / "student_runs" / "_live_app"),
        )
        model.load_lattice(use_sequence=False)
        model.apply_machine_state(machine_state)
        model.use_sequence()
        twiss_df = model.run_twiss()
        geometry = bpm_geometry_table(twiss_df)

    st.caption(f"Geometry OK -- {len(geometry)} real BPM positions from the lattice.")

    def pv_name_for_bpm(bpm_label, plane):
        return bpm_pv_pattern.format(bpm=str(bpm_label).upper(), plane=plane)

    try:
        measurements = get_bpm_measurements(geometry, make_fetch_value(), pv_name_for_bpm)
    except Exception as exc:  # noqa: BLE001
        st.error(f"Fetch failed: {exc}")
    else:
        st.dataframe(measurements, use_container_width=True)
