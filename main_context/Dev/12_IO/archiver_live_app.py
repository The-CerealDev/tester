"""
archiver_live_app.py

Tiny Streamlit app wiring the real EPICS archiver client
(optics_gui.io.epics_archiver_client) to the real corrector/BPM builders
(optics_gui.io.epics_live). This is the actual "live" tool -- unlike a
hosted artifact, it runs locally and can reach the ISIS network directly.

Run it once you're on the network and have the archiver's real host/port:

    streamlit run Dev/12_IO/archiver_live_app.py

Nothing here fabricates a connection: both archiver URL fields are blank
by default, and every fetch button fails with a clear message if you try
to run it without them, exactly like the underlying Python functions do.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import streamlit as st

from optics_gui.io import (
    archiver_fetch_value,
    archiver_list_available_times,
    bpm_geometry_table,
    get_bpm_measurements,
    get_corrector_settings,
)
from optics_gui.machine_state import MachineState
from optics_gui.cycle_time import RCSRamp
from optics_gui.madx_model import MadxModel

st.set_page_config(page_title="ISIS Archiver Bridge", page_icon="\U0001F50C", layout="wide")

st.title("ISIS Archiver → optics_gui bridge")
st.caption(
    "Runs the real functions built this session: epics_archiver_client.py + epics_live.py. "
    "This is a local tool, not a hosted page, so it can actually reach the archiver."
)

with st.sidebar:
    st.header("Archiver connection")
    mgmt_base_url = st.text_input(
        "Management base URL",
        value="",
        placeholder="http://<archiver-host>:17665/mgmt/bpl",
        help="Used by list_available_times (getAllPVs). Get the real host/port from staff.",
    )
    retrieval_base_url = st.text_input(
        "Retrieval base URL",
        value="",
        placeholder="http://<archiver-host>:17668/retrieval",
        help="Used by fetch_value (getData.json). Get the real host/port from staff.",
    )
    cycle_time_ms = st.number_input("cycle_time_ms", value=0.0, step=0.5)
    st.caption("Leave the URLs blank to see the exact error the real functions raise instead of a fake connection.")


def make_fetch_value():
    def fetch_value(pv_name, as_of=None):
        return archiver_fetch_value(pv_name, as_of=as_of, retrieval_base_url=retrieval_base_url)

    return fetch_value


def make_list_available_times():
    def list_available_times(device, family):
        return archiver_list_available_times(device, family, mgmt_base_url=mgmt_base_url)

    return list_available_times


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
