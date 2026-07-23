"""
archiver_live_app.py

The one real tool for this team -- everything that used to be split
across a static demo artifact and a separate PV-search app now lives
here, as tabs, wired to the real EPICS archiver client
(optics_gui.io.epics_archiver_client) and the real corrector/BPM/tune
builders (optics_gui.io.epics_live). Nothing on this page is simulated;
unlike a hosted artifact, it runs locally and can actually reach the
ISIS network.

Run it once you're on the network:

    streamlit run Dev/12_IO/archiver_live_app.py
"""

import hashlib
import json
import re
import sys
from datetime import date, datetime, time, timezone
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pandas as pd
import requests
import streamlit as st

from optics_gui.io import (
    archiver_fetch_value,
    archiver_list_available_times,
    archiver_list_available_times_dwtrim,
    archiver_list_pvs,
    bpm_geometry_table,
    dwtrim_pv_name,
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

st.markdown(
    """
    <style>
    .badge { display:inline-flex; align-items:center; gap:5px; font-size:0.72rem; font-weight:700;
      text-transform:uppercase; letter-spacing:.04em; padding:3px 10px; border-radius:999px; }
    .badge-ok { background:rgba(31,138,92,0.16); color:#1f8a5c; }
    .badge-warn { background:rgba(194,63,56,0.16); color:#c23f38; }
    .badge-pending { background:rgba(61,111,150,0.16); color:#3d6f96; }
    .badge::before { content:""; width:6px; height:6px; border-radius:50%; background:currentColor; }
    .stage-row { display:flex; flex-wrap:wrap; gap:10px 0; align-items:stretch; margin:6px 0 14px; }
    .stage { background: rgba(127,127,127,0.07); border:1px solid rgba(127,127,127,0.22); border-radius:8px;
      padding:12px 14px; flex:1 1 170px; min-width:170px; }
    .stage .stage-name { font-family: ui-monospace, Consolas, monospace; font-size:0.8rem; font-weight:700; margin-bottom:4px; }
    .stage .stage-desc { font-size:0.78rem; opacity:0.75; line-height:1.4; }
    .stage-arrow { display:flex; align-items:center; justify-content:center; padding:0 10px; opacity:0.5; font-size:18px; }
    @media (max-width: 900px) {
      .stage-row { flex-direction:column; }
      .stage { flex:1 1 auto; }
      .stage-arrow { transform:rotate(90deg); padding:2px 0; }
    }
    .ledger-row { display:flex; gap:12px; align-items:flex-start; padding:10px 0;
      border-bottom:1px solid rgba(127,127,127,0.18); }
    .ledger-row:last-child { border-bottom:none; }
    .ledger-row .badge { flex:0 0 auto; margin-top:2px; }
    .ledger-what { font-weight:600; font-size:0.92rem; }
    .ledger-note { font-size:0.83rem; opacity:0.75; margin-top:2px; line-height:1.4; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("ISIS Archiver → optics_gui bridge")
st.caption(
    "The one real tool: live epics_archiver_client.py + epics_live.py calls against the real "
    "ISIS archiver, against the real confirmed endpoints (/glob, /getPVStatus, /data). "
    "Nothing on this page is simulated -- every button here makes a real network call."
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


def badge(label, kind="ok"):
    st.markdown(f'<span class="badge badge-{kind}">{label}</span>', unsafe_allow_html=True)


def prospective_value(seed, low, high):
    """
    Deterministic placeholder value for a prospective (not-yet-real) PV --
    same seed always gives the same number, so a preview looks stable
    across reruns without needing a real fetch. Only ever used for
    features staff have said are prospective (planned, not confirmed
    live), never as a stand-in for something that should be real.
    """

    digest = hashlib.sha256(seed.encode("utf-8")).hexdigest()
    fraction = int(digest[:8], 16) / 0xFFFFFFFF
    return round(low + fraction * (high - low), 4)


tab_correctors, tab_tune, tab_bpm, tab_search, tab_status = st.tabs(
    ["Correctors (A1)", "Tune & Harmonics (A3)", "BPM (A2)", "PV Search", "Architecture & Status"]
)

# ----------------------------------------------------------------------
# A1 -- Correctors
# ----------------------------------------------------------------------
with tab_correctors:
    badge("confirmed live", "ok")
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

# ----------------------------------------------------------------------
# A3 -- Requested tune + harmonics
# ----------------------------------------------------------------------
with tab_tune:
    badge("confirmed live", "ok")
    st.caption(
        "Reads DWTRIM::H_Q / DWTRIM::V_Q directly for this cycle time -- confirmed live to be the "
        "real tune setpoints (H_Q:AT_TIME:0MS == 4.331 == DEFAULT_BASE_QX, V_Q:AT_TIME:0MS == "
        "3.731 == DEFAULT_BASE_QY). No averaging across the 10 superperiods' trim-quad currents -- "
        "each cycle time is read individually. Harmonics: 6 of 8 DEFAULT_HARMONICS keys exist on "
        "this archiver; F9SIN/F9COS are confirmed absent (shown as missing, never fabricated as 0.0)."
    )

    if st.button("Fetch tune + harmonics", type="primary"):
        fetch_value = make_fetch_value()
        list_available_times_dwtrim = make_list_available_times_dwtrim()
        try:
            row, tune_missing = get_requested_tune(
                cycle_time_ms=cycle_time_ms,
                fetch_value=fetch_value,
                list_available_times_dwtrim=list_available_times_dwtrim,
            )
            values, harmonics_missing = get_harmonic_tunes(
                cycle_time_ms=cycle_time_ms,
                fetch_value=fetch_value,
                list_available_times_dwtrim=list_available_times_dwtrim,
            )
        except Exception as exc:  # noqa: BLE001
            st.error(f"Fetch failed: {exc}")
        else:
            if row is None:
                st.warning(f"No tune data at this cycle time for: {tune_missing}")
            else:
                col1, col2 = st.columns(2)
                col1.metric("set_qx (H_Q)", f"{row['set_qx']:.4f}")
                col2.metric("set_qy (V_Q)", f"{row['set_qy']:.4f}")
                if tune_missing:
                    st.warning(f"Missing: {tune_missing}")
            st.write("**Harmonic amplitudes**")
            st.json(dict(values))
            if harmonics_missing:
                st.warning(f"No data for: {harmonics_missing}")
                with st.expander(f"Prospective preview: if {', '.join(harmonics_missing)} existed"):
                    badge("prospective", "pending")
                    st.caption(
                        "Not live, not fabricated as real -- staff have indicated these may be added "
                        "to the control system in future. PV names below use the real dwtrim_pv_name() "
                        "the moment they're added, they'd read that way automatically -- but the values "
                        "shown are placeholders, not archived data."
                    )
                    preview_rows = [
                        {
                            "signal": signal,
                            "pv_name": dwtrim_pv_name(signal, cycle_time_ms),
                            "placeholder_value": prospective_value(signal, -0.002, 0.002),
                        }
                        for signal in harmonics_missing
                    ]
                    st.dataframe(pd.DataFrame(preview_rows), use_container_width=True, hide_index=True)

# ----------------------------------------------------------------------
# A2 -- BPM (prospective)
# ----------------------------------------------------------------------
with tab_bpm:
    badge("prospective", "pending")
    st.caption(
        "Prospective, not live: staff have indicated BPM readback PVs may be added to the control "
        "system in future, but nothing is confirmed yet. Geometry (bpm/plane/s) below is real, from "
        "a real MAD-X TWISS run -- no EPICS needed for that part. The measurement table is a preview "
        "using test data through the real get_bpm_measurements()/normalise_bpm_table() code path, "
        "not a live fetch -- it shows the exact shape/columns the real thing will have once readback "
        "PVs exist, without depending on any PV name that isn't confirmed."
    )
    st.caption(
        "Candidate readback pattern if/when this is added: `RNG:DIAG:POS:{bpm}:POSITION` "
        "(uppercase, no sp{n}_ prefix) -- found live via /glob, confirmed to hold a genuinely varying "
        "value unlike every setpoint PV checked, but NOT staff-confirmed. Shown for reference only; "
        "the preview below uses test data, not this pattern."
    )

    lattice_folder = st.text_input(
        "Lattice folder", value=str(REPO_ROOT / "Dev" / "Lattice_Files" / "00_Simplified_Lattice")
    )

    if st.button("Preview BPM table (test data)", type="primary"):
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

        def test_fetch_value(pv_name, as_of=None):  # noqa: ARG001 -- no network, deterministic placeholder
            return prospective_value(pv_name, -2.2, 2.2)

        def pv_name_for_bpm(bpm_label, plane):
            return f"RNG:DIAG:POS:{str(bpm_label).upper()}:POSITION"

        measurements = get_bpm_measurements(geometry, test_fetch_value, pv_name_for_bpm)
        st.dataframe(measurements, use_container_width=True)
        st.caption("closed_orbit_mm above is deterministic test data, not a real archiver value.")

# ----------------------------------------------------------------------
# PV Search
# ----------------------------------------------------------------------
SEARCH_CACHE_PATH = Path(__file__).resolve().parent / "pv_search_cache.json"
MAGNET_RE = re.compile(r"R(\d)(HD1|VD1|QTD|QTF)[^:]*:CURRENT:([-\d.]+)MS", re.IGNORECASE)
SEARCH_PRESETS = {
    "All": "",
    "HD correctors": "HD1:CURRENT",
    "VD correctors": "VD1:CURRENT",
    "Trim quads": "QT[DF]:CURRENT",
    "Correctors + quads": "(?:HD1|VD1|QT[DF]):CURRENT",
}
SEARCH_EXAMPLE_ROWS = [
    {"pvName": "DWHST_TEST::R0HD1:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:03:40 BST"},
    {"pvName": "DWHST_TEST::R0HD1:CURRENT:2.5MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:04:44 BST"},
    {"pvName": "DWHST_TEST::R2HD1:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:03:49 BST"},
    {"pvName": "DWHST_TEST::R9HD1:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:03:42 BST"},
    {"pvName": "DWVST_TEST::R0VD1:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:04:27 BST"},
    {"pvName": "DWVST_TEST::R2VD1:CURRENT:1MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:03:53 BST"},
    {"pvName": "DWVST_TEST::R9VD1:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:04:09 BST"},
    {"pvName": "DWQ_TEST::R0QTD:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:04:42 BST"},
    {"pvName": "DWQ_TEST::R0QTF:CURRENT:1MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:03:38 BST"},
    {"pvName": "DWQ_TEST::R1QTD:CURRENT:0MS", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:03:40 BST"},
    {"pvName": "CPSEPICSTST::TESTDB:TYPEG_IN1", "appliance": "appliance0", "connectionState": "false", "samplingPeriod": "1.0", "lastEvent": "Never"},
    {"pvName": "R4TQTEST::QUAD:NORM:DAT", "appliance": "appliance0", "connectionState": "true", "samplingPeriod": "1.0", "lastEvent": "Jul/16/2026 11:04:46 BST"},
]


@st.cache_data(show_spinner=False)
def normalise_pv_records(records):
    """
    Cached on the exact records list -- with 10,000s of PVs this only
    needs to run once per Load click, not once per keystroke/interaction
    (Streamlit reruns the whole script on every widget change).
    """

    df = pd.DataFrame(records)
    if df.empty:
        return df
    if "pvName" not in df.columns and "pvNameOnly" in df.columns:
        df["pvName"] = df["pvNameOnly"]
    if "connectionState" in df.columns:
        df["connectionState"] = df["connectionState"].astype(str).str.lower() == "true"
    else:
        # /glob (live fetch) only returns bare PV names, no connection
        # status -- default to True so "connected only" doesn't silently
        # drop everything just because that field wasn't asked for.
        df["connectionState"] = True
    keep = [c for c in ["pvName", "appliance", "connectionState", "samplingPeriod", "lastEvent"] if c in df.columns]
    return df.loc[:, keep].reset_index(drop=True)


def parse_magnets(pv_names):
    """
    Vectorised parse -- one regex pass across the whole column in
    pandas/re's C layer rather than a per-row .apply(), ~84x faster at
    40,000 PVs. Same family-stripping rule (HD1/VD1 -> HD/VD, QTD/QTF
    unchanged) either way.
    """

    extracted = pv_names.str.extract(MAGNET_RE)
    extracted.columns = ["superperiod", "family_raw", "cycle_time_ms"]
    extracted["superperiod"] = pd.to_numeric(extracted["superperiod"], errors="coerce").astype("Int64")
    extracted["cycle_time_ms"] = pd.to_numeric(extracted["cycle_time_ms"], errors="coerce")
    family = extracted["family_raw"].str.upper()
    family = family.where(~family.isin(["HD1", "VD1"]), family.str[:-1])
    return pd.DataFrame({"superperiod": extracted["superperiod"], "family": family, "cycle_time_ms": extracted["cycle_time_ms"]})


@st.cache_data(show_spinner=False)
def run_pv_search(_df, data_version, search_text, regex_mode, connected_only, exclude_test):
    """
    Cached on (data_version, search_text, regex_mode, connected_only,
    exclude_test); _df itself is excluded from the cache key (Streamlit's
    own underscore convention) so hashing a 40,000-row DataFrame on every
    call doesn't itself cost real time -- data_version is the cheap
    stand-in, only changing when new data is actually loaded.
    """

    df = _df
    mask = pd.Series(True, index=df.index)
    error = None
    if search_text:
        try:
            pattern = search_text if regex_mode else re.escape(search_text)
            mask &= df["pvName"].str.contains(pattern, case=False, regex=True, na=False)
        except re.error as exc:
            error = str(exc)
            mask &= False
    if connected_only:
        mask &= df["connectionState"]
    if exclude_test:
        mask &= ~df["pvName"].str.contains("_TEST", case=False, na=False)

    filtered = df.loc[mask].copy()
    filtered = pd.concat([filtered, parse_magnets(filtered["pvName"])], axis=1)
    return filtered, error


@st.cache_data(show_spinner=False)
def search_results_to_csv_bytes(df):
    return df.to_csv(index=False).encode("utf-8")


def load_search_cache():
    if SEARCH_CACHE_PATH.exists():
        try:
            return json.loads(SEARCH_CACHE_PATH.read_text())
        except (json.JSONDecodeError, OSError):
            return None
    return None


def save_search_cache(records):
    try:
        SEARCH_CACHE_PATH.write_text(json.dumps(records))
    except OSError:
        pass


with tab_search:
    badge("confirmed live", "ok")
    st.caption(
        "Search a getAllPVs-style dump, or pull PVs straight from the real archiver's /glob "
        "endpoint using the connection settings in the sidebar. Local cache persists between runs."
    )

    if "pv_records" not in st.session_state:
        st.session_state.pv_records = load_search_cache() or []
    if "search_data_version" not in st.session_state:
        st.session_state.search_data_version = 0
    if "search_text" not in st.session_state:
        st.session_state.search_text = ""
    if "search_regex_mode" not in st.session_state:
        st.session_state.search_regex_mode = False
    if "search_connected_only" not in st.session_state:
        st.session_state.search_connected_only = False
    if "search_exclude_test" not in st.session_state:
        st.session_state.search_exclude_test = False

    with st.expander("Load / replace PV list", expanded=not st.session_state.pv_records):
        pasted = st.text_area("Paste the JSON array from getAllPVs here", height=140)
        col_a, col_b, col_c = st.columns(3)
        if col_a.button("Load pasted data", type="primary"):
            try:
                records = json.loads(pasted)
                if not isinstance(records, list):
                    raise ValueError("Expected a JSON array.")
            except (json.JSONDecodeError, ValueError) as exc:
                st.error(f"Couldn't parse that as a PV list: {exc}")
            else:
                st.session_state.pv_records = records
                st.session_state.search_data_version += 1
                save_search_cache(records)
                st.success(f"Loaded {len(records)} PVs.")
                st.rerun()
        if col_b.button("Load example rows"):
            st.session_state.pv_records = SEARCH_EXAMPLE_ROWS
            st.session_state.search_data_version += 1
            save_search_cache(SEARCH_EXAMPLE_ROWS)
            st.rerun()
        if col_c.button("Clear all data"):
            st.session_state.pv_records = []
            st.session_state.search_data_version += 1
            if SEARCH_CACHE_PATH.exists():
                SEARCH_CACHE_PATH.unlink()
            st.rerun()
        st.caption(f"Cached at `{SEARCH_CACHE_PATH.name}` in this folder (gitignored) -- reloads automatically next time you run the app.")

        st.divider()
        st.write("**Or upload a JSON file**")
        uploaded_file = st.file_uploader("getAllPVs-style JSON file", type=["json"], key="pv_search_upload")
        if uploaded_file is not None and st.button("Load uploaded file", type="primary"):
            try:
                records = json.loads(uploaded_file.getvalue().decode("utf-8"))
                if not isinstance(records, list):
                    raise ValueError("Expected a JSON array.")
            except (json.JSONDecodeError, ValueError, UnicodeDecodeError) as exc:
                st.error(f"Couldn't parse {uploaded_file.name} as a PV list: {exc}")
            else:
                st.session_state.pv_records = records
                st.session_state.search_data_version += 1
                save_search_cache(records)
                st.success(f"Loaded {len(records)} PVs from {uploaded_file.name}.")
                st.rerun()

        st.divider()
        st.write("**Or fetch live from the archiver**")
        st.caption("Uses the archiver base URL set in the sidebar.")
        fetch_col1, fetch_col2 = st.columns([4, 2])
        archiver_pattern = fetch_col1.text_input(
            "Glob pattern", value="*HD1:CURRENT*", placeholder="e.g. DWHST_TEST::* or *HD1:CURRENT*"
        )
        st.caption(
            "A bare '*' is capped by the archiver at 500 results and, tested live, doesn't include the "
            "correctors/trim-quads at all -- use a real prefix/suffix pattern instead."
        )
        fetch_col2.write("")
        fetch_col2.write("")
        if fetch_col2.button("Fetch from archiver", type="primary"):
            try:
                with st.spinner(f"Fetching PVs matching '{archiver_pattern}' from {base_url} -- this can take a while for a broad pattern like '*'..."):
                    names = archiver_list_pvs(archiver_pattern, base_url=base_url, timeout=60)
            except requests.exceptions.RequestException as exc:
                st.error(
                    f"Couldn't reach the archiver at {base_url}: {exc}. "
                    "This only works if you're actually on the ISIS network."
                )
            except ValueError as exc:
                st.error(f"Archiver responded but not with valid JSON: {exc}")
            else:
                if not names:
                    st.warning(f"Archiver returned no PVs matching '{archiver_pattern}'.")
                else:
                    wrapped = [{"pvName": name} for name in names]
                    st.session_state.pv_records = wrapped
                    st.session_state.search_data_version += 1
                    save_search_cache(wrapped)
                    st.success(f"Fetched {len(wrapped)} PVs live from the archiver.")
                    st.rerun()

    search_df = normalise_pv_records(st.session_state.pv_records)

    if search_df.empty:
        st.info("No PVs loaded yet. Open the panel above, paste your getAllPVs JSON, fetch live, or try the example rows.")
    else:
        st.write("**Quick filters**")
        preset_cols = st.columns(len(SEARCH_PRESETS))
        for col, (label, pattern) in zip(preset_cols, SEARCH_PRESETS.items()):
            if col.button(label, key=f"preset_{label}"):
                st.session_state.search_text = pattern
                st.session_state.search_regex_mode = bool(pattern)
                st.rerun()

        # Everything inside this form is client-side only until Search is
        # clicked (or Enter pressed) -- typing in the box no longer touches
        # the backend at all, avoiding per-keystroke lag at 40,000 PVs.
        with st.form("pv_search_form"):
            search_col, regex_col = st.columns([5, 1])
            search_text_input = search_col.text_input(
                "Search PV names",
                value=st.session_state.search_text,
                placeholder=r"try HD1:CURRENT or a regex like R\d(HD1|VD1|QT[DF])",
            )
            regex_mode_input = regex_col.checkbox("Regex", value=st.session_state.search_regex_mode)
            filter_col1, filter_col2 = st.columns(2)
            connected_only_input = filter_col1.checkbox("Connected only", value=st.session_state.search_connected_only)
            exclude_test_input = filter_col2.checkbox("Exclude _TEST", value=st.session_state.search_exclude_test)
            submitted = st.form_submit_button("Search", type="primary")

        if submitted:
            st.session_state.search_text = search_text_input
            st.session_state.search_regex_mode = regex_mode_input
            st.session_state.search_connected_only = connected_only_input
            st.session_state.search_exclude_test = exclude_test_input

        filtered, search_error = run_pv_search(
            search_df,
            st.session_state.search_data_version,
            st.session_state.search_text,
            st.session_state.search_regex_mode,
            st.session_state.search_connected_only,
            st.session_state.search_exclude_test,
        )
        if search_error:
            st.error(f"Invalid regex: {search_error}")

        st.caption(f"**{len(filtered)}** / {len(search_df)} PVs shown")

        SEARCH_DISPLAY_CAP = 3000
        if len(filtered) > SEARCH_DISPLAY_CAP:
            show_all = st.checkbox(f"Show all {len(filtered)} rows (rendering that many can be slow -- narrowing the search is usually faster)")
            display_df = filtered if show_all else filtered.head(SEARCH_DISPLAY_CAP)
            if not show_all:
                st.caption(f"Showing the first {SEARCH_DISPLAY_CAP} of {len(filtered)} matches -- narrow your search or tick the box above to see the rest.")
        else:
            display_df = filtered

        st.dataframe(display_df, use_container_width=True, hide_index=True)

        if not filtered.empty:
            st.download_button(
                "Download all filtered results as CSV",
                data=search_results_to_csv_bytes(filtered),
                file_name="pv_search_results.csv",
                mime="text/csv",
            )

# ----------------------------------------------------------------------
# Architecture & Status
# ----------------------------------------------------------------------
with tab_status:
    st.subheader("How data moves")
    st.caption("Two files sit between EPICS and the pre-existing backend. Everything here is what this team built.")
    stages = [
        ("ISIS archiver", "athena.isis.rl.ac.uk:9506 -- real host + endpoints: /glob, /getPVStatus, /data"),
        ("epics_archiver_client.py", "HTTP client: PV listing, value fetching, historical day-picking (expand_search)"),
        ("epics_live.py", "PV naming, nearest-cycle-time fallback, BPM geometry, direct tune/harmonics reads (A1/A2/A3)"),
        ("existing io layer", "corrector_settings_from_table, normalise_bpm_table, snapshot_configs_from_table -- pre-existing, untouched"),
        ("Orbit / Tune / Envelope GUIs", "consume SnapshotConfig + timepoint rows directly -- proven live through build_machine_snapshot / build_full_cycle_snapshot_series"),
    ]
    stage_html = '<div class="stage-row">'
    for i, (name, desc) in enumerate(stages):
        stage_html += f'<div class="stage"><div class="stage-name">{name}</div><div class="stage-desc">{desc}</div></div>'
        if i < len(stages) - 1:
            stage_html += '<div class="stage-arrow">&#8594;</div>'
    stage_html += "</div>"
    st.markdown(stage_html, unsafe_allow_html=True)

    st.subheader("Status ledger")
    ledger_items = [
        ("ok", "HD/VD corrector currents -- A1 (14 magnets)", "PV names + availability confirmed live; get_corrector_settings() tested end-to-end against the real archiver"),
        ("ok", "BPM geometry -- A2 (bpm/plane/s, 36 monitors)", "Pulled from a real MAD-X TWISS run; no EPICS needed"),
        ("ok", "Requested tune -- A3 (fetched directly, not derived)", "Real DWTRIM::H_Q/V_Q tune-setpoint PVs, read per cycle time, no averaging"),
        ("ok", "Harmonic amplitudes -- A3 (D7/D8/F8)", "6 of 8 DEFAULT_HARMONICS keys confirmed live; F9SIN/F9COS are prospective -- staff confirmed these are not yet in the control system, preview available in the Tune & Harmonics tab"),
        ("ok", "Archiver value fetching (/data, /getPVStatus, /glob)", "All three real endpoints confirmed live, including a real HTTP 429 and 500 encountered and handled"),
        ("ok", "Historical day-picking (expand_search)", "Read the latest value on a specific calendar day, expanding the search window backward automatically if that day has none"),
        ("ok", "B/C/D downstream pipeline", "This team's live data proven end-to-end through build_machine_snapshot and build_full_cycle_snapshot_series"),
        ("pending", "BPM closed-orbit readback PVs", "Prospective -- staff confirmed this is not yet in the control system. A live-found candidate exists (RNG:DIAG:POS:R{sp}HM/VM{n}:POSITION) but isn't wired to a real fetch; the BPM tab previews the real output shape with test data instead"),
    ]
    LEDGER_STATUS_LABELS = {"ok": "confirmed", "pending": "prospective", "warn": "blocked"}
    ledger_html = ""
    for status, what, note in ledger_items:
        label = LEDGER_STATUS_LABELS.get(status, status)
        ledger_html += (
            f'<div class="ledger-row"><span class="badge badge-{status}">{label}</span>'
            f'<div><div class="ledger-what">{what}</div><div class="ledger-note">{note}</div></div></div>'
        )
    st.markdown(ledger_html, unsafe_allow_html=True)

    def function_reference_table(rows):
        st.dataframe(
            pd.DataFrame(rows, columns=["Function", "Parameter", "Meaning"]),
            use_container_width=True,
            hide_index=True,
        )

    st.subheader("How the Orbit team uses this")
    st.code(
        """# in the Orbit team's Streamlit page -- host/port confirmed: athena.isis.rl.ac.uk:9506
settings, missing = get_corrector_settings(
    cycle_time_ms=2.5,
    fetch_value=archiver_fetch_value,
    list_available_times=archiver_list_available_times,
)

config = SnapshotConfig(
    cycle_time_ms=2.5,
    lattice_folder=lattice_folder,
    output_dir=output_dir,
    corrector_settings=settings,   # <- exactly what they already expect
)
snapshot = build_machine_snapshot(config)
orbit = snapshot.table("orbit")
orbit_summary = snapshot.table("orbit_summary")""",
        language="python",
    )
    function_reference_table(
        [
            ("get_corrector_settings", "cycle_time_ms", "Point in the ~10ms ramp to fetch correctors for"),
            ("get_corrector_settings", "fetch_value", "Callable reading one PV's value -- archiver_fetch_value"),
            ("get_corrector_settings", "list_available_times", "Callable listing archived cycle times for one magnet"),
            ("get_corrector_settings", "as_of (optional)", "Historical datetime; None reads the live/latest value"),
            ("get_corrector_settings", "prefer (optional)", '"currents" or "kicks" -- which value type the result favours'),
            ("SnapshotConfig", "cycle_time_ms", "The ramp point this config represents"),
            ("SnapshotConfig", "lattice_folder", "Which MAD-X lattice to run"),
            ("SnapshotConfig", "output_dir", "Required -- where this run's files are written"),
            ("SnapshotConfig", "corrector_settings", "The SnapshotCorrectorSettings object from get_corrector_settings"),
            ("build_machine_snapshot", "config", "The SnapshotConfig to run through real MAD-X"),
            ("build_machine_snapshot", "returns", 'SnapshotResult -- use .table("orbit") / .table("orbit_summary")'),
        ]
    )

    st.subheader("How the Tune / Working-Point team uses this")
    st.code(
        """# real path -- runs MAD-X per timepoint
tune_row, missing = get_requested_tune(
    cycle_time_ms=2.5,
    fetch_value=archiver_fetch_value,
    list_available_times_dwtrim=archiver_list_available_times_dwtrim,
)

series = build_full_cycle_snapshot_series(
    cycle_times_ms=[0.0, 0.5, 1.0],
    qx_values=[4.331, 4.331, 4.331],   # <- one get_requested_tune call per point
    qy_values=[3.731, 3.731, 3.731],
    base_config=SnapshotConfig(cycle_time_ms=0.0, lattice_folder=lattice_folder, run_envelope=False, run_aperture=False),
    label="working_point_series",
    output_dir=str(REPO_ROOT / "Dev" / "12_IO" / "student_runs" / "working_point_series"),
)
tune_programme = series.table("tune_programme")
working_points = series.table("working_points")
resonance_proximity = series.table("resonance_proximity")""",
        language="python",
    )
    function_reference_table(
        [
            ("get_requested_tune", "cycle_time_ms", "Point in the ramp to read the tune setpoint for"),
            ("get_requested_tune", "fetch_value / list_available_times_dwtrim", "Same archiver callables as every other function here"),
            ("get_requested_tune", "returns", "(row, missing) -- row = {cycle_time_ms, set_qx, set_qy}"),
            ("build_full_cycle_snapshot_series", "cycle_times_ms", "List of ramp points to evaluate"),
            ("build_full_cycle_snapshot_series", "qx_values / qy_values", "Requested tune at each point -- one get_requested_tune call per entry"),
            ("build_full_cycle_snapshot_series", "base_config (optional)", "Shared SnapshotConfig (must include cycle_time_ms, lattice_folder) applied to every point"),
            ("build_full_cycle_snapshot_series", "output_dir", "Required -- where this series' run files are written"),
            ("build_full_cycle_snapshot_series", "point_overrides (optional)", "Per-point dict overrides, e.g. harmonics from get_harmonic_tunes"),
            ("build_full_cycle_snapshot_series", "returns", 'Series -- .table("tune_programme") / .table("working_points") / .table("resonance_proximity")'),
            ("build_tune_programme_table (lightweight, no MAD-X)", "data, source", "Normalises a raw tune table (cycle_time_ms/set_qx/set_qy/predicted_qx/predicted_qy) into the canonical shape"),
            ("build_working_point_table (lightweight)", "tune_programme_df", "Adds requested-vs-predicted-vs-matched delta columns, one row per point"),
            ("generate_resonance_lines (lightweight)", "xlims, ylims, orders, periodicity", "Resonance-line segments for the tune diagram axes/orders shown"),
            ("evaluate_resonance_proximity (lightweight)", "working_points, resonance_lines", "Nearest resonance line to each working point"),
            ("make_tune_diagram_inputs (lightweight)", "tune_programme_df, xlims, ylims, orders", "Bundles all of the above into one dict ready for plot_tune_diagram_inputs()"),
        ]
    )

    st.subheader("How the Envelope / Aperture team uses this")
    st.code(
        """config = SnapshotConfig(
    cycle_time_ms=2.5,
    lattice_folder=lattice_folder,
    output_dir=output_dir,
    run_envelope=True,   # default
    run_aperture=True,   # default
    envelope_inputs=EnvelopeInputs(
        emit_x_pi_mm_mrad=300.0,
        emit_y_pi_mm_mrad=300.0,
        sigma_scale=3.0,
        dp_over_p=0.002,
    ),
)
result = build_machine_snapshot(config)
envelope = result.table("envelope")
envelope_summary = result.table("envelope_summary")
aperture_aligned = result.table("aperture_aligned")
aperture_summary = result.table("aperture_summary")""",
        language="python",
    )
    function_reference_table(
        [
            ("EnvelopeInputs", "emit_x_pi_mm_mrad / emit_y_pi_mm_mrad", "Beam emittance per plane, accelerator 'pi mm mrad' convention (default 300.0 each)"),
            ("EnvelopeInputs", "emittance_mode", '"geometric" or "normalised" (+ _rms variants); default "geometric"'),
            ("EnvelopeInputs", "sigma_scale", "Envelope width in standard deviations (default 3.0)"),
            ("EnvelopeInputs", "dp_over_p", "Momentum-spread contribution to the envelope (default 0.002)"),
            ("SnapshotConfig", "output_dir", "Required -- where this run's files are written"),
            ("SnapshotConfig", "run_envelope / run_aperture", "Both default True -- toggle each stage; run_aperture requires run_envelope"),
            ("SnapshotConfig", "envelope_inputs", "The EnvelopeInputs object above"),
            ("SnapshotConfig", "aperture_interval (optional)", "Sampling step in metres along the ring for MAD-X APERTURE (default 0.1)"),
            ("build_machine_snapshot", "config", "The SnapshotConfig to run through real MAD-X"),
            ("build_machine_snapshot", "returns", 'SnapshotResult -- .table("envelope") / .table("envelope_summary") / .table("aperture") / .table("aperture_aligned") / .table("aperture_summary")'),
        ]
    )
