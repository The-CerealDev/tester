"""
pv_search_app.py

Local Streamlit port of the PV Search artifact -- same regex-friendly
search over a getAllPVs dump, same quick-filter presets, same
superperiod/family/cycle-time parsing. This version persists what you
load to a local cache file (pv_search_cache.json, gitignored) instead of
browser localStorage, so it survives between runs on this machine.

Run it with:

    streamlit run Dev/12_IO/pv_search_app.py
"""

import json
import re
from pathlib import Path

import pandas as pd
import requests
import streamlit as st

CACHE_PATH = Path(__file__).resolve().parent / "pv_search_cache.json"
DEFAULT_ARCHIVER_BASE_URL = "http://athena.isis.rl.ac.uk:9506"

MAGNET_RE = re.compile(r"R(\d)(HD1|VD1|QTD|QTF)[^:]*:CURRENT:([-\d.]+)MS", re.IGNORECASE)

PRESETS = {
    "All": "",
    "HD correctors": "HD1:CURRENT",
    "VD correctors": "VD1:CURRENT",
    "Trim quads": "QT[DF]:CURRENT",
    "Correctors + quads": "(?:HD1|VD1|QT[DF]):CURRENT",
}

EXAMPLE_ROWS = [
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
def normalise(records):
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
    Vectorised replacement for a per-row .apply(parse_magnet) -- at 40,000
    PVs the old per-row version took ~5 seconds per keystroke (one Python
    function call + one pandas Series built per row); this version runs
    the regex once across the whole column in pandas/re's C layer, ~0.06s
    for the same data. Same output columns, same family-stripping rule
    (HD1/VD1 -> HD/VD, QTD/QTF unchanged) as the row-wise version had.
    """

    extracted = pv_names.str.extract(MAGNET_RE)
    extracted.columns = ["superperiod", "family_raw", "cycle_time_ms"]
    extracted["superperiod"] = pd.to_numeric(extracted["superperiod"], errors="coerce").astype("Int64")
    extracted["cycle_time_ms"] = pd.to_numeric(extracted["cycle_time_ms"], errors="coerce")
    family = extracted["family_raw"].str.upper()
    family = family.where(~family.isin(["HD1", "VD1"]), family.str[:-1])
    return pd.DataFrame({"superperiod": extracted["superperiod"], "family": family, "cycle_time_ms": extracted["cycle_time_ms"]})


@st.cache_data(show_spinner=False)
def run_search(_df, data_version, search_text, regex_mode, connected_only, exclude_test):
    """
    Cached on (data_version, search_text, regex_mode, connected_only,
    exclude_test) -- reruns triggered by something unrelated to the
    search itself (e.g. toggling the "show all" checkbox below) reuse
    this instantly instead of redoing the mask/parse/copy work.

    _df is deliberately excluded from the cache key (the leading
    underscore is Streamlit's own convention for this) -- hashing a
    40,000-row DataFrame on every call just to check the cache would
    itself cost real time. data_version is the cheap stand-in: it only
    changes when new data is actually loaded.
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
def to_csv_bytes(df):
    return df.to_csv(index=False).encode("utf-8")


def fetch_pvs_from_archiver(pattern, base_url, timeout=60):
    """
    Pull PV records live from the real ISIS archiver's /glob endpoint
    (not /getPVStatus -- that's for checking one specific PV, /glob is
    the actual "find all matching PVs" search). Plain Python requests,
    so there's no CORS/CSP involved -- those are browser-only concerns
    that would block this same call from a webpage. This just needs real
    network access to the archiver host, e.g. being on the ISIS network.
    """

    response = requests.get(f"{base_url}/glob", params={"pv": pattern}, timeout=timeout)
    response.raise_for_status()
    return response.json()


def load_cache():
    if CACHE_PATH.exists():
        try:
            return json.loads(CACHE_PATH.read_text())
        except (json.JSONDecodeError, OSError):
            return None
    return None


def save_cache(records):
    try:
        CACHE_PATH.write_text(json.dumps(records))
    except OSError:
        pass


st.set_page_config(page_title="PV Search", page_icon="\U0001F50D", layout="wide")

if "pv_records" not in st.session_state:
    st.session_state.pv_records = load_cache() or []
if "data_version" not in st.session_state:
    st.session_state.data_version = 0
if "search_text" not in st.session_state:
    st.session_state.search_text = ""
if "regex_mode" not in st.session_state:
    st.session_state.regex_mode = False

st.title("PV Search")
st.caption("Local port of the PV Search artifact -- same regex search, same quick filters, persists to a local cache file instead of the browser.")

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
            st.session_state.data_version += 1
            save_cache(records)
            st.success(f"Loaded {len(records)} PVs.")
            st.rerun()
    if col_b.button("Load example rows"):
        st.session_state.pv_records = EXAMPLE_ROWS
        st.session_state.data_version += 1
        save_cache(EXAMPLE_ROWS)
        st.rerun()
    if col_c.button("Clear all data"):
        st.session_state.pv_records = []
        st.session_state.data_version += 1
        if CACHE_PATH.exists():
            CACHE_PATH.unlink()
        st.rerun()
    st.caption(f"Cached at `{CACHE_PATH.name}` in this folder (gitignored) -- reloads automatically next time you run the app.")

    st.divider()
    st.write("**Or fetch live from the archiver**")
    fetch_col1, fetch_col2, fetch_col3 = st.columns([3, 3, 2])
    archiver_pattern = fetch_col1.text_input(
        "Glob pattern", value="*HD1:CURRENT*", placeholder="e.g. DWHST_TEST::* or *HD1:CURRENT*"
    )
    st.caption(
        "A bare '*' is capped by the archiver at 500 results and, tested live, doesn't include the "
        "correctors/trim-quads at all -- use a real prefix/suffix pattern instead."
    )
    archiver_base_url = fetch_col2.text_input("Archiver base URL", value=DEFAULT_ARCHIVER_BASE_URL)
    fetch_col3.write("")
    fetch_col3.write("")
    if fetch_col3.button("Fetch from archiver", type="primary"):
        try:
            with st.spinner(f"Fetching PVs matching '{archiver_pattern}' from {archiver_base_url} -- this can take a while for a broad pattern like '*'..."):
                records = fetch_pvs_from_archiver(archiver_pattern, base_url=archiver_base_url)
        except requests.exceptions.RequestException as exc:
            st.error(
                f"Couldn't reach the archiver at {archiver_base_url}: {exc}. "
                "This only works if you're actually on the ISIS network -- it won't work from an unrelated machine or a hosted page."
            )
        except ValueError as exc:
            st.error(f"Archiver responded but not with valid JSON: {exc}")
        else:
            if not records:
                st.warning(f"Archiver returned no PVs matching '{archiver_pattern}'.")
            else:
                # /glob returns a flat array of plain PV name strings (confirmed live),
                # not the {"pvName": ...} records normalise() expects -- wrap them.
                wrapped = [{"pvName": name} for name in records]
                st.session_state.pv_records = wrapped
                st.session_state.data_version += 1
                save_cache(wrapped)
                st.success(f"Fetched {len(wrapped)} PVs live from the archiver.")
                st.rerun()
    st.caption("This is a real HTTP call to /glob, not a browser fetch -- CORS/CSP restrictions that would block this from a webpage don't apply here.")

df = normalise(st.session_state.pv_records)

if df.empty:
    st.info("No PVs loaded yet. Open the panel above, paste your getAllPVs JSON, and click Load -- or try the example rows.")
    st.stop()

st.divider()

if "connected_only" not in st.session_state:
    st.session_state.connected_only = False
if "exclude_test" not in st.session_state:
    st.session_state.exclude_test = False

st.write("**Quick filters**")
preset_cols = st.columns(len(PRESETS))
for col, (label, pattern) in zip(preset_cols, PRESETS.items()):
    if col.button(label):
        st.session_state.search_text = pattern
        st.session_state.regex_mode = bool(pattern)
        st.rerun()

# Everything inside this form is client-side only until Search is clicked
# (or Enter is pressed) -- typing in the box below no longer touches the
# backend at all, which is what was causing the per-keystroke lag at
# 40,000 PVs even after the search itself got fast.
with st.form("search_form"):
    search_col, regex_col = st.columns([5, 1])
    search_text_input = search_col.text_input(
        "Search PV names",
        value=st.session_state.search_text,
        placeholder=r"try HD1:CURRENT or a regex like R\d(HD1|VD1|QT[DF])",
    )
    regex_mode_input = regex_col.checkbox("Regex", value=st.session_state.regex_mode)
    filter_col1, filter_col2 = st.columns(2)
    connected_only_input = filter_col1.checkbox("Connected only", value=st.session_state.connected_only)
    exclude_test_input = filter_col2.checkbox("Exclude _TEST", value=st.session_state.exclude_test)
    submitted = st.form_submit_button("Search", type="primary")

if submitted:
    st.session_state.search_text = search_text_input
    st.session_state.regex_mode = regex_mode_input
    st.session_state.connected_only = connected_only_input
    st.session_state.exclude_test = exclude_test_input

search_text = st.session_state.search_text
regex_mode = st.session_state.regex_mode
connected_only = st.session_state.connected_only
exclude_test = st.session_state.exclude_test

filtered, search_error = run_search(df, st.session_state.data_version, search_text, regex_mode, connected_only, exclude_test)
if search_error:
    st.error(f"Invalid regex: {search_error}")

st.caption(f"**{len(filtered)}** / {len(df)} PVs shown")

DISPLAY_CAP = 3000
if len(filtered) > DISPLAY_CAP:
    show_all = st.checkbox(f"Show all {len(filtered)} rows (rendering that many can be slow -- narrowing the search is usually faster)")
    display_df = filtered if show_all else filtered.head(DISPLAY_CAP)
    if not show_all:
        st.caption(f"Showing the first {DISPLAY_CAP} of {len(filtered)} matches -- narrow your search or tick the box above to see the rest.")
else:
    display_df = filtered

st.dataframe(display_df, use_container_width=True, hide_index=True)

if not filtered.empty:
    st.download_button(
        "Download all filtered results as CSV",
        data=to_csv_bytes(filtered),
        file_name="pv_search_results.csv",
        mime="text/csv",
    )
