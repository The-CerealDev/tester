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
import streamlit as st

CACHE_PATH = Path(__file__).resolve().parent / "pv_search_cache.json"

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


def normalise(records):
    df = pd.DataFrame(records)
    if df.empty:
        return df
    if "pvName" not in df.columns and "pvNameOnly" in df.columns:
        df["pvName"] = df["pvNameOnly"]
    df["connectionState"] = df["connectionState"].astype(str).str.lower() == "true"
    keep = [c for c in ["pvName", "appliance", "connectionState", "samplingPeriod", "lastEvent"] if c in df.columns]
    return df.loc[:, keep].reset_index(drop=True)


def parse_magnet(pv_name):
    match = MAGNET_RE.search(str(pv_name))
    if not match:
        return pd.Series({"superperiod": None, "family": None, "cycle_time_ms": None})
    sp, family, cycle = match.groups()
    family = family.upper()
    if family in ("HD1", "VD1"):
        family = family[:-1]  # the regex group includes the trailing "1" only for HD1/VD1, not QTD/QTF
    return pd.Series({"superperiod": int(sp), "family": family, "cycle_time_ms": float(cycle)})


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
            save_cache(records)
            st.success(f"Loaded {len(records)} PVs.")
            st.rerun()
    if col_b.button("Load example rows"):
        st.session_state.pv_records = EXAMPLE_ROWS
        save_cache(EXAMPLE_ROWS)
        st.rerun()
    if col_c.button("Clear all data"):
        st.session_state.pv_records = []
        if CACHE_PATH.exists():
            CACHE_PATH.unlink()
        st.rerun()
    st.caption(f"Cached at `{CACHE_PATH.name}` in this folder (gitignored) -- reloads automatically next time you run the app.")

df = normalise(st.session_state.pv_records)

if df.empty:
    st.info("No PVs loaded yet. Open the panel above, paste your getAllPVs JSON, and click Load -- or try the example rows.")
    st.stop()

st.divider()

search_col, regex_col = st.columns([5, 1])
search_text = search_col.text_input(
    "Search PV names",
    value=st.session_state.search_text,
    placeholder=r"try HD1:CURRENT or a regex like R\d(HD1|VD1|QT[DF])",
)
regex_mode = regex_col.checkbox("Regex", value=st.session_state.regex_mode)

st.write("**Quick filters**")
preset_cols = st.columns(len(PRESETS))
for col, (label, pattern) in zip(preset_cols, PRESETS.items()):
    if col.button(label):
        st.session_state.search_text = pattern
        st.session_state.regex_mode = bool(pattern)
        st.rerun()

filter_col1, filter_col2 = st.columns(2)
connected_only = filter_col1.checkbox("Connected only")
exclude_test = filter_col2.checkbox("Exclude _TEST")

mask = pd.Series(True, index=df.index)
if search_text:
    try:
        pattern = search_text if regex_mode else re.escape(search_text)
        mask &= df["pvName"].str.contains(pattern, case=False, regex=True, na=False)
    except re.error as exc:
        st.error(f"Invalid regex: {exc}")
        mask &= False
if connected_only:
    mask &= df["connectionState"]
if exclude_test:
    mask &= ~df["pvName"].str.contains("_TEST", case=False, na=False)

filtered = df.loc[mask].copy()
parsed = filtered["pvName"].apply(parse_magnet)
filtered = pd.concat([filtered, parsed], axis=1)

st.caption(f"**{len(filtered)}** / {len(df)} PVs shown")
st.dataframe(filtered, use_container_width=True, hide_index=True)

if not filtered.empty:
    st.download_button(
        "Download filtered as CSV",
        data=filtered.to_csv(index=False).encode("utf-8"),
        file_name="pv_search_results.csv",
        mime="text/csv",
    )
