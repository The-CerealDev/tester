"""
Plotting helpers for MAD-X error-table misalignments.

This module keeps error-table plotting separate from the pure table utilities
in errors.py. It consumes MAD-X-style error tables and converts supported ISIS
main-magnet rows into marker-style offset segments versus ring S.
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd

from errors import normalise_error_table_columns, read_error_table


RING_CIRCUMFERENCE_M = 163.36282
DIPOLE_LENGTH_M = 5.18
STRAIGHT_SECTION_LENGTH_M = 11.156282
SUPERPERIOD_COUNT = 10

MAGNET_LENGTHS_M = {
    # These are survey marker-to-marker lengths used for plotting offsets.
    # They deliberately differ from the full magnet design lengths for quads.
    "DIP": DIPOLE_LENGTH_M,
    "QD": 0.739,
    "QF": 0.722,
    "QDS": 0.383,
}
MAGNET_OFFSETS_FROM_D_M = {
    "DIP": -0.5 * DIPOLE_LENGTH_M,
    "QD": 3.401,
    "QF": 4.846,
    "QDS": 10.95,
}
MAGNET_LABELS = {
    "DIP": "Dipole",
    "QD": "QD",
    "QF": "QF",
    "QDS": "QC",
}
MAGNET_COLOURS = {
    "Dipole": "#FFCC00",
    "QD": "red",
    "QF": "blue",
    "QC": "green",
}

SUPPORTED_ERROR_NAME_RE = re.compile(
    r"^SP(?P<period>\d+)_(?P<kind>DIP|QD|QF|QDS)$",
    re.IGNORECASE,
)


def _load_error_table_like(error_table):
    if isinstance(error_table, (str, Path)):
        return read_error_table(error_table)
    if isinstance(error_table, pd.DataFrame):
        return normalise_error_table_columns(error_table, fill_missing=True, copy=True)
    raise TypeError("error_table must be a pandas DataFrame or file path.")


def _normalise_plane(plane):
    plane = str(plane).strip().lower()
    aliases = {
        "x": "x",
        "h": "x",
        "horizontal": "x",
        "y": "y",
        "v": "y",
        "vertical": "y",
    }
    if plane not in aliases:
        raise ValueError("plane must be 'x'/'horizontal' or 'y'/'vertical'.")
    return aliases[plane]


def _error_name_parts(name):
    match = SUPPORTED_ERROR_NAME_RE.match(str(name).strip())
    if match is None:
        return None
    period = int(match.group("period"))
    kind = match.group("kind").upper()
    if period < 0 or period >= SUPERPERIOD_COUNT:
        return None
    return period, kind


def _continuous_s_window(start, centre, end):
    if start < 0.0:
        return (
            start + RING_CIRCUMFERENCE_M,
            centre + RING_CIRCUMFERENCE_M,
            end + RING_CIRCUMFERENCE_M,
        )
    return start, centre, end


def _default_magnet_s(period, kind):
    d_marker_s = period * (DIPOLE_LENGTH_M + STRAIGHT_SECTION_LENGTH_M)
    centre = d_marker_s + MAGNET_OFFSETS_FROM_D_M[kind]
    length = MAGNET_LENGTHS_M[kind]
    start = centre - 0.5 * length
    end = centre + 0.5 * length

    return _continuous_s_window(start, centre, end)


def _clean_lattice_name(name):
    return str(name).strip().strip("\"'").split(":")[0].lower()


def _twiss_magnet_s(name, twiss_df):
    if twiss_df is None:
        return None
    if not isinstance(twiss_df, pd.DataFrame):
        raise TypeError("twiss_df must be a pandas DataFrame.")
    if "name" not in twiss_df.columns or "s" not in twiss_df.columns:
        raise ValueError("twiss_df must contain 'name' and 's' columns.")

    lookup_key = _clean_lattice_name(name)
    local = twiss_df.copy()
    local["_clean_name"] = local["name"].map(_clean_lattice_name)
    rows = local.loc[local["_clean_name"] == lookup_key]
    if rows.empty:
        return None

    centre = float(pd.to_numeric(rows["s"], errors="coerce").dropna().iloc[0])
    length = None
    if "l" in rows.columns:
        lengths = pd.to_numeric(rows["l"], errors="coerce").dropna()
        if not lengths.empty and float(lengths.iloc[0]) > 0.0:
            length = float(lengths.iloc[0])
    if length is None:
        parts = _error_name_parts(name)
        if parts is None:
            return None
        _, kind = parts
        length = MAGNET_LENGTHS_M[kind]
    start = centre - 0.5 * length
    end = centre + 0.5 * length
    return _continuous_s_window(start, centre, end)


def error_table_to_misalignment_offsets(error_table, plane="x", strict=False, twiss_df=None):
    """
    Convert an ISIS MAD-X error table into marker-style misalignment rows.

    Returned offsets are in millimetres and angles are in milliradians, matching
    the plotting convention used by the survey-to-misalignment notebooks.
    """

    df = _load_error_table_like(error_table)
    df = normalise_error_table_columns(df, fill_missing=True, copy=True)
    plane = _normalise_plane(plane)

    rows = []
    skipped = []
    for _, row in df.iterrows():
        name = str(row["name"])
        parts = _error_name_parts(name)
        if parts is None:
            skipped.append(name)
            continue

        period, kind = parts
        s_values = _twiss_magnet_s(name, twiss_df)
        if s_values is None:
            s_values = _default_magnet_s(period, kind)
        s_start, s_centre, s_end = s_values

        if plane == "x":
            offset = 1.0e3 * float(row["dx"])
            angle = 1.0e3 * float(row["dtheta"])
        else:
            offset = 1.0e3 * float(row["dy"])
            angle = -1.0e3 * float(row["dphi"])

        label = f"{MAGNET_LABELS[kind]} {period}"
        rows.append(
            {
                "name": name,
                "magnet": label,
                "magnet_type": MAGNET_LABELS[kind],
                "S_start": s_start,
                "S_centre": s_centre,
                "S_end": s_end,
                "offset_centre": offset,
                "angle": angle,
                "plane": plane,
            }
        )

    if strict and skipped:
        raise ValueError(f"Unsupported error-table element names: {skipped}")

    out = pd.DataFrame(
        rows,
        columns=[
            "name",
            "magnet",
            "magnet_type",
            "S_start",
            "S_centre",
            "S_end",
            "offset_centre",
            "angle",
            "plane",
        ],
    )
    if out.empty:
        raise ValueError("No supported ISIS main-magnet error rows were found.")

    numeric_columns = ("S_start", "S_centre", "S_end", "offset_centre", "angle")
    for column in numeric_columns:
        out[column] = pd.to_numeric(out[column], errors="coerce").astype(float)
    if not np.isfinite(out[list(numeric_columns)].to_numpy()).all():
        raise ValueError("Misalignment offset table contains non-finite values.")

    return out.sort_values("S_centre").reset_index(drop=True)


def _offset_column(df, corrected):
    if corrected and "offset_corrected" in df.columns:
        return "offset_corrected"
    return "offset_centre"


def _symmetric_integer_ylim(values):
    numeric = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna()
    if numeric.empty:
        return -1, 1
    limit = max(1, int(np.ceil(numeric.abs().max())))
    return -limit, limit


def _marker_segment_offsets(df, offset_col):
    values = []
    for _, row in df.iterrows():
        half_length = (row["S_end"] - row["S_start"]) / 2
        offset = row[offset_col]
        angle_rad = row["angle"] / 1000.0
        values.extend(
            [
                offset - half_length * 1000.0 * np.tan(angle_rad),
                offset + half_length * 1000.0 * np.tan(angle_rad),
            ]
        )
    return values


def _plot_offset_centres(ax, df, offset_col):
    ax.scatter(df["S_centre"], df[offset_col], s=30)
    for _, row in df.iterrows():
        ax.text(row["S_centre"], row[offset_col], str(row["magnet"]), fontsize=6, rotation=45)


def _plot_magnet_segments(ax, df, offset_col):
    for _, row in df.iterrows():
        magnet = str(row["magnet"])
        magnet_type = magnet.split()[0]
        colour = MAGNET_COLOURS.get(magnet_type, "black")

        half_length = (row["S_end"] - row["S_start"]) / 2
        offset = row[offset_col]
        angle_rad = row["angle"] / 1000.0

        s_values = [row["S_start"], row["S_end"]]
        offset_values = [
            offset - half_length * 1000.0 * np.tan(angle_rad),
            offset + half_length * 1000.0 * np.tan(angle_rad),
        ]
        ax.plot(s_values, offset_values, color=colour, marker="o", linewidth=1.4)
        ax.text(row["S_centre"], offset + 0.35, magnet, ha="center", va="bottom", fontsize=7, color=colour)


def plot_misalignment_offsets(
    df,
    *,
    corrected=False,
    show_segments=True,
    mode=None,
    title="Misalignment Offsets",
    savename=None,
):
    """
    Plot magnet misalignments against S.

    The input is the DataFrame returned by error_table_to_misalignment_offsets()
    or an equivalent survey-style table with S_start/S_end/S_centre,
    offset_centre and angle columns. Angles are expected in milliradians.
    """

    import matplotlib.pyplot as plt

    if mode is not None:
        if mode not in {"segments", "offsets"}:
            raise ValueError("mode must be 'segments' or 'offsets'.")
        show_segments = mode == "segments"

    required = {"magnet", "S_start", "S_end", "S_centre", "offset_centre", "angle"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Misalignment DataFrame is missing columns: {missing}")

    offset_col = _offset_column(df, corrected)
    fig, ax = plt.subplots(figsize=(12, 5))
    if show_segments:
        _plot_magnet_segments(ax, df, offset_col)
        ax.set_ylim(*_symmetric_integer_ylim(_marker_segment_offsets(df, offset_col)))
        legend_handles = [
            plt.Line2D([0], [0], color=colour, marker="o", linestyle="-", label=label)
            for label, colour in MAGNET_COLOURS.items()
        ]
        ax.legend(handles=legend_handles, loc="best", fontsize=8, title="Magnet")
    else:
        _plot_offset_centres(ax, df, offset_col)
        ax.set_ylim(*_symmetric_integer_ylim(df[offset_col]))

    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlim(0.0, RING_CIRCUMFERENCE_M)
    ax.set_xlabel("S [m]")
    ax.set_ylabel("Offset [mm]")
    ax.set_title(title)
    ax.grid(True, linestyle=":", linewidth=0.5)
    if savename is not None:
        Path(savename).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savename, bbox_inches="tight", dpi=200)
    return fig, ax


def plot_error_table_misalignment_offsets(
    error_table,
    *,
    plane="x",
    strict=False,
    twiss_df=None,
    title=None,
    savename=None,
    show_segments=True,
    mode=None,
):
    """
    Plot MAD-X error-table misalignment offsets versus ring S.
    """

    plane = _normalise_plane(plane)
    offsets = error_table_to_misalignment_offsets(
        error_table,
        plane=plane,
        strict=strict,
        twiss_df=twiss_df,
    )
    if title is None:
        title = f"{'Horizontal' if plane == 'x' else 'Vertical'} Misalignment Offsets"
    return plot_misalignment_offsets(
        offsets,
        title=title,
        savename=savename,
        show_segments=show_segments,
        mode=mode,
    )
