"""
tune.py

GUI-ready tune analysis helpers for the ISIS RCS optics backend.

This layer turns MAD-X summaries and existing tune workflow outputs into
DataFrames suitable for tune diagrams. It does not own low-level MAD-X
execution; workflows that need real model tunes should use MadxModel or the
existing tune_matching wrappers.
"""

import math

import numpy as np
import pandas as pd


TUNE_SUMMARY_ALIASES = {
    "qx": ("qx", "q1", "predicted_qx", "actual_qx", "matched_qx", "set_qx", "requested_qx"),
    "qy": ("qy", "q2", "predicted_qy", "actual_qy", "matched_qy", "set_qy", "requested_qy"),
    "dqx": ("dqx", "dq1", "predicted_dqx", "actual_dqx", "matched_dqx"),
    "dqy": ("dqy", "dq2", "predicted_dqy", "actual_dqy", "matched_dqy"),
}

def _first_present(mapping, names, default=None):
    for name in names:
        if name in mapping and pd.notna(mapping[name]):
            return float(mapping[name])
    return default


def _as_mapping(summary):
    if isinstance(summary, pd.DataFrame):
        if summary.empty:
            raise ValueError("summary DataFrame is empty.")
        return {
            str(key).lower(): value
            for key, value in summary.iloc[0].to_dict().items()
        }

    if hasattr(summary, "summary_df") and summary.summary_df is not None:
        return _as_mapping(summary.summary_df)

    if hasattr(summary, "to_summary_dataframe"):
        return _as_mapping(summary.to_summary_dataframe())

    if isinstance(summary, dict):
        return {str(key).lower(): value for key, value in summary.items()}

    raise TypeError("summary must be a dict, DataFrame, or tune workflow result.")


def extract_tune_summary(summary, source=None):
    """
    Extract a compact tune summary from MAD-X or workflow summary data.
    """

    mapping = _as_mapping(summary)
    return {
        "qx": _first_present(mapping, TUNE_SUMMARY_ALIASES["qx"]),
        "qy": _first_present(mapping, TUNE_SUMMARY_ALIASES["qy"]),
        "dqx": _first_present(mapping, TUNE_SUMMARY_ALIASES["dqx"]),
        "dqy": _first_present(mapping, TUNE_SUMMARY_ALIASES["dqy"]),
        "source": source,
    }


def build_tune_programme_table(data, source="programme"):
    """
    Normalise tune programme outputs into a stable GUI-facing table.
    """

    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    if data.empty:
        raise ValueError("Tune programme data is empty.")

    rows = []
    for _, row in data.iterrows():
        mapping = {str(key).lower(): value for key, value in row.to_dict().items()}
        rows.append(
            {
                "cycle_time_ms": _first_present(mapping, ("cycle_time_ms", "time_ms", "t_ms")),
                "set_qx": _first_present(mapping, ("set_qx", "requested_qx", "di_wright_qx", "qx")),
                "set_qy": _first_present(mapping, ("set_qy", "requested_qy", "di_wright_qy", "qy")),
                "predicted_qx": _first_present(mapping, ("predicted_qx", "actual_qx", "q1")),
                "predicted_qy": _first_present(mapping, ("predicted_qy", "actual_qy", "q2")),
                "matched_qx": _first_present(mapping, ("matched_qx",)),
                "matched_qy": _first_present(mapping, ("matched_qy",)),
                "dqx": _first_present(mapping, ("predicted_dqx", "actual_dqx", "matched_dqx", "dqx", "dq1")),
                "dqy": _first_present(mapping, ("predicted_dqy", "actual_dqy", "matched_dqy", "dqy", "dq2")),
                "iqtf_A": _first_present(mapping, ("iqtf_a",)),
                "iqtd_A": _first_present(mapping, ("iqtd_a",)),
                "kqtf": _first_present(mapping, ("kqtf",)),
                "kqtd": _first_present(mapping, ("kqtd",)),
                "source": mapping.get("source", source),
            }
        )

    programme = pd.DataFrame(rows)
    programme["actual_qx"] = programme["predicted_qx"]
    programme["actual_qy"] = programme["predicted_qy"]
    return programme


def build_working_point_table(
    tune_programme_df,
    use_actual=True,
    use_matched=True,
):
    """
    Build one row per tune point with requested, actual and matched deltas.
    """

    programme = build_tune_programme_table(tune_programme_df)
    rows = []

    for index, row in programme.iterrows():
        set_qx = row["set_qx"]
        set_qy = row["set_qy"]
        predicted_qx = row["predicted_qx"] if use_actual else None
        predicted_qy = row["predicted_qy"] if use_actual else None
        matched_qx = row["matched_qx"] if use_matched else None
        matched_qy = row["matched_qy"] if use_matched else None

        rows.append(
            {
                "index": int(index),
                "cycle_time_ms": row["cycle_time_ms"],
                "set_qx": set_qx,
                "set_qy": set_qy,
                "predicted_qx": predicted_qx,
                "predicted_qy": predicted_qy,
                "actual_qx": predicted_qx,
                "actual_qy": predicted_qy,
                "matched_qx": matched_qx,
                "matched_qy": matched_qy,
                "predicted_minus_set_qx": None if pd.isna(predicted_qx) or pd.isna(set_qx) else predicted_qx - set_qx,
                "predicted_minus_set_qy": None if pd.isna(predicted_qy) or pd.isna(set_qy) else predicted_qy - set_qy,
                "actual_minus_set_qx": None if pd.isna(predicted_qx) or pd.isna(set_qx) else predicted_qx - set_qx,
                "actual_minus_set_qy": None if pd.isna(predicted_qy) or pd.isna(set_qy) else predicted_qy - set_qy,
                "matched_minus_set_qx": None if pd.isna(matched_qx) or pd.isna(set_qx) else matched_qx - set_qx,
                "matched_minus_set_qy": None if pd.isna(matched_qy) or pd.isna(set_qy) else matched_qy - set_qy,
                "source": row["source"],
            }
        )

    return pd.DataFrame(rows)


def _line_rectangle_segment(nx, ny, resonance_sum, xlims, ylims, tolerance=1.0e-12):
    x_min, x_max = map(float, xlims)
    y_min, y_max = map(float, ylims)
    points = []

    if ny != 0:
        for x in (x_min, x_max):
            y = (float(resonance_sum) - float(nx) * x) / float(ny)
            if y_min - tolerance <= y <= y_max + tolerance:
                points.append((x, min(max(y, y_min), y_max)))

    if nx != 0:
        for y in (y_min, y_max):
            x = (float(resonance_sum) - float(ny) * y) / float(nx)
            if x_min - tolerance <= x <= x_max + tolerance:
                points.append((min(max(x, x_min), x_max), y))

    unique = []
    for point in points:
        if not any(abs(point[0] - other[0]) < tolerance and abs(point[1] - other[1]) < tolerance for other in unique):
            unique.append(point)

    if len(unique) < 2:
        return None

    return unique[0], unique[1]


def generate_resonance_lines(
    xlims=(4.0, 4.5),
    ylims=(3.5, 4.0),
    orders=(1, 2, 3, 4),
    periodicity=10,
):
    """
    Generate visible resonance-line segments for a tune diagram.
    """

    x_min, x_max = map(float, xlims)
    y_min, y_max = map(float, ylims)
    rows = []

    for order in tuple(int(order) for order in orders):
        terms = np.arange(-order, order + 1)
        nx_values = order - np.abs(terms)
        ny_values = terms

        for nx, ny in zip(nx_values, ny_values):
            nx = int(nx)
            ny = int(ny)
            if nx == 0 and ny == 0:
                continue

            corners = np.array(
                [
                    nx * math.floor(x_min) + ny * math.floor(y_min),
                    nx * math.ceil(x_max) + ny * math.floor(y_min),
                    nx * math.floor(x_min) + ny * math.ceil(y_max),
                    nx * math.ceil(x_max) + ny * math.ceil(y_max),
                ],
                dtype=int,
            )

            for resonance_sum in range(int(corners.min()), int(corners.max()) + 1):
                segment = _line_rectangle_segment(nx, ny, resonance_sum, (x_min, x_max), (y_min, y_max))
                if segment is None:
                    continue

                (x0, y0), (x1, y1) = segment
                systematic = resonance_sum % int(periodicity) == 0
                skew = bool(ny % 2)
                sign = "+" if ny >= 0 else "-"
                label = f"{nx}Qx {sign} {abs(ny)}Qy = {resonance_sum}"

                rows.append(
                    {
                        "order": order,
                        "nx": nx,
                        "ny": ny,
                        "resonance_sum": int(resonance_sum),
                        "systematic": systematic,
                        "skew": skew,
                        "x0": float(x0),
                        "y0": float(y0),
                        "x1": float(x1),
                        "y1": float(y1),
                        "label": label,
                    }
                )

    return pd.DataFrame(rows)


def _point_line_distance(qx, qy, nx, ny, resonance_sum):
    norm = math.hypot(float(nx), float(ny))
    if norm == 0.0:
        return None, None
    signed = (float(nx) * float(qx) + float(ny) * float(qy) - float(resonance_sum)) / norm
    return signed, abs(signed)


def evaluate_resonance_proximity(
    working_points,
    resonance_lines=None,
    xlims=(4.0, 4.5),
    ylims=(3.5, 4.0),
    orders=(1, 2, 3, 4),
    periodicity=10,
    qx_column="predicted_qx",
    qy_column="predicted_qy",
):
    """
    Find the nearest generated resonance for each working point.
    """

    if not isinstance(working_points, pd.DataFrame):
        working_points = pd.DataFrame(working_points)

    if resonance_lines is None:
        resonance_lines = generate_resonance_lines(
            xlims=xlims,
            ylims=ylims,
            orders=orders,
            periodicity=periodicity,
        )

    if qx_column not in working_points.columns or qy_column not in working_points.columns:
        raise ValueError(f"working_points must contain {qx_column!r} and {qy_column!r}.")

    rows = []
    for point_index, point in working_points.iterrows():
        qx = point[qx_column]
        qy = point[qy_column]
        if pd.isna(qx) or pd.isna(qy):
            continue

        candidates = []
        for _, line in resonance_lines.iterrows():
            signed, absolute = _point_line_distance(
                qx=qx,
                qy=qy,
                nx=line["nx"],
                ny=line["ny"],
                resonance_sum=line["resonance_sum"],
            )
            candidates.append((absolute, signed, line))

        if not candidates:
            continue

        _, signed, nearest = min(candidates, key=lambda item: item[0])
        rows.append(
            {
                "point_index": int(point_index),
                "cycle_time_ms": point.get("cycle_time_ms", None),
                "qx": float(qx),
                "qy": float(qy),
                "nearest_order": int(nearest["order"]),
                "nearest_nx": int(nearest["nx"]),
                "nearest_ny": int(nearest["ny"]),
                "nearest_resonance_sum": int(nearest["resonance_sum"]),
                "nearest_label": nearest["label"],
                "nearest_systematic": bool(nearest["systematic"]),
                "signed_distance": float(signed),
                "absolute_distance": abs(float(signed)),
            }
        )

    return pd.DataFrame(rows)


def make_tune_diagram_inputs(
    tune_programme_df,
    xlims=(4.0, 4.5),
    ylims=(3.5, 4.0),
    orders=(1, 2, 3, 4),
    periodicity=10,
    proximity_qx_column="predicted_qx",
    proximity_qy_column="predicted_qy",
):
    """
    Build all DataFrame inputs needed for a tune diagram.
    """

    programme = build_tune_programme_table(tune_programme_df)
    working_points = build_working_point_table(programme)
    resonance_lines = generate_resonance_lines(
        xlims=xlims,
        ylims=ylims,
        orders=orders,
        periodicity=periodicity,
    )
    proximity = evaluate_resonance_proximity(
        working_points,
        resonance_lines=resonance_lines,
        qx_column=proximity_qx_column,
        qy_column=proximity_qy_column,
    )

    return {
        "programme": programme,
        "working_points": working_points,
        "resonance_lines": resonance_lines,
        "resonance_proximity": proximity,
        "metadata": {
            "xlims": tuple(float(value) for value in xlims),
            "ylims": tuple(float(value) for value in ylims),
            "orders": tuple(int(order) for order in orders),
            "periodicity": int(periodicity),
        },
    }
