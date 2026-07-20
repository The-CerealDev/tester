"""
aperture.py

Read-only aperture and beam-clearance utilities for the ISIS RCS optics GUI
backend.

This layer consumes real MAD-X aperture/TWISS/envelope outputs and the ISIS
source aperture spreadsheet. It does not generate model data itself.
"""

from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd

from .envelope import EnvelopeResult
from .madx_model import MadxModel


DEFAULT_SOURCE_APERTURE_CSV = (
    Path(__file__).resolve().parent / "data" / "aperture" / "jvt_synch_aperture.csv"
)
APERTURE_EPS = 1.0e-12


@dataclass
class ApertureResult:
    """
    Result container for one aperture-margin evaluation.
    """

    aperture_df: pd.DataFrame
    aligned_df: pd.DataFrame
    summary_df: pd.DataFrame
    source: str
    metadata: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def to_summary_dataframe(self):
        return self.summary_df.copy()

    def to_plot_dataframe(self):
        return self.aligned_df.copy()


def _normalise_name(value):
    return str(value).strip().lower().split(":")[0]


def _coerce_numeric(df, columns):
    out = df.copy()
    for col in columns:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def read_source_aperture_csv(path=None):
    """
    Read the ISIS source aperture spreadsheet exported as CSV.
    """

    if path is None:
        resource = resources.files(__package__).joinpath(
            "data",
            "aperture",
            "jvt_synch_aperture.csv",
        )
        with resources.as_file(resource) as resource_path:
            df = pd.read_csv(resource_path)
            return normalise_source_aperture_table(df, source_path=resource_path)
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing source aperture CSV: {path}")

    df = pd.read_csv(path)
    return normalise_source_aperture_table(df, source_path=path)


def normalise_source_aperture_table(df, source_path=None):
    """
    Normalise JVT/ISIS source aperture spreadsheet columns.

    The source spreadsheet stores semi-apertures in mm. This returns metres.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    out = df.copy()
    out.columns = [str(col).strip().lower() for col in out.columns]

    required = ("dist_datum_d", "semi_ap_h", "semi_ap_v")
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f"Source aperture table is missing columns: {missing}")

    if "element" in out.columns and "name" not in out.columns:
        out["name"] = out["element"].astype(str)

    out = _coerce_numeric(out, ["dist_datum_d", "semi_ap_h", "semi_ap_v"])
    out["s"] = out["dist_datum_d"].astype(float)
    out["aperture_x_m"] = 1.0e-3 * out["semi_ap_h"].astype(float)
    out["aperture_y_m"] = 1.0e-3 * out["semi_ap_v"].astype(float)
    out["aperture_x_mm"] = 1.0e3 * out["aperture_x_m"]
    out["aperture_y_mm"] = 1.0e3 * out["aperture_y_m"]
    out["aperture_source"] = "source_spreadsheet"
    out["apertype"] = "RECTANGLE"
    out["name_key"] = out.get("name", pd.Series(index=out.index, dtype=str)).map(_normalise_name)

    if source_path is not None:
        out["source_path"] = str(source_path)

    return _validate_normalised_aperture(out, source="source spreadsheet")


def normalise_madx_aperture_table(df, drop_zero_apertures=True):
    """
    Normalise a MAD-X APERTURE table.

    MAD-X reports zero aperture for drifts that cannot hold aperture metadata.
    Those rows are invalid for margin calculations and are dropped by default.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    out = df.copy()
    out.columns = [str(col).strip().lower() for col in out.columns]

    if "s" not in out.columns:
        raise ValueError("MAD-X aperture table is missing column: s")
    if "name" not in out.columns:
        raise ValueError("MAD-X aperture table is missing column: name")

    out["apertype"] = out.get("apertype", "UNKNOWN")
    out["apertype"] = out["apertype"].astype(str).str.strip().str.strip('"').str.upper()

    if {"aper_1", "aper_2"}.issubset(out.columns):
        out = _coerce_numeric(out, ["s", "aper_1", "aper_2"])
        out["aperture_x_m"] = out["aper_1"].astype(float)
        out["aperture_y_m"] = out["aper_2"].astype(float)
        circle_mask = out["apertype"] == "CIRCLE"
        if circle_mask.any():
            radius = out.loc[circle_mask, "aper_1"].astype(float)
            out.loc[circle_mask, "aperture_x_m"] = radius
            out.loc[circle_mask, "aperture_y_m"] = radius
    elif {"n1x_m", "n1y_m"}.issubset(out.columns):
        out = _coerce_numeric(out, ["s", "n1x_m", "n1y_m"])
        out["aperture_x_m"] = out["n1x_m"].astype(float)
        out["aperture_y_m"] = out["n1y_m"].astype(float)
    else:
        raise ValueError(
            "MAD-X aperture table must contain either aper_1/aper_2 or n1x_m/n1y_m."
        )

    if drop_zero_apertures:
        valid = (out["aperture_x_m"] > APERTURE_EPS) & (out["aperture_y_m"] > APERTURE_EPS)
        out = out.loc[valid].copy()

    unsupported = sorted(
        {
            str(value).upper()
            for value in out["apertype"].dropna().unique()
            if str(value).upper() not in ("RECTANGLE", "CIRCLE", "UNKNOWN")
        }
    )
    if unsupported:
        raise ValueError(f"Unsupported aperture types: {unsupported}")

    out["aperture_x_mm"] = 1.0e3 * out["aperture_x_m"]
    out["aperture_y_mm"] = 1.0e3 * out["aperture_y_m"]
    out["aperture_source"] = "madx_aperture"
    out["name_key"] = out["name"].map(_normalise_name)

    return _validate_normalised_aperture(out, source="MAD-X aperture")


def _validate_normalised_aperture(df, source):
    required = ("s", "aperture_x_m", "aperture_y_m", "aperture_source")
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{source} table is missing columns: {missing}")

    out = _coerce_numeric(df, ["s", "aperture_x_m", "aperture_y_m"])
    numeric = out[["s", "aperture_x_m", "aperture_y_m"]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError(f"{source} table contains non-finite aperture values.")

    if ((out["aperture_x_m"] <= 0.0) | (out["aperture_y_m"] <= 0.0)).any():
        raise ValueError(f"{source} table contains zero or negative aperture values.")

    return out.sort_values("s").reset_index(drop=True)


def extract_madx_aperture_from_machine_state(
    machine_state,
    lattice_folder="../Lattice_Files/00_Simplified_Lattice",
    aperture_file="ISIS.aperture",
    output_dir="./aperture_tests/madx_aperture",
    sequence_name="synchrotron",
    interval=0.1,
):
    """
    Run MAD-X for a MachineState and return the normalised APERTURE table.
    """

    model = MadxModel(
        lattice_folder=lattice_folder,
        sequence_name=sequence_name,
        aperture_file=aperture_file,
        output_dir=output_dir,
    )
    model.load_lattice(use_sequence=False)
    model.apply_machine_state(machine_state)
    model.use_sequence()
    aperture_df = model.run_aperture(interval=interval)
    normalised = normalise_madx_aperture_table(aperture_df)

    return ApertureResult(
        aperture_df=normalised,
        aligned_df=pd.DataFrame(),
        summary_df=pd.DataFrame(),
        source="madx_aperture",
        metadata=model.get_metadata(),
        warnings=[],
    )


def _envelope_dataframe(envelope):
    if isinstance(envelope, EnvelopeResult):
        return envelope.envelope_df.copy()
    if isinstance(envelope, pd.DataFrame):
        return envelope.copy()
    raise TypeError("envelope must be an EnvelopeResult or pandas DataFrame.")


def align_aperture_to_envelope(aperture_df, envelope):
    """
    Interpolate aperture half-widths onto the envelope/TWISS grid.
    """

    aperture = _validate_normalised_aperture(aperture_df, source="aperture")
    env = _envelope_dataframe(envelope)
    env.columns = [str(col).lower() for col in env.columns]

    required_env = (
        "name",
        "s",
        "envelope_x_plus_m",
        "envelope_x_minus_m",
        "envelope_y_plus_m",
        "envelope_y_minus_m",
        "x",
        "y",
    )
    missing = [col for col in required_env if col not in env.columns]
    if missing:
        raise ValueError(f"Envelope table is missing columns: {missing}")

    env = _coerce_numeric(
        env,
        [
            "s",
            "x",
            "y",
            "envelope_x_plus_m",
            "envelope_x_minus_m",
            "envelope_y_plus_m",
            "envelope_y_minus_m",
        ],
    )
    env["name_key"] = env["name"].map(_normalise_name)

    ap = aperture.drop_duplicates("name_key", keep="first")
    merged = env.merge(
        ap[["name_key", "aperture_x_m", "aperture_y_m", "aperture_source"]],
        on="name_key",
        how="left",
    )
    merged["aperture_alignment"] = np.where(
        merged["aperture_x_m"].notna() & merged["aperture_y_m"].notna(),
        "name",
        "s_interpolated",
    )

    source_s = aperture["s"].to_numpy(dtype=float)
    source_x = aperture["aperture_x_m"].to_numpy(dtype=float)
    source_y = aperture["aperture_y_m"].to_numpy(dtype=float)
    order = np.argsort(source_s)
    source_s = source_s[order]
    source_x = source_x[order]
    source_y = source_y[order]

    missing_mask = merged["aperture_x_m"].isna() | merged["aperture_y_m"].isna()
    if missing_mask.any():
        target_s = merged.loc[missing_mask, "s"].to_numpy(dtype=float)
        merged.loc[missing_mask, "aperture_x_m"] = np.interp(target_s, source_s, source_x)
        merged.loc[missing_mask, "aperture_y_m"] = np.interp(target_s, source_s, source_y)
        merged.loc[missing_mask, "aperture_source"] = aperture["aperture_source"].iloc[0]

    merged["aperture_x_plus_m"] = merged["aperture_x_m"]
    merged["aperture_x_minus_m"] = -merged["aperture_x_m"]
    merged["aperture_y_plus_m"] = merged["aperture_y_m"]
    merged["aperture_y_minus_m"] = -merged["aperture_y_m"]

    aperture_mm_columns = {
        "aperture_x_m": "aperture_x_mm",
        "aperture_y_m": "aperture_y_mm",
        "aperture_x_plus_m": "aperture_x_plus_mm",
        "aperture_x_minus_m": "aperture_x_minus_mm",
        "aperture_y_plus_m": "aperture_y_plus_mm",
        "aperture_y_minus_m": "aperture_y_minus_mm",
    }
    for metres_col, millimetres_col in aperture_mm_columns.items():
        merged[millimetres_col] = 1.0e3 * merged[metres_col]

    return merged


def compute_aperture_margins(aligned_df):
    """
    Compute beam-envelope clearance against rectangular aperture.
    """

    df = aligned_df.copy()
    required = (
        "aperture_x_plus_m",
        "aperture_x_minus_m",
        "aperture_y_plus_m",
        "aperture_y_minus_m",
        "envelope_x_plus_m",
        "envelope_x_minus_m",
        "envelope_y_plus_m",
        "envelope_y_minus_m",
    )
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Aligned table is missing columns: {missing}")

    df["margin_x_plus_m"] = df["aperture_x_plus_m"] - df["envelope_x_plus_m"]
    df["margin_x_minus_m"] = df["envelope_x_minus_m"] - df["aperture_x_minus_m"]
    df["margin_y_plus_m"] = df["aperture_y_plus_m"] - df["envelope_y_plus_m"]
    df["margin_y_minus_m"] = df["envelope_y_minus_m"] - df["aperture_y_minus_m"]
    df["margin_x_m"] = df[["margin_x_plus_m", "margin_x_minus_m"]].min(axis=1)
    df["margin_y_m"] = df[["margin_y_plus_m", "margin_y_minus_m"]].min(axis=1)
    df["margin_min_m"] = df[["margin_x_m", "margin_y_m"]].min(axis=1)

    side_cols = [
        "margin_x_plus_m",
        "margin_x_minus_m",
        "margin_y_plus_m",
        "margin_y_minus_m",
    ]
    df["limiting_side"] = df[side_cols].idxmin(axis=1).str.replace("margin_", "").str.replace("_m", "")
    df["limiting_plane"] = df["limiting_side"].str[0]

    margin_mm_columns = {
        "margin_x_plus_m": "margin_x_plus_mm",
        "margin_x_minus_m": "margin_x_minus_mm",
        "margin_y_plus_m": "margin_y_plus_mm",
        "margin_y_minus_m": "margin_y_minus_mm",
        "margin_x_m": "margin_x_mm",
        "margin_y_m": "margin_y_mm",
        "margin_min_m": "margin_min_mm",
    }
    for metres_col, millimetres_col in margin_mm_columns.items():
        df[millimetres_col] = 1.0e3 * df[metres_col]

    return df


def summarise_aperture_margins(margin_df, label=None):
    """
    Return the limiting aperture locations for GUI display.
    """

    if "margin_min_m" not in margin_df.columns:
        raise ValueError("margin_df must contain margin_min_m.")

    idx = margin_df["margin_min_m"].idxmin()
    row = margin_df.loc[idx]
    return pd.DataFrame(
        [
            {
                "label": label,
                "name": row["name"],
                "s_m": float(row["s"]),
                "limiting_plane": row["limiting_plane"],
                "limiting_side": row["limiting_side"],
                "margin_min_m": float(row["margin_min_m"]),
                "margin_min_mm": float(row["margin_min_mm"]),
                "aperture_source": row.get("aperture_source", None),
                "aperture_alignment": row.get("aperture_alignment", None),
            }
        ]
    )


def evaluate_aperture_margins(aperture_df, envelope, label=None, source="aperture"):
    """
    Align aperture to an envelope result and calculate clearance margins.
    """

    aligned = align_aperture_to_envelope(aperture_df, envelope)
    margins = compute_aperture_margins(aligned)
    summary = summarise_aperture_margins(margins, label=label)
    warnings = []
    if (margins["margin_min_m"] < 0.0).any():
        warnings.append("Envelope exceeds aperture at one or more locations.")

    return ApertureResult(
        aperture_df=aperture_df.copy(),
        aligned_df=margins,
        summary_df=summary,
        source=source,
        metadata={"label": label, "source": source, "n_rows": len(margins)},
        warnings=warnings,
    )


def compare_aperture_summaries(results):
    rows = [result.summary_df for result in results]
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def plot_aperture_envelope(result, plane="x", ax=None, title=None, show_orbit=True):
    """
    Plot aperture and envelope bounds for one plane.
    """

    import matplotlib.pyplot as plt

    plane = str(plane).lower()
    if plane not in {"x", "y"}:
        raise ValueError("plane must be 'x' or 'y'.")

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    df = result.aligned_df
    s = df["s"]

    aperture_color = "black"
    envelope_color = "tab:blue"
    source_color = "tab:orange"

    ax.plot(
        s,
        df[f"aperture_{plane}_plus_mm"],
        color=aperture_color,
        linewidth=1.3,
        linestyle="-",
        label="MAD-X/source aperture +",
    )
    ax.plot(
        s,
        df[f"aperture_{plane}_minus_mm"],
        color=aperture_color,
        linewidth=1.3,
        linestyle="-",
        label="MAD-X/source aperture -",
    )
    ax.plot(
        s,
        df[f"envelope_{plane}_plus_mm"],
        color=envelope_color,
        linewidth=1.0,
        linestyle="-",
        label="envelope +",
    )
    ax.plot(
        s,
        df[f"envelope_{plane}_minus_mm"],
        color=envelope_color,
        linewidth=1.0,
        linestyle="-",
        label="envelope -",
    )
    if show_orbit:
        ax.plot(s, df[f"orbit_{plane}_mm"], color=source_color, linewidth=1.0, label="closed orbit")

    ax.set_xlabel("s [m]")
    ax.set_ylabel(f"{plane} [mm]")
    ax.set_title(title or f"{plane.upper()} aperture and envelope")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return ax


def plot_aperture_source_overlay(madx_result, source_result, plane="x", ax=None, title=None):
    """
    Plot MAD-X extracted aperture and source-spreadsheet aperture together.
    """

    import matplotlib.pyplot as plt

    plane = str(plane).lower()
    if plane not in {"x", "y"}:
        raise ValueError("plane must be 'x' or 'y'.")

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    for result, color, label in (
        (madx_result, "black", "MAD-X APERTURE"),
        (source_result, "tab:orange", "source spreadsheet"),
    ):
        df = result.aligned_df
        ax.plot(df["s"], df[f"aperture_{plane}_plus_mm"], color=color, linewidth=1.1, label=f"{label} +")
        ax.plot(df["s"], df[f"aperture_{plane}_minus_mm"], color=color, linewidth=1.1, label=f"{label} -")

    ax.set_xlabel("s [m]")
    ax.set_ylabel(f"{plane} aperture [mm]")
    ax.set_title(title or f"{plane.upper()} MAD-X vs source aperture")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return ax


def plot_margin(result, ax=None, title=None):
    """
    Plot horizontal, vertical and limiting aperture margin.
    """

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    df = result.aligned_df
    ax.plot(df["s"], df["margin_x_mm"], label="x margin")
    ax.plot(df["s"], df["margin_y_mm"], label="y margin")
    ax.plot(df["s"], df["margin_min_mm"], color="black", linewidth=1.3, label="limiting margin")
    ax.axhline(0.0, color="red", linewidth=1.0, linestyle="--", label="zero clearance")
    ax.set_xlabel("s [m]")
    ax.set_ylabel("margin [mm]")
    ax.set_title(title or "Aperture margin")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return ax


def plot_aperture_envelope_with_margin(
    result,
    plane="x",
    ax=None,
    title=None,
    show_orbit=True,
    margin_planes=None,
):
    """
    Plot aperture, filled beam envelope, closed orbit and aperture margins.

    Visual convention:
      - aperture: black boundary lines
      - envelope: translucent filled beam envelope
      - closed orbit: orange line
      - margin: purple dashed lines
      - aperture-hit parts of the margin: red dashed segments

    By default, only the margin for the plotted plane is shown:
      - plane="x" shows only horizontal/x margins
      - plane="y" shows only vertical/y margins

    The + margin line is plotted as +margin_plus.
    The - margin line is plotted as -margin_minus for visual symmetry.

    Red logic:
      - + margin line is red where the plotted + margin is below zero.
      - - margin line is red where the plotted - margin is above zero.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    plane = str(plane).lower()
    if plane not in {"x", "y"}:
        raise ValueError("plane must be 'x' or 'y'.")

    if margin_planes is None:
        margin_planes = (plane,)

    margin_planes = tuple(str(item).lower() for item in margin_planes)
    invalid_margin_planes = sorted(set(margin_planes) - {"x", "y"})
    if invalid_margin_planes:
        raise ValueError("margin_planes entries must be 'x' or 'y'.")

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    df = result.aligned_df

    required = [
        "s",
        f"aperture_{plane}_plus_mm",
        f"aperture_{plane}_minus_mm",
        f"envelope_{plane}_plus_mm",
        f"envelope_{plane}_minus_mm",
    ]
    if show_orbit:
        required.append(f"orbit_{plane}_mm")

    for margin_plane in margin_planes:
        required.extend(
            [
                f"margin_{margin_plane}_plus_mm",
                f"margin_{margin_plane}_minus_mm",
            ]
        )

    missing = [col for col in required if col not in df.columns]
    if missing:
        raise KeyError(f"Required columns not found: {missing}")

    s = df["s"].to_numpy(dtype=float)

    aperture_plus = df[f"aperture_{plane}_plus_mm"].to_numpy(dtype=float)
    aperture_minus = df[f"aperture_{plane}_minus_mm"].to_numpy(dtype=float)
    envelope_plus = df[f"envelope_{plane}_plus_mm"].to_numpy(dtype=float)
    envelope_minus = df[f"envelope_{plane}_minus_mm"].to_numpy(dtype=float)

    aperture_color = "black"
    envelope_color = "tab:blue"
    orbit_color = "tab:orange"
    margin_color = "tab:purple"
    hit_color = "tab:red"

    def _red_segments_where(ax, s_values, y_values, hit_when, *, linewidth, linestyle):
        """
        Overlay red line segments where the plotted margin violates aperture.

        hit_when="below" colours only y < 0.
        hit_when="above" colours only y > 0.

        Segments are clipped to y=0 by linear interpolation.
        """

        s_values = np.asarray(s_values, dtype=float)
        y_values = np.asarray(y_values, dtype=float)

        if hit_when not in {"below", "above"}:
            raise ValueError("hit_when must be 'below' or 'above'.")

        def is_hit(y):
            if hit_when == "below":
                return y < 0.0
            return y > 0.0

        segments = []

        for i in range(len(s_values) - 1):
            x0 = s_values[i]
            x1 = s_values[i + 1]
            y0 = y_values[i]
            y1 = y_values[i + 1]

            if not (
                np.isfinite(x0)
                and np.isfinite(x1)
                and np.isfinite(y0)
                and np.isfinite(y1)
            ):
                continue

            hit0 = is_hit(y0)
            hit1 = is_hit(y1)

            if hit0 and hit1:
                segments.append([[x0, y0], [x1, y1]])
                continue

            if not hit0 and not hit1:
                continue

            if y1 == y0:
                continue

            t = -y0 / (y1 - y0)
            x_cross = x0 + t * (x1 - x0)

            if hit0 and not hit1:
                segments.append([[x0, y0], [x_cross, 0.0]])
            elif not hit0 and hit1:
                segments.append([[x_cross, 0.0], [x1, y1]])

        if not segments:
            return False

        collection = LineCollection(
            segments,
            colors=hit_color,
            linewidths=linewidth,
            linestyles=linestyle,
            zorder=10,
        )
        ax.add_collection(collection)

        return True

    ax.plot(
        s,
        aperture_plus,
        color=aperture_color,
        linewidth=1.3,
        linestyle="-",
        label="aperture",
    )
    ax.plot(
        s,
        aperture_minus,
        color=aperture_color,
        linewidth=1.3,
        linestyle="-",
        label="_nolegend_",
    )

    ax.fill_between(
        s,
        envelope_minus,
        envelope_plus,
        color=envelope_color,
        alpha=0.20,
        linewidth=0.0,
        label="envelope",
    )

    if show_orbit:
        ax.plot(
            s,
            df[f"orbit_{plane}_mm"],
            color=orbit_color,
            linewidth=1.2,
            linestyle="-",
            label="closed orbit",
        )

    margin_linewidth = 1.8
    margin_hit_linewidth = 2.4
    margin_linestyle = "--"
    margin_label_used = False

    for margin_plane in margin_planes:
        plus_margin = df[f"margin_{margin_plane}_plus_mm"].to_numpy(dtype=float)
        minus_margin = df[f"margin_{margin_plane}_minus_mm"].to_numpy(dtype=float)

        plus_y = plus_margin
        minus_y = -minus_margin

        ax.plot(
            s,
            plus_y,
            color=margin_color,
            linewidth=margin_linewidth,
            linestyle=margin_linestyle,
            label="margin" if not margin_label_used else "_nolegend_",
        )
        margin_label_used = True

        ax.plot(
            s,
            minus_y,
            color=margin_color,
            linewidth=margin_linewidth,
            linestyle=margin_linestyle,
            label="_nolegend_",
        )

        _red_segments_where(
            ax,
            s,
            plus_y,
            hit_when="below",
            linewidth=margin_hit_linewidth,
            linestyle=margin_linestyle,
        )

        _red_segments_where(
            ax,
            s,
            minus_y,
            hit_when="above",
            linewidth=margin_hit_linewidth,
            linestyle=margin_linestyle,
        )

    ax.axhline(
        0.0,
        color="0.4",
        linewidth=0.9,
        linestyle=":",
        label="_nolegend_",
    )

    ax.set_xlabel("s [m]")
    ax.set_ylabel(f"{plane} [mm]")
    ax.set_title(title or f"{plane.upper()} aperture, envelope and margin")
    ax.grid(True, alpha=0.3)

    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0.0,
    )

    ax.figure.subplots_adjust(right=0.78)

    return ax
