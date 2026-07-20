"""
tune_plots.py

Plotting helpers for tune matching and harmonic tune development.

The functions return Matplotlib figure/axes objects and do not save files by
default, which keeps them usable from notebooks and future Streamlit callbacks.
"""

from collections import OrderedDict
import re

import matplotlib.cm as cm
import matplotlib.lines as mlines
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from machine_state_defaults import DEFAULT_HARMONICS, DEFAULT_TQ_GCAL


DEFAULT_HQCAL = 1.25
DEFAULT_INJECTION_BRHO_TM = 1.23


def _as_float_array(values):
    return np.asarray(values, dtype=float)


def _normalise_harmonic_values(harmonics=None):
    values = OrderedDict(DEFAULT_HARMONICS)
    lookup = {
        re.sub(r"[^A-Z0-9]", "", key.upper()): key
        for key in values
    }
    if harmonics is not None:
        for key, value in dict(harmonics).items():
            lookup_key = re.sub(r"[^A-Z0-9]", "", str(key).upper())
            if lookup_key not in lookup:
                raise KeyError(f"Unknown harmonic tune variable: {key}")
            values[lookup[lookup_key]] = float(value)
    return values


def _harmonic_scale(brho_Tm=DEFAULT_INJECTION_BRHO_TM, hqcal=DEFAULT_HQCAL, tqgcal=DEFAULT_TQ_GCAL):
    return float(hqcal) * float(tqgcal) / float(brho_Tm)


def _harmonic_wave(
    harmonics=None,
    superperiod_coordinate=None,
    brho_Tm=DEFAULT_INJECTION_BRHO_TM,
    hqcal=DEFAULT_HQCAL,
    tqgcal=DEFAULT_TQ_GCAL,
):
    values = _normalise_harmonic_values(harmonics)
    if superperiod_coordinate is None:
        superperiod_coordinate = np.linspace(0.0, 10.0, 1001)

    superperiod_coordinate = np.asarray(superperiod_coordinate, dtype=float)
    phase = (superperiod_coordinate / 10.0) * 2.0 * np.pi
    scale = _harmonic_scale(brho_Tm=brho_Tm, hqcal=hqcal, tqgcal=tqgcal)

    qtd_k = (
        values["D7COS"] * scale * np.cos(7.0 * phase)
        + values["D7SIN"] * scale * np.sin(7.0 * phase)
        + values["D8COS"] * scale * np.cos(8.0 * phase)
        + values["D8SIN"] * scale * np.sin(8.0 * phase)
    )
    qtf_k = (
        values["F8COS"] * scale * np.cos(8.0 * phase)
        + values["F8SIN"] * scale * np.sin(8.0 * phase)
        + values["F9COS"] * scale * np.cos(9.0 * phase)
        + values["F9SIN"] * scale * np.sin(9.0 * phase)
    )

    return superperiod_coordinate, qtd_k, qtf_k


def calculate_harmonic_trim_quad_pattern(
    harmonics=None,
    brho_Tm=DEFAULT_INJECTION_BRHO_TM,
    hqcal=DEFAULT_HQCAL,
    tqgcal=DEFAULT_TQ_GCAL,
):
    """
    Calculate expected QTD/QTF harmonic K offsets by superperiod.
    """

    superperiod = np.arange(10, dtype=int)
    _, qtd_k, qtf_k = _harmonic_wave(
        harmonics=harmonics,
        superperiod_coordinate=superperiod,
        brho_Tm=brho_Tm,
        hqcal=hqcal,
        tqgcal=tqgcal,
    )

    return pd.DataFrame(
        {
            "superperiod": superperiod,
            "qtd_expected_k": qtd_k,
            "qtf_expected_k": qtf_k,
            "qtd_label": [f"R{i}_QTD" for i in superperiod],
            "qtf_label": [f"R{i}_QTF" for i in superperiod],
        }
    )


def extract_programmed_trim_quad_pattern(
    twiss_df,
    base_kqtd=None,
    base_kqtf=None,
):
    """
    Extract programmed QTD/QTF K values from a MAD-X TWISS table.

    The returned k1 values are local quadrupole strengths. When base_kqtd or
    base_kqtf are supplied, delta_k columns expose the harmonic offsets.
    """

    if not isinstance(twiss_df, pd.DataFrame):
        raise TypeError("twiss_df must be a pandas DataFrame.")

    required = ("name", "l", "k1l")
    missing = [column for column in required if column not in twiss_df.columns]
    if missing:
        raise ValueError(f"TWISS DataFrame is missing columns: {missing}")

    rows = []
    for _, row in twiss_df.iterrows():
        name = str(row["name"]).lower()
        match = re.search(r"(?:r|sp)([0-9])_?q(t[fd])", name)
        if match is None:
            continue

        length = float(row["l"])
        k1l = float(row["k1l"])
        k1 = k1l / length if length != 0.0 else k1l
        family = match.group(2)
        base = base_kqtd if family == "td" else base_kqtf
        rows.append(
            {
                "superperiod": int(match.group(1)),
                "family": "qtd" if family == "td" else "qtf",
                "name": str(row["name"]),
                "l": length,
                "k1l": k1l,
                "k1": k1,
                "delta_k": None if base is None else k1 - float(base),
            }
        )

    pattern = pd.DataFrame(rows)
    if pattern.empty:
        raise ValueError("No RnQTD/RnQTF trim quadrupoles found in TWISS table.")

    return pattern.sort_values(["family", "superperiod"]).reset_index(drop=True)


def _trim_quad_positions_by_family(twiss_df):
    if not isinstance(twiss_df, pd.DataFrame) or "name" not in twiss_df.columns or "s" not in twiss_df.columns:
        return {}

    rows = []
    for _, row in twiss_df.iterrows():
        name = str(row["name"]).lower()
        match = re.search(r"(?:r|sp)([0-9])_?q(t[fd])", name)
        if match is None:
            continue

        rows.append(
            {
                "superperiod": int(match.group(1)),
                "family": "qtd" if match.group(2) == "td" else "qtf",
                "s": float(row["s"]),
            }
        )

    if not rows:
        return {}

    positions = pd.DataFrame(rows).sort_values(["family", "superperiod"])
    return {
        family: family_df[["superperiod", "s"]].reset_index(drop=True)
        for family, family_df in positions.groupby("family")
    }


def _plot_family_wave_on_s_axis(
    ax,
    family_positions,
    harmonics,
    family,
    color,
    circumference,
    brho_Tm=DEFAULT_INJECTION_BRHO_TM,
    hqcal=DEFAULT_HQCAL,
    tqgcal=DEFAULT_TQ_GCAL,
    smooth_points=1001,
):
    superperiod = family_positions["superperiod"].astype(float).to_numpy()
    s_positions = family_positions["s"].astype(float).to_numpy()
    if len(superperiod) < 2:
        return

    dense_superperiod = np.linspace(0.0, 10.0, int(smooth_points))
    _, qtd_wave, qtf_wave = _harmonic_wave(
        harmonics=harmonics,
        superperiod_coordinate=dense_superperiod,
        brho_Tm=brho_Tm,
        hqcal=hqcal,
        tqgcal=tqgcal,
    )
    dense_values = qtd_wave if family == "qtd" else qtf_wave

    phase_points = np.append(superperiod, 10.0)
    s_points = np.append(s_positions, s_positions[0] + circumference)
    dense_s = np.interp(dense_superperiod, phase_points, s_points)
    dense_s_wrapped = np.mod(dense_s, circumference)
    wrap_indices = np.where(np.diff(dense_s_wrapped) < 0.0)[0] + 1
    starts = np.r_[0, wrap_indices]
    stops = np.r_[wrap_indices, len(dense_s_wrapped)]

    for start, stop in zip(starts, stops):
        label = f"Expected {family.upper()} K" if start == 0 else None
        ax.plot(
            dense_s_wrapped[start:stop],
            dense_values[start:stop],
            color=color,
            lw=1.2,
            label=label,
        )

    _, qtd_at_elements, qtf_at_elements = _harmonic_wave(
        harmonics=harmonics,
        superperiod_coordinate=superperiod,
        brho_Tm=brho_Tm,
        hqcal=hqcal,
        tqgcal=tqgcal,
    )
    element_values = qtd_at_elements if family == "qtd" else qtf_at_elements
    ax.scatter(s_positions, element_values, color=color, s=24, zorder=3)


class ResonanceLines:
    """
    Draw tune resonance lines on an existing Matplotlib axis.
    """

    def __init__(self, qx_range, qy_range, orders=(1, 2, 3, 4), periodicity=10):
        self.qx_min = float(np.min(qx_range))
        self.qx_max = float(np.max(qx_range))
        self.qy_min = float(np.min(qy_range))
        self.qy_max = float(np.max(qy_range))
        self.orders = tuple(int(order) for order in orders)
        self.periodicity = int(periodicity)

    def plot(self, ax):
        ax.set_xlim(self.qx_min, self.qx_max)
        ax.set_ylim(self.qy_min, self.qy_max)
        ax.set_xlabel(r"Horizontal Tune $Q_x$")
        ax.set_ylabel(r"Vertical Tune $Q_y$")

        for order in self.orders:
            terms = np.arange(-order, order + 1)
            nx_values = order - np.abs(terms)
            ny_values = terms

            for nx, ny in zip(nx_values, ny_values):
                if nx == 0 and ny == 0:
                    continue

                corners = np.array(
                    [
                        nx * self.qx_min + ny * self.qy_min,
                        nx * self.qx_max + ny * self.qy_min,
                        nx * self.qx_min + ny * self.qy_max,
                        nx * self.qx_max + ny * self.qy_max,
                    ]
                )

                for resonance_sum in range(int(np.floor(corners.min())), int(np.ceil(corners.max())) + 1):
                    systematic = resonance_sum % self.periodicity == 0
                    color = "red" if systematic else "blue"
                    width = 1.0 if systematic else 0.5
                    style = "--" if ny % 2 else "-"

                    if ny == 0:
                        if nx == 0:
                            continue
                        qx = float(resonance_sum) / float(nx)
                        ax.plot([qx, qx], [self.qy_min, self.qy_max], color=color, lw=width, ls=style, alpha=0.8)
                    else:
                        y0 = (resonance_sum - nx * self.qx_min) / ny
                        y1 = (resonance_sum - nx * self.qx_max) / ny
                        ax.plot([self.qx_min, self.qx_max], [y0, y1], color=color, lw=width, ls=style, alpha=0.8)

        return ax


def plot_resonance_working_points(
    requested_qx,
    requested_qy,
    actual_qx=None,
    actual_qy=None,
    matched_qx=None,
    matched_qy=None,
    cycle_times=None,
    xlims=(4.0, 4.5),
    ylims=(3.5, 4.0),
    orders=(1, 2, 3, 4),
    periodicity=10,
    ax=None,
):
    """
    Plot requested, actual and matched tune working points with resonances.
    """

    requested_qx = _as_float_array(requested_qx)
    requested_qy = _as_float_array(requested_qy)
    if cycle_times is None:
        cycle_times = np.arange(len(requested_qx), dtype=float)
    cycle_times = _as_float_array(cycle_times)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
    else:
        fig = ax.figure

    ResonanceLines(xlims, ylims, orders=orders, periodicity=periodicity).plot(ax)

    unique_times = np.unique(cycle_times)
    cmap = ListedColormap(cm.rainbow(np.linspace(0.0, 1.0, len(unique_times))))
    bounds = np.append(unique_times - 0.5, unique_times[-1] + 0.5)
    norm = BoundaryNorm(bounds, cmap.N)

    scatter = ax.scatter(requested_qx, requested_qy, c=cycle_times, cmap=cmap, norm=norm, marker="o", s=45, zorder=3)
    ax.plot(requested_qx, requested_qy, ls=":", lw=0.8, color="black")
    legend_handles = [
        mlines.Line2D(
            [],
            [],
            color="black",
            marker="o",
            linestyle=":",
            markersize=7,
            lw=0.8,
            label="Set",
        )
    ]

    if actual_qx is not None and actual_qy is not None:
        ax.scatter(_as_float_array(actual_qx), _as_float_array(actual_qy), c=cycle_times, cmap=cmap, norm=norm, marker="+", s=60, zorder=4)
        ax.plot(_as_float_array(actual_qx), _as_float_array(actual_qy), ls="--", lw=0.8, color="black")
        legend_handles.append(
            mlines.Line2D(
                [],
                [],
                color="black",
                marker="+",
                linestyle="--",
                markersize=9,
                lw=0.8,
                label="Measured",
            )
        )

    if matched_qx is not None and matched_qy is not None:
        ax.scatter(_as_float_array(matched_qx), _as_float_array(matched_qy), c=cycle_times, cmap=cmap, norm=norm, marker="x", s=55, zorder=5)
        ax.plot(_as_float_array(matched_qx), _as_float_array(matched_qy), ls="-.", lw=0.8, color="black")
        legend_handles.append(
            mlines.Line2D(
                [],
                [],
                color="black",
                marker="x",
                linestyle="-.",
                markersize=8,
                lw=0.8,
                label="Matched",
            )
        )

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.grid(which="both", ls=":", lw=0.5, color="grey", alpha=0.7)

    cbar = fig.colorbar(scatter, ax=ax, ticks=unique_times)
    cbar.set_label("Cycle Time [ms]")
    ax.legend(handles=legend_handles, loc="best")
    return fig, ax


def plot_beta_variation(
    reference_twiss_df,
    comparison_twiss_df,
    harmonics=None,
    brho_Tm=DEFAULT_INJECTION_BRHO_TM,
    hqcal=DEFAULT_HQCAL,
    tqgcal=DEFAULT_TQ_GCAL,
    smooth_points=1001,
    ax=None,
):
    """
    Plot beta functions, local beta-function changes and optional harmonic wave.
    """

    for label, twiss_df in (("reference", reference_twiss_df), ("comparison", comparison_twiss_df)):
        missing = [column for column in ("s", "betx", "bety") if column not in twiss_df.columns]
        if missing:
            raise ValueError(f"{label} TWISS DataFrame is missing columns: {missing}")

    include_harmonics = harmonics is not None

    if ax is None:
        if include_harmonics:
            fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True, tight_layout=True)
            ax_betx, ax_bety, ax_delta, ax_wave = axes
        else:
            fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True, tight_layout=True)
            ax_betx, ax_bety, ax_delta = axes
            ax_wave = None
    else:
        if include_harmonics:
            ax_betx, ax_bety, ax_delta, ax_wave = ax
            axes = (ax_betx, ax_bety, ax_delta, ax_wave)
        else:
            ax_betx, ax_bety, ax_delta = ax
            ax_wave = None
            axes = (ax_betx, ax_bety, ax_delta)
        fig = ax_betx.figure

    s = reference_twiss_df["s"].astype(float)
    ax_betx.plot(s, reference_twiss_df["betx"].astype(float), label=r"Reference $\beta_x$")
    ax_betx.plot(s, comparison_twiss_df["betx"].astype(float), label=r"Comparison $\beta_x$", ls="--")
    ax_betx.set_ylabel(r"$\beta_x$ [m]")
    ax_betx.grid(ls=":", lw=0.5, color="grey", which="both")
    ax_betx.legend()

    ax_bety.plot(s, reference_twiss_df["bety"].astype(float), label=r"Reference $\beta_y$")
    ax_bety.plot(s, comparison_twiss_df["bety"].astype(float), label=r"Comparison $\beta_y$", ls="--")
    ax_bety.set_ylabel(r"$\beta_y$ [m]")
    ax_bety.grid(ls=":", lw=0.5, color="grey", which="both")
    ax_bety.legend()

    ax_delta.plot(s, reference_twiss_df["betx"].astype(float) - comparison_twiss_df["betx"].astype(float), label=r"delta $\beta_x$")
    ax_delta.plot(s, reference_twiss_df["bety"].astype(float) - comparison_twiss_df["bety"].astype(float), label=r"delta $\beta_y$")
    if not include_harmonics:
        ax_delta.set_xlabel("S [m]")
    ax_delta.set_ylabel(r"delta $\beta$ [m]")
    ax_delta.grid(ls=":", lw=0.5, color="grey", which="both")
    ax_delta.legend()

    if include_harmonics:
        circumference = float(np.nanmax(s))
        trim_positions = _trim_quad_positions_by_family(comparison_twiss_df)
        if trim_positions:
            if "qtd" in trim_positions:
                _plot_family_wave_on_s_axis(
                    ax_wave,
                    trim_positions["qtd"],
                    harmonics=harmonics,
                    family="qtd",
                    color="tab:blue",
                    circumference=circumference,
                    brho_Tm=brho_Tm,
                    hqcal=hqcal,
                    tqgcal=tqgcal,
                    smooth_points=smooth_points,
                )
            if "qtf" in trim_positions:
                _plot_family_wave_on_s_axis(
                    ax_wave,
                    trim_positions["qtf"],
                    harmonics=harmonics,
                    family="qtf",
                    color="tab:orange",
                    circumference=circumference,
                    brho_Tm=brho_Tm,
                    hqcal=hqcal,
                    tqgcal=tqgcal,
                    smooth_points=smooth_points,
                )
        else:
            smooth_s = np.linspace(float(np.nanmin(s)), circumference, int(smooth_points))
            smooth_superperiod = smooth_s / circumference * 10.0
            _, qtd_wave, qtf_wave = _harmonic_wave(
                harmonics=harmonics,
                superperiod_coordinate=smooth_superperiod,
                brho_Tm=brho_Tm,
                hqcal=hqcal,
                tqgcal=tqgcal,
            )
            ax_wave.plot(smooth_s, qtd_wave, color="tab:blue", lw=1.2, label="Expected QTD K")
            ax_wave.plot(smooth_s, qtf_wave, color="tab:orange", lw=1.2, label="Expected QTF K")

        ax_wave.axhline(0.0, color="black", lw=0.8)
        ax_wave.set_xlabel("S [m]")
        ax_wave.set_ylabel("Expected K [m^-2]")
        ax_wave.grid(ls=":", lw=0.5, color="grey", which="both")
        ax_wave.legend(loc="best")

    return fig, tuple(axes)


def plot_harmonic_trim_quad_pattern(
    expected_df,
    programmed_df=None,
    value_column="delta_k",
    harmonics=None,
    brho_Tm=DEFAULT_INJECTION_BRHO_TM,
    hqcal=DEFAULT_HQCAL,
    tqgcal=DEFAULT_TQ_GCAL,
    smooth_points=1001,
    show_phase_axis=True,
    shared_ylim=True,
    zero_tolerance=1.0e-12,
    ax=None,
):
    """
    Plot expected harmonic K offsets and optional programmed values.
    """

    required = ("superperiod", "qtd_expected_k", "qtf_expected_k")
    missing = [column for column in required if column not in expected_df.columns]
    if missing:
        raise ValueError(f"Expected pattern DataFrame is missing columns: {missing}")

    if ax is None:
        fig, (ax_qtd, ax_qtf) = plt.subplots(2, 1, figsize=(10, 7), sharex=True, tight_layout=True)
    else:
        ax_qtd, ax_qtf = ax
        fig = ax_qtd.figure

    superperiod = expected_df["superperiod"].astype(int).to_numpy()
    dense_superperiod = np.linspace(0.0, 10.0, int(smooth_points))
    if harmonics is not None:
        _, dense_qtd, dense_qtf = _harmonic_wave(
            harmonics=harmonics,
            superperiod_coordinate=dense_superperiod,
            brho_Tm=brho_Tm,
            hqcal=hqcal,
            tqgcal=tqgcal,
        )
    else:
        dense_qtd = np.interp(
            dense_superperiod,
            np.append(superperiod, 10.0),
            np.append(expected_df["qtd_expected_k"].astype(float).to_numpy(), expected_df["qtd_expected_k"].astype(float).to_numpy()[0]),
        )
        dense_qtf = np.interp(
            dense_superperiod,
            np.append(superperiod, 10.0),
            np.append(expected_df["qtf_expected_k"].astype(float).to_numpy(), expected_df["qtf_expected_k"].astype(float).to_numpy()[0]),
        )

    width = 0.35

    for axis, family, column, color in (
        (ax_qtd, "qtd", "qtd_expected_k", "tab:blue"),
        (ax_qtf, "qtf", "qtf_expected_k", "tab:orange"),
    ):
        dense_values = dense_qtd if family == "qtd" else dense_qtf
        axis.plot(dense_superperiod, dense_values, color=color, lw=1.2, label="Expected wave")
        axis.plot(
            superperiod,
            expected_df[column].astype(float),
            marker="o",
            ls="None",
            color=color,
            label="Expected at trim quads",
        )

        if programmed_df is not None:
            family_df = programmed_df[programmed_df["family"].astype(str).str.lower() == family]
            if value_column not in family_df.columns:
                raise ValueError(f"Programmed DataFrame is missing column: {value_column}")
            axis.bar(
                family_df["superperiod"].astype(int).to_numpy(),
                family_df[value_column].astype(float).to_numpy(),
                width=width,
                alpha=0.45,
                color=color,
                label=f"Programmed {value_column}",
            )

        axis.axhline(0.0, color="black", lw=0.8)
        axis.set_ylabel(f"{family.upper()} K [m^-2]")
        axis.set_xlim(-0.2, 10.2)
        axis.grid(ls=":", lw=0.5, color="grey", which="both")
        axis.legend(loc="best")

    ax_qtf.set_xlabel("Superperiod")
    ax_qtf.set_xticks(np.arange(0, 11, 1))
    if show_phase_axis:
        secax = ax_qtd.secondary_xaxis(
            "top",
            functions=(lambda x: x * 36.0, lambda degrees: degrees / 36.0),
        )
        secax.set_xlabel("Ring phase [deg]")

    if shared_ylim:
        y_values = [
            np.asarray(dense_qtd, dtype=float),
            np.asarray(dense_qtf, dtype=float),
            expected_df["qtd_expected_k"].astype(float).to_numpy(),
            expected_df["qtf_expected_k"].astype(float).to_numpy(),
        ]
        if programmed_df is not None and value_column in programmed_df.columns:
            y_values.append(programmed_df[value_column].astype(float).to_numpy())

        finite_values = np.concatenate([values[np.isfinite(values)] for values in y_values])
        finite_values[np.abs(finite_values) < float(zero_tolerance)] = 0.0
        max_abs = float(np.max(np.abs(finite_values))) if len(finite_values) else 0.0
        if max_abs < float(zero_tolerance):
            max_abs = float(zero_tolerance)

        margin = 1.1
        ax_qtd.set_ylim(-margin * max_abs, margin * max_abs)
        ax_qtf.set_ylim(-margin * max_abs, margin * max_abs)

    return fig, (ax_qtd, ax_qtf)


def plot_trim_quad_currents(currents, ax=None):
    """
    Plot QTF/QTD trim quad currents from a result or mapping.
    """

    if hasattr(currents, "trim_quad_display_dict"):
        values = currents.trim_quad_display_dict()
    else:
        values = dict(currents)

    labels = ["iqtf_A", "iqtd_A"]
    amps = [float(values[label]) for label in labels]

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 3), tight_layout=True)
    else:
        fig = ax.figure

    ax.bar(labels, amps, color=["tab:orange", "tab:blue"])
    ax.set_ylabel("Current [A]")
    ax.grid(axis="y", ls=":", lw=0.5, color="grey")
    return fig, ax
