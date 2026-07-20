"""
envelope.py

Beam-envelope utilities for the ISIS RCS optics GUI backend.

This layer consumes real MAD-X TWISS tables and explicit beam-size assumptions.
It does not generate optics data itself; MAD-X execution stays in MadxModel or
workflow wrappers.
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .madx_model import MadxModel


DEFAULT_EMITTANCE_PI_MM_MRAD = 300.0
DEFAULT_SIGMA_SCALE = 3.0
DEFAULT_DP_OVER_P = 0.002

GEOMETRIC_EMITTANCE_MODES = {"geometric", "geometric_rms"}
NORMALISED_EMITTANCE_MODES = {
    "normalised",
    "normalised_rms",
    "normalized",
    "normalized_rms",
}
SUPPORTED_EMITTANCE_MODES = GEOMETRIC_EMITTANCE_MODES | NORMALISED_EMITTANCE_MODES

REQUIRED_TWISS_COLUMNS = ("name", "s", "betx", "bety", "x", "y", "dx", "dy")


@dataclass
class EnvelopeInputs:
    """
    User-facing beam-envelope assumptions.

    Emittances are accepted using the accelerator convention "pi mm mrad".
    The numerical conversion used here is value * 1e-6 m rad; the pi is treated
    as part of the unit convention, not as an extra numerical multiplier.
    """

    emit_x_pi_mm_mrad: float = DEFAULT_EMITTANCE_PI_MM_MRAD
    emit_y_pi_mm_mrad: float = DEFAULT_EMITTANCE_PI_MM_MRAD
    emittance_mode: str = "geometric"
    sigma_scale: float = DEFAULT_SIGMA_SCALE
    dp_over_p: float = DEFAULT_DP_OVER_P
    label: str = "envelope"
    metadata: dict = field(default_factory=dict)

    def validated_mode(self):
        mode = str(self.emittance_mode).strip().lower()
        if mode not in SUPPORTED_EMITTANCE_MODES:
            valid = ", ".join(sorted(SUPPORTED_EMITTANCE_MODES))
            raise ValueError(f"Unsupported emittance_mode {self.emittance_mode!r}. Valid modes: {valid}")
        return mode

    def validate(self):
        self.validated_mode()
        _require_positive("emit_x_pi_mm_mrad", self.emit_x_pi_mm_mrad)
        _require_positive("emit_y_pi_mm_mrad", self.emit_y_pi_mm_mrad)
        _require_positive("sigma_scale", self.sigma_scale)
        if not np.isfinite(float(self.dp_over_p)):
            raise ValueError("dp_over_p must be finite.")
        return self

    @property
    def emit_x_m_rad_input(self):
        return pi_mm_mrad_to_m_rad(self.emit_x_pi_mm_mrad)

    @property
    def emit_y_m_rad_input(self):
        return pi_mm_mrad_to_m_rad(self.emit_y_pi_mm_mrad)


@dataclass
class EnvelopeResult:
    """
    Result container for one envelope evaluation.
    """

    inputs: EnvelopeInputs
    envelope_df: pd.DataFrame
    summary_df: pd.DataFrame
    beam_state: object = None
    source: str = "twiss"
    metadata: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def to_summary_dataframe(self):
        return self.summary_df.copy()

    def to_plot_dataframe(self):
        return self.envelope_df.copy()


def pi_mm_mrad_to_m_rad(value):
    """
    Convert a value quoted in pi mm mrad to m rad.

    ISIS operational notation normally treats "pi mm mrad" as the unit
    convention. Therefore 300 pi mm mrad is represented internally as 300e-6
    m rad, not 300 * pi * 1e-6.
    """

    value = float(value)
    if not np.isfinite(value):
        raise ValueError("Emittance must be finite.")
    return value * 1.0e-6


def _require_positive(name, value):
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite.")
    return value


def _require_beam_state_for_normalised(inputs, beam_state):
    mode = inputs.validated_mode()
    if mode in NORMALISED_EMITTANCE_MODES and beam_state is None:
        raise ValueError(
            "Normalised emittance mode requires a beam_state with beta and gamma."
        )


def _beam_beta_gamma(beam_state):
    try:
        beta = float(beam_state.beta)
        gamma = float(beam_state.gamma)
    except AttributeError as exc:
        raise ValueError("beam_state must expose beta and gamma attributes.") from exc

    _require_positive("beam_state.beta", beta)
    _require_positive("beam_state.gamma", gamma)
    return beta, gamma


def geometric_emittances_m_rad(inputs, beam_state=None):
    """
    Return horizontal and vertical geometric RMS emittances in m rad.
    """

    inputs.validate()
    mode = inputs.validated_mode()
    emit_x = inputs.emit_x_m_rad_input
    emit_y = inputs.emit_y_m_rad_input

    if mode in NORMALISED_EMITTANCE_MODES:
        _require_beam_state_for_normalised(inputs, beam_state)
        beta, gamma = _beam_beta_gamma(beam_state)
        scale = beta * gamma
        emit_x /= scale
        emit_y /= scale

    return float(emit_x), float(emit_y)


def validate_twiss_for_envelope(twiss_df):
    """
    Return a cleaned TWISS table after checking required envelope columns.
    """

    if twiss_df is None:
        raise ValueError("twiss_df is required.")

    df = twiss_df.copy()
    df.columns = [str(col).lower() for col in df.columns]

    missing = [col for col in REQUIRED_TWISS_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"TWISS table is missing required envelope columns: {missing}")

    for col in ("s", "betx", "bety", "x", "y", "dx", "dy"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    numeric_cols = ["s", "betx", "bety", "x", "y", "dx", "dy"]
    bad_cols = [col for col in numeric_cols if not np.isfinite(df[col]).all()]
    if bad_cols:
        raise ValueError(f"TWISS table contains non-finite values in columns: {bad_cols}")

    if (df["betx"] < 0.0).any() or (df["bety"] < 0.0).any():
        raise ValueError("TWISS beta functions must be non-negative.")

    return df


def evaluate_envelope_from_twiss(
    twiss_df,
    inputs=None,
    beam_state=None,
    source="twiss",
    metadata=None,
):
    """
    Add beam-envelope columns to a real MAD-X TWISS DataFrame.
    """

    if inputs is None:
        inputs = EnvelopeInputs()
    inputs.validate()
    _require_beam_state_for_normalised(inputs, beam_state)

    df = validate_twiss_for_envelope(twiss_df)
    emit_x_m_rad, emit_y_m_rad = geometric_emittances_m_rad(inputs, beam_state=beam_state)

    sigma_scale = float(inputs.sigma_scale)
    dp_over_p = float(inputs.dp_over_p)

    sigma_x_beta_m = np.sqrt(df["betx"] * emit_x_m_rad)
    sigma_y_beta_m = np.sqrt(df["bety"] * emit_y_m_rad)
    sigma_x_disp_m = np.abs(df["dx"] * dp_over_p)
    sigma_y_disp_m = np.abs(df["dy"] * dp_over_p)

    df["sigma_x_beta_m"] = sigma_x_beta_m
    df["sigma_y_beta_m"] = sigma_y_beta_m
    df["sigma_x_disp_m"] = sigma_x_disp_m
    df["sigma_y_disp_m"] = sigma_y_disp_m
    df["sigma_x_m"] = np.sqrt(sigma_x_beta_m**2 + sigma_x_disp_m**2)
    df["sigma_y_m"] = np.sqrt(sigma_y_beta_m**2 + sigma_y_disp_m**2)
    df["sigma_x_mm"] = df["sigma_x_m"] * 1.0e3
    df["sigma_y_mm"] = df["sigma_y_m"] * 1.0e3

    df["orbit_x_mm"] = df["x"] * 1.0e3
    df["orbit_y_mm"] = df["y"] * 1.0e3
    df["dispersion_x_mm"] = df["dx"] * dp_over_p * 1.0e3
    df["dispersion_y_mm"] = df["dy"] * dp_over_p * 1.0e3

    df["envelope_x_plus_m"] = df["x"] + sigma_scale * df["sigma_x_m"]
    df["envelope_x_minus_m"] = df["x"] - sigma_scale * df["sigma_x_m"]
    df["envelope_y_plus_m"] = df["y"] + sigma_scale * df["sigma_y_m"]
    df["envelope_y_minus_m"] = df["y"] - sigma_scale * df["sigma_y_m"]

    envelope_mm_columns = {
        "envelope_x_plus_m": "envelope_x_plus_mm",
        "envelope_x_minus_m": "envelope_x_minus_mm",
        "envelope_y_plus_m": "envelope_y_plus_mm",
        "envelope_y_minus_m": "envelope_y_minus_mm",
    }
    for metres_col, millimetres_col in envelope_mm_columns.items():
        df[millimetres_col] = df[metres_col] * 1.0e3

    summary_df = summarise_envelope(df, inputs, emit_x_m_rad, emit_y_m_rad)

    result_metadata = {
        "source": source,
        "n_rows": len(df),
        "emittance_mode": inputs.validated_mode(),
        "emit_x_geometric_m_rad": emit_x_m_rad,
        "emit_y_geometric_m_rad": emit_y_m_rad,
        "sigma_scale": sigma_scale,
        "dp_over_p": dp_over_p,
    }
    if metadata:
        result_metadata.update(dict(metadata))

    return EnvelopeResult(
        inputs=inputs,
        envelope_df=df,
        summary_df=summary_df,
        beam_state=beam_state,
        source=source,
        metadata=result_metadata,
        warnings=[],
    )


def summarise_envelope(envelope_df, inputs, emit_x_m_rad, emit_y_m_rad):
    """
    Build a compact summary table for notebook and GUI display.
    """

    df = envelope_df
    rows = []
    definitions = [
        ("x", "sigma_x_mm", "envelope_x_plus_mm", "envelope_x_minus_mm"),
        ("y", "sigma_y_mm", "envelope_y_plus_mm", "envelope_y_minus_mm"),
    ]

    for plane, sigma_col, plus_col, minus_col in definitions:
        sigma_idx = df[sigma_col].idxmax()
        plus_idx = df[plus_col].idxmax()
        minus_idx = df[minus_col].idxmin()
        rows.append(
            {
                "label": inputs.label,
                "plane": plane,
                "sigma_scale": float(inputs.sigma_scale),
                "dp_over_p": float(inputs.dp_over_p),
                "emittance_mode": inputs.validated_mode(),
                "emit_x_geometric_m_rad": float(emit_x_m_rad),
                "emit_y_geometric_m_rad": float(emit_y_m_rad),
                "max_sigma_mm": float(df.loc[sigma_idx, sigma_col]),
                "max_sigma_name": df.loc[sigma_idx, "name"],
                "max_sigma_s_m": float(df.loc[sigma_idx, "s"]),
                "max_plus_mm": float(df.loc[plus_idx, plus_col]),
                "max_plus_name": df.loc[plus_idx, "name"],
                "max_plus_s_m": float(df.loc[plus_idx, "s"]),
                "min_minus_mm": float(df.loc[minus_idx, minus_col]),
                "min_minus_name": df.loc[minus_idx, "name"],
                "min_minus_s_m": float(df.loc[minus_idx, "s"]),
            }
        )

    return pd.DataFrame(rows)


def evaluate_envelope_from_machine_state(
    machine_state,
    lattice_folder="../Lattice_Files/00_Simplified_Lattice",
    output_dir="./envelope_tests/madx",
    sequence_name="synchrotron",
    inputs=None,
    twiss_columns=None,
):
    """
    Run MAD-X for a MachineState and evaluate the resulting envelope.
    """

    model = MadxModel(
        lattice_folder=lattice_folder,
        sequence_name=sequence_name,
        output_dir=output_dir,
    )
    model.load_lattice(use_sequence=False)
    model.apply_machine_state(machine_state)
    model.use_sequence()
    twiss_df = model.run_twiss(columns=twiss_columns)

    metadata = model.get_metadata()
    metadata["output_dir"] = str(Path(output_dir))

    result = evaluate_envelope_from_twiss(
        twiss_df,
        inputs=inputs,
        beam_state=machine_state.beam_state,
        source="machine_state",
        metadata=metadata,
    )
    result.metadata["madx_summary"] = model.get_summary_dict()
    return result


def compare_envelopes(results):
    """
    Return a compact comparison table for multiple EnvelopeResult objects.
    """

    rows = []
    for result in results:
        summary = result.summary_df.copy()
        rows.append(summary)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def plot_envelope(result, plane="x", ax=None, title=None):
    """
    Plot orbit and plus/minus envelope in one plane.
    """

    import matplotlib.pyplot as plt

    plane = str(plane).lower()
    if plane not in {"x", "y"}:
        raise ValueError("plane must be 'x' or 'y'.")

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    df = result.envelope_df
    s = df["s"]
    orbit_col = f"orbit_{plane}_mm"
    plus_col = f"envelope_{plane}_plus_mm"
    minus_col = f"envelope_{plane}_minus_mm"

    envelope_color = "tab:blue"
    envelope_linewidth = 1.0
    envelope_linestyle = "-"

    ax.plot(s, df[orbit_col], color="black", linewidth=1.4, label=f"{plane} orbit")
    ax.plot(
        s,
        df[plus_col],
        color=envelope_color,
        linewidth=envelope_linewidth,
        linestyle=envelope_linestyle,
        label=f"+{result.inputs.sigma_scale:g} sigma",
    )
    ax.plot(
        s,
        df[minus_col],
        color=envelope_color,
        linewidth=envelope_linewidth,
        linestyle=envelope_linestyle,
        label=f"-{result.inputs.sigma_scale:g} sigma",
    )
    ax.fill_between(s, df[minus_col], df[plus_col], color=envelope_color, alpha=0.15)
    ax.set_xlabel("s [m]")
    ax.set_ylabel(f"{plane} [mm]")
    ax.set_title(title or f"{result.inputs.label}: {plane.upper()} envelope")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return ax


def plot_sigma(result, ax=None, title=None):
    """
    Plot horizontal and vertical RMS beam size.
    """

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    df = result.envelope_df
    ax.plot(df["s"], df["sigma_x_mm"], label="sigma x", color="tab:blue")
    ax.plot(df["s"], df["sigma_y_mm"], label="sigma y", color="tab:orange")
    ax.set_xlabel("s [m]")
    ax.set_ylabel("sigma [mm]")
    ax.set_title(title or f"{result.inputs.label}: RMS beam size")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return ax


def plot_envelope_comparison(results, plane="x", ax=None, title=None):
    """
    Plot plus/minus envelopes from several results on one axis.
    """

    import matplotlib.pyplot as plt

    plane = str(plane).lower()
    if plane not in {"x", "y"}:
        raise ValueError("plane must be 'x' or 'y'.")

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 4))

    for result in results:
        df = result.envelope_df
        label = result.inputs.label
        plus_line = ax.plot(
            df["s"],
            df[f"envelope_{plane}_plus_mm"],
            linewidth=1.0,
            linestyle="-",
            label=f"{label} +",
        )[0]
        ax.plot(
            df["s"],
            df[f"envelope_{plane}_minus_mm"],
            color=plus_line.get_color(),
            linewidth=plus_line.get_linewidth(),
            linestyle=plus_line.get_linestyle(),
            label=f"{label} -",
        )

    ax.set_xlabel("s [m]")
    ax.set_ylabel(f"{plane} envelope [mm]")
    ax.set_title(title or f"{plane.upper()} envelope comparison")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    return ax
