"""
orbit_branch.py

Independent orbit branch runner for the ISIS RCS optics GUI backend.

This layer owns workflow-level branch execution and comparison. MAD-X execution
remains inside MadxModel.
"""

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np
import pandas as pd

from errors import combine_error_tables, write_error_table
from machine_state_writer import write_machine_state_file
from madx_model import MadxModel


ORBIT_COLUMNS = (
    "name",
    "keyword",
    "s",
    "x",
    "y",
    "px",
    "py",
)


@dataclass
class OrbitBranchConfig:
    """
    Configuration for one independent MAD-X orbit branch.
    """

    name: str
    lattice_folder: str
    sequence_name: str = "synchrotron"
    machine_state: object = None
    machine_state_file: str = None
    aperture_file: str = None
    error_table_paths: list = field(default_factory=list)
    output_dir: str = None
    twiss_columns: list = None
    metadata: dict = field(default_factory=dict)


@dataclass
class OrbitBranchResult:
    """
    Result from one orbit branch run.
    """

    name: str
    twiss_df: pd.DataFrame
    orbit_df: pd.DataFrame
    summary_df: pd.DataFrame
    metadata: dict
    machine_state_file: str = None
    error_table_paths: list = field(default_factory=list)


def _safe_name(name):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name).strip())
    safe = safe.strip("._")
    return safe or "orbit_branch"


def extract_orbit_df(twiss_df):
    """
    Extract GUI-ready closed-orbit columns from a TWISS DataFrame.
    """

    if not isinstance(twiss_df, pd.DataFrame):
        raise TypeError("twiss_df must be a pandas DataFrame.")

    missing = [column for column in ORBIT_COLUMNS if column not in twiss_df.columns]
    if missing:
        raise ValueError(f"TWISS DataFrame is missing orbit columns: {missing}")

    orbit_df = twiss_df.loc[:, ORBIT_COLUMNS].copy()

    for column in ("s", "x", "y", "px", "py"):
        orbit_df[column] = pd.to_numeric(orbit_df[column], errors="coerce").astype(float)
        if not np.isfinite(orbit_df[column].to_numpy()).all():
            raise ValueError(f"Orbit column {column!r} contains non-finite values.")

    orbit_df["x_mm"] = 1.0e3 * orbit_df["x"]
    orbit_df["y_mm"] = 1.0e3 * orbit_df["y"]
    orbit_df["orbit_radius_m"] = np.sqrt(
        orbit_df["x"] * orbit_df["x"] + orbit_df["y"] * orbit_df["y"]
    )
    orbit_df["orbit_radius_mm"] = 1.0e3 * orbit_df["orbit_radius_m"]

    return orbit_df


def summarise_orbit_df(orbit_df):
    """
    Return compact orbit metrics suitable for notebook and GUI display.
    """

    if not isinstance(orbit_df, pd.DataFrame):
        raise TypeError("orbit_df must be a pandas DataFrame.")

    required = ("x", "y", "orbit_radius_m")
    missing = [column for column in required if column not in orbit_df.columns]
    if missing:
        raise ValueError(f"Orbit DataFrame is missing summary columns: {missing}")

    x = orbit_df["x"].astype(float).to_numpy()
    y = orbit_df["y"].astype(float).to_numpy()
    radius = orbit_df["orbit_radius_m"].astype(float).to_numpy()

    return {
        "n_rows": int(len(orbit_df)),
        "max_abs_x_m": float(np.nanmax(np.abs(x))) if len(x) else 0.0,
        "max_abs_y_m": float(np.nanmax(np.abs(y))) if len(y) else 0.0,
        "rms_x_m": float(np.sqrt(np.nanmean(x * x))) if len(x) else 0.0,
        "rms_y_m": float(np.sqrt(np.nanmean(y * y))) if len(y) else 0.0,
        "max_orbit_radius_m": float(np.nanmax(radius)) if len(radius) else 0.0,
        "max_abs_x_mm": float(1.0e3 * np.nanmax(np.abs(x))) if len(x) else 0.0,
        "max_abs_y_mm": float(1.0e3 * np.nanmax(np.abs(y))) if len(y) else 0.0,
        "rms_x_mm": float(1.0e3 * np.sqrt(np.nanmean(x * x))) if len(x) else 0.0,
        "rms_y_mm": float(1.0e3 * np.sqrt(np.nanmean(y * y))) if len(y) else 0.0,
        "max_orbit_radius_mm": float(1.0e3 * np.nanmax(radius)) if len(radius) else 0.0,
    }


def compare_orbit_results(reference, other):
    """
    Compare two OrbitBranchResult objects or orbit DataFrames.
    """

    reference_orbit = reference.orbit_df if isinstance(reference, OrbitBranchResult) else reference
    other_orbit = other.orbit_df if isinstance(other, OrbitBranchResult) else other

    if not isinstance(reference_orbit, pd.DataFrame):
        raise TypeError("reference must be an OrbitBranchResult or pandas DataFrame.")
    if not isinstance(other_orbit, pd.DataFrame):
        raise TypeError("other must be an OrbitBranchResult or pandas DataFrame.")

    for label, orbit_df in (("reference", reference_orbit), ("other", other_orbit)):
        missing = [column for column in ("name", "s", "x", "y") if column not in orbit_df.columns]
        if missing:
            raise ValueError(f"{label} orbit DataFrame is missing columns: {missing}")

    if len(reference_orbit) != len(other_orbit):
        raise ValueError("Orbit DataFrames have different row counts.")

    reference_names = reference_orbit["name"].astype(str).reset_index(drop=True)
    other_names = other_orbit["name"].astype(str).reset_index(drop=True)
    if not reference_names.equals(other_names):
        raise ValueError("Orbit DataFrames have different row order or element names.")

    comparison = pd.DataFrame(
        {
            "name": reference_names,
            "keyword": reference_orbit.get("keyword", pd.Series(index=reference_orbit.index)).to_numpy(),
            "s": reference_orbit["s"].astype(float).to_numpy(),
            "x_reference": reference_orbit["x"].astype(float).to_numpy(),
            "y_reference": reference_orbit["y"].astype(float).to_numpy(),
            "x_other": other_orbit["x"].astype(float).to_numpy(),
            "y_other": other_orbit["y"].astype(float).to_numpy(),
        }
    )

    comparison["delta_x"] = comparison["x_other"] - comparison["x_reference"]
    comparison["delta_y"] = comparison["y_other"] - comparison["y_reference"]
    comparison["delta_orbit_radius_m"] = np.sqrt(
        comparison["delta_x"] * comparison["delta_x"]
        + comparison["delta_y"] * comparison["delta_y"]
    )
    comparison["delta_x_mm"] = 1.0e3 * comparison["delta_x"]
    comparison["delta_y_mm"] = 1.0e3 * comparison["delta_y"]
    comparison["delta_orbit_radius_mm"] = 1.0e3 * comparison["delta_orbit_radius_m"]

    return comparison


def summarise_orbit_difference(comparison_df):
    """
    Return compact metrics for an orbit comparison DataFrame.
    """

    if not isinstance(comparison_df, pd.DataFrame):
        raise TypeError("comparison_df must be a pandas DataFrame.")

    required = ("delta_x", "delta_y", "delta_orbit_radius_m")
    missing = [column for column in required if column not in comparison_df.columns]
    if missing:
        raise ValueError(f"Comparison DataFrame is missing columns: {missing}")

    dx = comparison_df["delta_x"].astype(float).to_numpy()
    dy = comparison_df["delta_y"].astype(float).to_numpy()
    radius = comparison_df["delta_orbit_radius_m"].astype(float).to_numpy()

    return {
        "n_rows": int(len(comparison_df)),
        "max_abs_delta_x_m": float(np.nanmax(np.abs(dx))) if len(dx) else 0.0,
        "max_abs_delta_y_m": float(np.nanmax(np.abs(dy))) if len(dy) else 0.0,
        "rms_delta_x_m": float(np.sqrt(np.nanmean(dx * dx))) if len(dx) else 0.0,
        "rms_delta_y_m": float(np.sqrt(np.nanmean(dy * dy))) if len(dy) else 0.0,
        "max_delta_orbit_radius_m": float(np.nanmax(radius)) if len(radius) else 0.0,
        "max_abs_delta_x_mm": float(1.0e3 * np.nanmax(np.abs(dx))) if len(dx) else 0.0,
        "max_abs_delta_y_mm": float(1.0e3 * np.nanmax(np.abs(dy))) if len(dy) else 0.0,
        "rms_delta_x_mm": float(1.0e3 * np.sqrt(np.nanmean(dx * dx))) if len(dx) else 0.0,
        "rms_delta_y_mm": float(1.0e3 * np.sqrt(np.nanmean(dy * dy))) if len(dy) else 0.0,
        "max_delta_orbit_radius_mm": float(1.0e3 * np.nanmax(radius)) if len(radius) else 0.0,
    }


class OrbitBranch:
    """
    Independent MAD-X orbit branch.
    """

    def __init__(self, config):
        if not isinstance(config, OrbitBranchConfig):
            raise TypeError("config must be an OrbitBranchConfig.")

        self.config = config

    def _output_dir(self):
        if self.config.output_dir is not None:
            return Path(self.config.output_dir)
        return Path("./orbit_branch_runs") / _safe_name(self.config.name)

    def _machine_state_file(self, output_dir):
        if self.config.machine_state is not None:
            state_dir = output_dir / "machine_states"
            filename = f"machine_state_{_safe_name(self.config.name)}.strength"
            return write_machine_state_file(
                self.config.machine_state,
                output_dir=state_dir,
                filename=filename,
            )

        if self.config.machine_state_file is not None:
            return str(self.config.machine_state_file)

        return None

    def _error_table_paths(self, output_dir):
        error_table_paths = list(self.config.error_table_paths or [])

        if len(error_table_paths) <= 1:
            return [str(path) for path in error_table_paths]

        combined = combine_error_tables(error_table_paths)
        error_dir = output_dir / "error_tables"
        error_path = error_dir / f"combined_errors_{_safe_name(self.config.name)}.tfs"

        return [
            write_error_table(
                combined,
                error_path,
                table_name=f"COMBINED_ERRORS_{_safe_name(self.config.name)}".upper(),
            )
        ]

    def run(self):
        """
        Run this branch and return an OrbitBranchResult.
        """

        output_dir = self._output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        machine_state_file = self._machine_state_file(output_dir)
        error_table_paths = self._error_table_paths(output_dir)

        model = MadxModel(
            lattice_folder=self.config.lattice_folder,
            sequence_name=self.config.sequence_name,
            machine_state_file=machine_state_file,
            aperture_file=self.config.aperture_file,
            output_dir=output_dir / "madx",
        )

        model.load_lattice(use_sequence=True)
        for error_table_path in error_table_paths:
            model.apply_error_table(error_table_path)

        twiss_df = model.run_twiss(columns=self.config.twiss_columns)
        orbit_df = extract_orbit_df(twiss_df)
        summary_df = model.get_summary_df()

        metadata = dict(self.config.metadata or {})
        metadata.update(model.get_metadata())
        metadata.update(
            {
                "branch_name": self.config.name,
                "errors_applied": bool(error_table_paths),
                "orbit_summary": summarise_orbit_df(orbit_df),
                "machine_state_file": machine_state_file,
                "error_table_paths": list(error_table_paths),
                "output_dir": str(output_dir),
            }
        )

        return OrbitBranchResult(
            name=self.config.name,
            twiss_df=twiss_df,
            orbit_df=orbit_df,
            summary_df=summary_df,
            metadata=metadata,
            machine_state_file=machine_state_file,
            error_table_paths=list(error_table_paths),
        )
