"""
Durable run-bundle helpers for snapshot and series results.
"""

from dataclasses import dataclass, field
import json
from pathlib import Path

import pandas as pd

from .. import __version__
from ..snapshot import SnapshotResult, SnapshotSeriesResult
from .configs import config_to_record, series_config_to_record


@dataclass
class ArchivedRun:
    """
    Lightweight representation of a saved run bundle.
    """

    root_dir: str
    manifest: dict
    metadata: dict = field(default_factory=dict)
    tables: dict = field(default_factory=dict)

    def available_tables(self):
        return list(self.tables)

    def table(self, name):
        key = str(name)
        if key not in self.tables:
            raise KeyError(f"Unknown archived table {name!r}. Available tables: {self.available_tables()}")
        return self.tables[key].copy()


def _json_safe(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    return str(value)


def _write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_json_safe(data), indent=2, sort_keys=True))
    return str(path)


def _write_tables(result, bundle_dir):
    tables_dir = Path(bundle_dir) / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    table_entries = []
    for table_name in result.available_tables():
        table = result.table(table_name)
        if not isinstance(table, pd.DataFrame):
            continue
        path = tables_dir / f"{table_name}.csv"
        table.to_csv(path, index=False)
        table_entries.append(
            {
                "name": table_name,
                "path": str(path.relative_to(bundle_dir)),
                "format": "csv",
                "rows": int(len(table)),
                "columns": list(table.columns),
            }
        )
    return table_entries


def _bundle_dir_for_result(result, output_dir):
    if output_dir is not None:
        return Path(output_dir)
    if getattr(result, "run_paths", None) is not None:
        run_paths = result.run_paths
        if getattr(run_paths, "snapshot_dir", None):
            return Path(run_paths.snapshot_dir) / "bundle"
        if getattr(run_paths, "series_dir", None):
            return Path(run_paths.series_dir) / "bundle"
    raise ValueError("output_dir is required when result.run_paths is unavailable.")


def write_snapshot_bundle(result, output_dir=None):
    """
    Save a SnapshotResult's config, metadata, manifest and available tables.
    """

    if not isinstance(result, SnapshotResult):
        raise TypeError("result must be a SnapshotResult.")
    bundle_dir = _bundle_dir_for_result(result, output_dir)
    bundle_dir.mkdir(parents=True, exist_ok=True)
    table_entries = _write_tables(result, bundle_dir)
    _write_json(bundle_dir / "config.json", config_to_record(result.config))
    _write_json(bundle_dir / "metadata.json", result.metadata)
    manifest = {
        "bundle_type": "snapshot",
        "package_version": __version__,
        "snapshot_id": result.snapshot_id,
        "label": result.label,
        "case": result.case,
        "tables": table_entries,
        "metadata_path": "metadata.json",
        "config_path": "config.json",
    }
    _write_json(bundle_dir / "run_manifest.json", manifest)
    return str(bundle_dir)


def write_series_bundle(result, output_dir=None):
    """
    Save a SnapshotSeriesResult's config, metadata, manifest and aggregate tables.
    """

    if not isinstance(result, SnapshotSeriesResult):
        raise TypeError("result must be a SnapshotSeriesResult.")
    bundle_dir = _bundle_dir_for_result(result, output_dir)
    bundle_dir.mkdir(parents=True, exist_ok=True)
    table_entries = _write_tables(result, bundle_dir)
    _write_json(bundle_dir / "config.json", series_config_to_record(result.config))
    _write_json(bundle_dir / "metadata.json", result.metadata)
    manifest = {
        "bundle_type": "snapshot_series",
        "package_version": __version__,
        "label": result.config.label,
        "n_snapshots": len(result.snapshots),
        "tables": table_entries,
        "metadata_path": "metadata.json",
        "config_path": "config.json",
    }
    _write_json(bundle_dir / "run_manifest.json", manifest)
    return str(bundle_dir)


def read_run_bundle(bundle_dir):
    """
    Read a saved run bundle without rerunning MAD-X.
    """

    bundle_dir = Path(bundle_dir)
    manifest = json.loads((bundle_dir / "run_manifest.json").read_text())
    metadata_path = manifest.get("metadata_path", "metadata.json")
    metadata = json.loads((bundle_dir / metadata_path).read_text())
    tables = {}
    for entry in manifest.get("tables", []):
        if entry.get("format") != "csv":
            raise ValueError(f"Unsupported table format {entry.get('format')!r}.")
        tables[entry["name"]] = pd.read_csv(bundle_dir / entry["path"])
    return ArchivedRun(
        root_dir=str(bundle_dir),
        manifest=manifest,
        metadata=metadata,
        tables=tables,
    )
