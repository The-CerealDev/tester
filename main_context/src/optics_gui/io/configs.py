"""
JSON-compatible config records for the optics GUI backend.
"""

from dataclasses import fields, is_dataclass
import json
from pathlib import Path

import pandas as pd

from ..envelope import EnvelopeInputs
from ..orbit_branch import OrbitBranchConfig
from ..snapshot import (
    SnapshotConfig,
    SnapshotCorrectorSettings,
    SnapshotOrbitCorrectionConfig,
    SnapshotPlotSaveConfig,
    SnapshotSeriesConfig,
)


_TYPE_REGISTRY = {
    "EnvelopeInputs": EnvelopeInputs,
    "OrbitBranchConfig": OrbitBranchConfig,
    "SnapshotConfig": SnapshotConfig,
    "SnapshotCorrectorSettings": SnapshotCorrectorSettings,
    "SnapshotOrbitCorrectionConfig": SnapshotOrbitCorrectionConfig,
    "SnapshotPlotSaveConfig": SnapshotPlotSaveConfig,
    "SnapshotSeriesConfig": SnapshotSeriesConfig,
}


def _json_safe(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_json_safe(item) for item in value]
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, pd.DataFrame):
        raise TypeError("DataFrames must be saved as tables, not embedded in JSON configs.")
    if is_dataclass(value):
        record = {
            field.name: _json_safe(getattr(value, field.name))
            for field in fields(value)
        }
        record["__type__"] = type(value).__name__
        return record
    raise TypeError(f"Value of type {type(value).__name__} is not JSON-compatible.")


def _strip_type(record):
    out = dict(record)
    out.pop("__type__", None)
    return out


def _construct_typed(record):
    if not isinstance(record, dict):
        return record
    type_name = record.get("__type__")
    if type_name is None:
        return {key: _construct_typed(value) for key, value in record.items()}
    cls = _TYPE_REGISTRY.get(type_name)
    if cls is None:
        raise ValueError(f"Unknown config record type {type_name!r}.")
    values = {key: _construct_typed(value) for key, value in _strip_type(record).items()}
    return cls(**values)


def config_to_record(config):
    """
    Convert a SnapshotConfig to a JSON-compatible dictionary.
    """

    if not isinstance(config, SnapshotConfig):
        raise TypeError("config must be a SnapshotConfig.")
    return _json_safe(config)


def config_from_record(record):
    """
    Build a SnapshotConfig from a record produced by config_to_record().
    """

    if not isinstance(record, dict):
        raise TypeError("record must be a dictionary.")
    record = dict(record)
    record.setdefault("__type__", "SnapshotConfig")
    config = _construct_typed(record)
    if not isinstance(config, SnapshotConfig):
        raise TypeError("record did not describe a SnapshotConfig.")
    return config


def series_config_to_record(config):
    """
    Convert a SnapshotSeriesConfig to a JSON-compatible dictionary.
    """

    if not isinstance(config, SnapshotSeriesConfig):
        raise TypeError("config must be a SnapshotSeriesConfig.")
    return _json_safe(config)


def series_config_from_record(record):
    """
    Build a SnapshotSeriesConfig from a JSON-compatible record.
    """

    if not isinstance(record, dict):
        raise TypeError("record must be a dictionary.")
    record = dict(record)
    record.setdefault("__type__", "SnapshotSeriesConfig")
    config = _construct_typed(record)
    if not isinstance(config, SnapshotSeriesConfig):
        raise TypeError("record did not describe a SnapshotSeriesConfig.")
    return config


def write_snapshot_config(config, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config_to_record(config), indent=2, sort_keys=True))
    return str(path)


def read_snapshot_config(path):
    return config_from_record(json.loads(Path(path).read_text()))


def write_snapshot_series_config(config, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(series_config_to_record(config), indent=2, sort_keys=True))
    return str(path)


def read_snapshot_series_config(path):
    return series_config_from_record(json.loads(Path(path).read_text()))
