from pathlib import Path

from optics_gui.aperture import read_source_aperture_csv
from optics_gui.io import config_from_record, config_to_record
from optics_gui.snapshot import SnapshotConfig


def test_snapshot_config_json_record_roundtrip():
    config = SnapshotConfig(
        cycle_time_ms=0.0,
        requested_qx=4.31,
        requested_qy=3.83,
        output_dir="./runs",
        run_envelope=False,
        run_aperture=False,
    )

    record = config_to_record(config)
    restored = config_from_record(record)

    assert restored.cycle_time_ms == config.cycle_time_ms
    assert restored.requested_qx == config.requested_qx
    assert restored.output_dir == config.output_dir


def test_packaged_source_aperture_data_loads():
    table = read_source_aperture_csv()

    assert not table.empty
    assert {"s", "aperture_x_m", "aperture_y_m"}.issubset(table.columns)


def test_repo_local_lattice_data_remains_external():
    assert Path("Dev/Lattice_Files").is_dir()
