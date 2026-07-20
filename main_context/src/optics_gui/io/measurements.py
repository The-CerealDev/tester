"""
Adapters from table-shaped external inputs to backend-ready objects.
"""

from ..orbit_correction import normalise_bpm_measurements, normalise_corrector_selection
from ..snapshot import corrector_settings_from_dataframe, snapshot_configs_from_timepoint_table


def corrector_settings_from_table(table, cycle_time_ms=None, prefer="currents", source="epics_archiver"):
    return corrector_settings_from_dataframe(
        table,
        cycle_time_ms=cycle_time_ms,
        prefer=prefer,
        source=source,
    )


def snapshot_configs_from_table(table, base_config, corrector_table=None):
    return snapshot_configs_from_timepoint_table(
        table,
        base_config=base_config,
        corrector_table=corrector_table,
    )


def normalise_bpm_table(table, enabled_default=True):
    return normalise_bpm_measurements(table, enabled_default=enabled_default)


def normalise_corrector_table(table=None, plane=None, enabled_default=True):
    return normalise_corrector_selection(
        table,
        plane=plane,
        enabled_default=enabled_default,
    )
