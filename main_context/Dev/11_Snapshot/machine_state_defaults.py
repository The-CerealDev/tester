"""
machine_state_defaults.py

Default values and standard variable names for the ISIS RCS machine-state layer.
"""

from collections import OrderedDict


# ----------------------------------------------------------------------
# Main magnet scaling presets
# ----------------------------------------------------------------------

RCS_BARE_SCALING = OrderedDict(
    [
        ("main_dipole_scale", 1.0),
        ("fringe_dipole_scale", 1.0),
        ("fringe_dipole_neg_scale", 1.0),
        ("qd_scale", 1.0),
        ("qdfr_scale", 1.0),
        ("qf_scale", 1.0),
        ("qffr_scale", 1.0),
        ("qds_scale", 1.0),
        ("qdsfr_scale", 1.0),
    ]
)

SRM_BARE_SCALING = OrderedDict(
    [
        ("main_dipole_scale", 1.0),
        ("fringe_dipole_scale", 1.0),
        ("fringe_dipole_neg_scale", 1.0),
        ("qd_scale", 0.984962103600024),
        ("qdfr_scale", 0.994661334722062),
        ("qf_scale", 0.9887119657096543),
        ("qffr_scale", 0.9959758637608666),
        ("qds_scale", 0.9931335289706212),
        ("qdsfr_scale", 0.9958189246720367),
    ]
)

MAIN_MAGNET_SCALING_PRESETS = {
    "rcs_bare": RCS_BARE_SCALING,
    "srm_bare": SRM_BARE_SCALING,
}


# ----------------------------------------------------------------------
# Tune defaults
# ----------------------------------------------------------------------

DEFAULT_QX = 4.331
DEFAULT_QY = 3.731

DEFAULT_BASE_QX = 4.331
DEFAULT_BASE_QY = 3.731

DEFAULT_TUNE_METHOD = "di_wright"

DEFAULT_KQTD = 0.0
DEFAULT_KQTF = 0.0

DEFAULT_TQ_GCAL = 1.997e-3

DEFAULT_DI_TUNE_COEFFICIENTS = [
    -4.73e-3,
    -5.99e-3,
    4.45e-3,
    2.40e-3,
]


# ----------------------------------------------------------------------
# Harmonic tune defaults
# ----------------------------------------------------------------------

DEFAULT_HARMONICS = OrderedDict(
    [
        ("D7SIN", 0.0),
        ("D7COS", 0.0),
        ("D8SIN", 0.0),
        ("D8COS", 0.0),
        ("F8SIN", 0.0),
        ("F8COS", 0.0),
        ("F9SIN", 0.0),
        ("F9COS", 0.0),
    ]
)


# ----------------------------------------------------------------------
# Operational HD/VD corrector names
# ----------------------------------------------------------------------

CORRECTOR_SUPERPERIODS = [0, 2, 3, 4, 5, 7, 9]

HD_CORRECTOR_NAMES = [
    "r0hd1_kick",
    "r2hd1_kick",
    "r3hd1_kick",
    "r4hd1_kick",
    "r5hd1_kick",
    "r7hd1_kick",
    "r9hd1_kick",
]

VD_CORRECTOR_NAMES = [
    "r0vd1_kick",
    "r2vd1_kick",
    "r3vd1_kick",
    "r4vd1_kick",
    "r5vd1_kick",
    "r7vd1_kick",
    "r9vd1_kick",
]

DEFAULT_HD_CORRECTOR_KICKS_RAD = OrderedDict(
    [(name, 0.0) for name in HD_CORRECTOR_NAMES]
)

DEFAULT_VD_CORRECTOR_KICKS_RAD = OrderedDict(
    [(name, 0.0) for name in VD_CORRECTOR_NAMES]
)

DEFAULT_HD_CORRECTOR_CURRENTS_A = OrderedDict(
    [(name, 0.0) for name in HD_CORRECTOR_NAMES]
)

DEFAULT_VD_CORRECTOR_CURRENTS_A = OrderedDict(
    [(name, 0.0) for name in VD_CORRECTOR_NAMES]
)


# ----------------------------------------------------------------------
# Corrector calibration
# ----------------------------------------------------------------------
# Units:
#   calibration value converts current [A] to integrated field/kick numerator
#   kick_mrad = current_A * calibration / brho_Tm

CORRECTOR_CALIBRATION = {
    "0H": 0.08350,
    "2H": 0.09121,
    "3H": 0.08000,
    "4H": 0.06600,
    "5H": 0.07780,
    "7H": 0.07580,
    "9H": 0.07660,
    "0V": 0.04620,
    "2V": 0.04330,
    "3V": 0.05210,
    "4V": 0.04770,
    "5V": 0.05400,
    "7V": 0.05220,
    "9V": 0.04510,
}
