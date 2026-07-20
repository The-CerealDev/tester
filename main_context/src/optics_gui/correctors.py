"""
correctors.py

HD/VD corrector helpers for the ISIS RCS machine-state layer.

The machine-state database stores both:
    - MAD-X kicks in radians
    - controls-style currents in amperes

Only the kicks are written as MAD-X assignments.
Currents are retained for GUI display and future correction suggestions.
"""

from collections import OrderedDict

from .machine_state_defaults import (
    CORRECTOR_CALIBRATION,
    HD_CORRECTOR_NAMES,
    VD_CORRECTOR_NAMES,
)


def default_hd_corrector_kicks_rad():
    return OrderedDict((name, 0.0) for name in HD_CORRECTOR_NAMES)


def default_vd_corrector_kicks_rad():
    return OrderedDict((name, 0.0) for name in VD_CORRECTOR_NAMES)


def default_hd_corrector_currents_A():
    return OrderedDict((name, 0.0) for name in HD_CORRECTOR_NAMES)


def default_vd_corrector_currents_A():
    return OrderedDict((name, 0.0) for name in VD_CORRECTOR_NAMES)


def default_corrector_kicks_rad():
    kicks = OrderedDict()
    kicks.update(default_hd_corrector_kicks_rad())
    kicks.update(default_vd_corrector_kicks_rad())
    return kicks


def default_corrector_currents_A():
    currents = OrderedDict()
    currents.update(default_hd_corrector_currents_A())
    currents.update(default_vd_corrector_currents_A())
    return currents


def _extract_superperiod_and_plane(corrector_name):
    """
    Extract superperiod and plane from names such as:
        r0hd1_kick
        r2vd1_kick
    """

    name = str(corrector_name).lower()

    if not name.startswith("r"):
        raise ValueError(f"Cannot parse corrector name: {corrector_name}")

    try:
        superperiod = int(name[1])
    except ValueError as exc:
        raise ValueError(f"Cannot parse superperiod from {corrector_name}") from exc

    if "hd" in name:
        plane = "H"
    elif "vd" in name:
        plane = "V"
    else:
        raise ValueError(f"Cannot parse plane from {corrector_name}")

    return superperiod, plane


def corrector_calibration_key(corrector_name):
    superperiod, plane = _extract_superperiod_and_plane(corrector_name)
    return f"{superperiod}{plane}"


def current_to_kick_rad(corrector_name, current_A, beam_state, calibration=None):
    """
    Convert one corrector current [A] to MAD-X kick [rad].
    """

    if calibration is None:
        calibration = CORRECTOR_CALIBRATION

    key = corrector_calibration_key(corrector_name)

    if key not in calibration:
        raise KeyError(f"No corrector calibration found for {corrector_name} ({key}).")

    brho_Tm = float(beam_state.brho_Tm)

    if brho_Tm == 0.0:
        raise ZeroDivisionError("beam_state.brho_Tm must be non-zero.")

    kick_mrad = float(current_A) * float(calibration[key]) / brho_Tm
    kick_rad = kick_mrad * 1.0e-3

    return float(kick_rad)


def kick_rad_to_current(corrector_name, kick_rad, beam_state, calibration=None):
    """
    Convert one MAD-X corrector kick [rad] to current [A].
    """

    if calibration is None:
        calibration = CORRECTOR_CALIBRATION

    key = corrector_calibration_key(corrector_name)

    if key not in calibration:
        raise KeyError(f"No corrector calibration found for {corrector_name} ({key}).")

    calibration_value = float(calibration[key])

    if calibration_value == 0.0:
        raise ZeroDivisionError(f"Calibration for {corrector_name} is zero.")

    brho_Tm = float(beam_state.brho_Tm)
    kick_mrad = float(kick_rad) * 1.0e3

    current_A = kick_mrad * brho_Tm / calibration_value

    return float(current_A)


def currents_to_kicks_rad(currents_A, beam_state, calibration=None):
    """
    Convert a dictionary of corrector currents [A] to kicks [rad].
    """

    return OrderedDict(
        (
            name,
            current_to_kick_rad(
                corrector_name=name,
                current_A=current,
                beam_state=beam_state,
                calibration=calibration,
            ),
        )
        for name, current in currents_A.items()
    )


def kicks_rad_to_currents_A(kicks_rad, beam_state, calibration=None):
    """
    Convert a dictionary of corrector kicks [rad] to currents [A].
    """

    return OrderedDict(
        (
            name,
            kick_rad_to_current(
                corrector_name=name,
                kick_rad=kick,
                beam_state=beam_state,
                calibration=calibration,
            ),
        )
        for name, kick in kicks_rad.items()
    )


def split_hd_vd(corrector_dict):
    """
    Split a corrector dictionary into HD and VD dictionaries.
    """

    hd = OrderedDict()
    vd = OrderedDict()

    for name, value in corrector_dict.items():
        lower_name = str(name).lower()

        if "hd" in lower_name:
            hd[name] = value
        elif "vd" in lower_name:
            vd[name] = value
        else:
            raise ValueError(f"Corrector name is neither HD nor VD: {name}")

    return hd, vd


def build_corrector_state(
    beam_state,
    hd_corrector_kicks_rad=None,
    vd_corrector_kicks_rad=None,
    hd_corrector_currents_A=None,
    vd_corrector_currents_A=None,
    prefer="kicks",
):
    """
    Build a complete corrector state.

    Parameters
    ----------
    prefer : {"kicks", "currents"}
        If both kicks and currents are supplied, this determines which values
        are treated as primary.
    """

    if prefer not in ("kicks", "currents"):
        raise ValueError("prefer must be either 'kicks' or 'currents'.")

    if hd_corrector_kicks_rad is None:
        hd_corrector_kicks_rad = default_hd_corrector_kicks_rad()

    if vd_corrector_kicks_rad is None:
        vd_corrector_kicks_rad = default_vd_corrector_kicks_rad()

    if hd_corrector_currents_A is None:
        hd_corrector_currents_A = default_hd_corrector_currents_A()

    if vd_corrector_currents_A is None:
        vd_corrector_currents_A = default_vd_corrector_currents_A()

    hd_corrector_kicks_rad = OrderedDict(hd_corrector_kicks_rad)
    vd_corrector_kicks_rad = OrderedDict(vd_corrector_kicks_rad)
    hd_corrector_currents_A = OrderedDict(hd_corrector_currents_A)
    vd_corrector_currents_A = OrderedDict(vd_corrector_currents_A)

    if prefer == "currents":
        hd_corrector_kicks_rad = currents_to_kicks_rad(
            hd_corrector_currents_A,
            beam_state=beam_state,
        )
        vd_corrector_kicks_rad = currents_to_kicks_rad(
            vd_corrector_currents_A,
            beam_state=beam_state,
        )
    else:
        hd_corrector_currents_A = kicks_rad_to_currents_A(
            hd_corrector_kicks_rad,
            beam_state=beam_state,
        )
        vd_corrector_currents_A = kicks_rad_to_currents_A(
            vd_corrector_kicks_rad,
            beam_state=beam_state,
        )

    return {
        "hd_corrector_kicks_rad": hd_corrector_kicks_rad,
        "vd_corrector_kicks_rad": vd_corrector_kicks_rad,
        "hd_corrector_currents_A": hd_corrector_currents_A,
        "vd_corrector_currents_A": vd_corrector_currents_A,
    }
