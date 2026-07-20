"""
tune_control.py

Tune-control helpers for the ISIS RCS machine-state layer.
"""

import numpy as np

from machine_state_defaults import (
    DEFAULT_BASE_QX,
    DEFAULT_BASE_QY,
    DEFAULT_DI_TUNE_COEFFICIENTS,
    DEFAULT_TQ_GCAL,
)


def tune_to_trim_quad_current_di(
    qx,
    qy,
    base_qx=DEFAULT_BASE_QX,
    base_qy=DEFAULT_BASE_QY,
    pn=1.0,
    coefficients=None,
):
    """
    Calculate QTF/QTD currents using Di Wright's tune-control equations.

    Parameters
    ----------
    qx, qy : float
        Requested horizontal and vertical tunes.

    base_qx, base_qy : float
        Base tunes.

    pn : float
        Normalised momentum.

    coefficients : sequence of 4 floats
        Di Wright tune-control coefficients.

    Returns
    -------
    iqtf_A, iqtd_A : float
        QTF and QTD currents in amperes.
    """

    if coefficients is None:
        coefficients = DEFAULT_DI_TUNE_COEFFICIENTS

    z1, z2, z3, z4 = np.asarray(coefficients, dtype=float)

    dqx = float(qx) - float(base_qx)
    dqy = float(qy) - float(base_qy)

    denominator = z1 * z4 - z2 * z3

    if denominator == 0.0:
        raise ZeroDivisionError("Tune-control coefficient denominator is zero.")

    iqtf = pn * (z1 * dqy - z3 * dqx) / denominator
    iqtd = pn * (z4 * dqx - z2 * dqy) / denominator

    return -float(iqtf), -float(iqtd)


def trim_quad_current_to_tune_di(
    iqtf_A,
    iqtd_A,
    base_qx=DEFAULT_BASE_QX,
    base_qy=DEFAULT_BASE_QY,
    pn=1.0,
    coefficients=None,
):
    """
    Convert QTF/QTD currents back to the tune implied by Di Wright's equations.
    """

    if coefficients is None:
        coefficients = DEFAULT_DI_TUNE_COEFFICIENTS

    z1, z2, z3, z4 = np.asarray(coefficients, dtype=float)

    if pn == 0.0:
        raise ZeroDivisionError("Normalised momentum pn must be non-zero.")

    dqx = (iqtf_A * z2 + z1 * iqtd_A) / pn
    dqy = (iqtf_A * z4 + z3 * iqtd_A) / pn

    qx = float(base_qx) + dqx
    qy = float(base_qy) + dqy

    return float(qx), float(qy)


def current_to_strength(current_A, gcal=DEFAULT_TQ_GCAL, brho_Tm=1.23, pn=1.0):
    """
    Convert trim-quadrupole current to MAD-X strength.

    The convention follows the existing ISIS helper functions:

        k = I * Gcal / Brho / pn
    """

    if brho_Tm == 0.0:
        raise ZeroDivisionError("brho_Tm must be non-zero.")

    if pn == 0.0:
        raise ZeroDivisionError("pn must be non-zero.")

    return float(current_A) * float(gcal) / float(brho_Tm) / float(pn)


def strength_to_current(strength, gcal=DEFAULT_TQ_GCAL, brho_Tm=1.23, pn=1.0):
    """
    Convert MAD-X trim-quadrupole strength to current.
    """

    if gcal == 0.0:
        raise ZeroDivisionError("gcal must be non-zero.")

    return float(strength) * float(brho_Tm) * float(pn) / float(gcal)


def calculate_kqtf_kqtd_di(
    qx,
    qy,
    beam_state,
    base_qx=DEFAULT_BASE_QX,
    base_qy=DEFAULT_BASE_QY,
    coefficients=None,
    gcal=DEFAULT_TQ_GCAL,
):
    """
    Calculate kqtf and kqtd from requested tunes using Di Wright's equations.

    Parameters
    ----------
    qx, qy : float
        Requested tunes.

    beam_state : RCSState-like object
        Must provide brho_Tm and normalised_momentum.

    Returns
    -------
    result : dict
        Dictionary containing currents and strengths:
            iqtf_A, iqtd_A, kqtf, kqtd
    """

    pn = float(beam_state.normalised_momentum)
    brho_Tm = float(beam_state.brho_Tm)

    iqtf_A, iqtd_A = tune_to_trim_quad_current_di(
        qx=qx,
        qy=qy,
        base_qx=base_qx,
        base_qy=base_qy,
        pn=pn,
        coefficients=coefficients,
    )

    kqtf = current_to_strength(
        current_A=iqtf_A,
        gcal=gcal,
        brho_Tm=brho_Tm,
        pn=pn,
    )

    kqtd = current_to_strength(
        current_A=iqtd_A,
        gcal=gcal,
        brho_Tm=brho_Tm,
        pn=pn,
    )

    return {
        "iqtf_A": float(iqtf_A),
        "iqtd_A": float(iqtd_A),
        "kqtf": float(kqtf),
        "kqtd": float(kqtd),
    }
