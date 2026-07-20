########################################
# STFC ISIS Synchrotron Tune Controls
# Python 3.x
# 07.06.23
# Haroon Rafique
# ISIS Accelerator Physics Group
########################################

import numpy as np

# Start with original PGH functions

def tune_to_trim_quad_current_pgh(Qh=4.331, Qv=3.731,
                                  h=np.array([ 0.79319992,  1.99412645,  1.57813958, -1.83046275]),
                                  v=np.array([ 0.99954791,  1.18077384,  2.28017058, -1.4796457 ]),
                                  superperiods=10, h_period=0, h_sign=1, v_period=0, v_sign=1,
                                  Gcal=np.array([1.997e-3, 1.997e-3]), Brho=1.23, pn=1.0):
    '''
    Calculates the trim quad currents required to obtain the given tune under
    PGH's proposed tune control system.

    Parameters
    ----------
    Qh : Float, optional
        Horizontal tune required.
        Can be a float array of tunes. Must be of same length as Qv.
        The default is 4.331.
    Qv : Float, optional
        Vertical tune required.
        Can be a float array of tunes. Must be of same length as Qh
        The default is 3.731.
    h : Float array of length 4, optional
        Horizontal one turn transfer matrix trace coefficients.
        The default is np.array([0.79716154, 2.007355  , 1.58717351, -1.82842749]).
    v : Float array of length 4, optional
        Vertical one turn transfer matrix trace coefficients.
        The default is np.array([0.99538634, 1.17390209, 2.26760524, -1.4836992 ]).
    superperiods : Int, optional
        Number of superperiods in the machine.
        The default is 10.
    h_period : Int, optional
        arccos period of horizontal phase advance.
        The default is 0.
    h_sign : Int, optional
        arccos sign of horizontal phase advance, should be +/- 1.
        The default is 1.
    v_period : Int, optional
        arccos period of vertical phase advance.
        The default is 0.
    v_sign : Int, optional
        arccos sign of vertical phase advance, should be +/- 1.
        The default is 1.
    Gcal : Float array of length 2, optional
        QTF, QTD calibrations.
        The default is np.array([1.997e-3, 1.997e-3]).
    Brho : Float, optional
        Magnetic rigidity.
        The default is 1.23.
    pn : Float, optional
        Normalised momentum.
        Can be a float array of momenta. Must be of same length as Qh, Qv.
        The default is 1.0.

    Returns
    -------
    I_F : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTF to obtain Qh, Qv according to PGH controls.
    I_D : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTD to obtain Qh, Qv according to PGH controls.
    '''

    # Calculate the required trace of the one turn transfer matrix
    phih = tune_to_trace(Qy=Qh, superperiods=superperiods, arccos_sign=h_sign, arccos_period=h_period)
    phiv = tune_to_trace(Qy=Qv, superperiods=superperiods, arccos_sign=v_sign, arccos_period=v_period)
    
    # Calculate the coefficients of the quadratic equation that gives k_QTF
    a, b, c = trace_to_coeffs(phih=phih, phiv=phiv, h=h, v=v)

    # Solve the simultaneous equations
    k_F, k_D = coeffs_to_strengths(a=a, b=b, c=c, phih=phih, h=h)
    
    # Convert the strengths into currents
    Gcal_F, Gcal_D = Gcal
    I_F = strength_to_current(k=k_F, Gcal=Gcal_F, Brho=Brho, pn=pn)
    I_D = strength_to_current(k=k_D, Gcal=Gcal_D, Brho=Brho, pn=pn)

    # Return the currents
    return I_F, I_D


def trim_quad_current_to_tune_pgh(I_F=0.0, I_D=0.0,
                                  h=np.array([ 0.79319992,  1.99412645,  1.57813958, -1.83046275]),
                                  v=np.array([ 0.99954791,  1.18077384,  2.28017058, -1.4796457 ]),
                                  superperiods=10, h_period=0, h_sign=1, v_period=0, v_sign=1,
                                  Gcal=np.array([1.997e-3, 1.997e-3]), Brho=1.23, pn=1.0):
    '''
    For given trim quad currents, calculates the tunes reported by
    PGH's proposed tune control system.

    Parameters
    ----------
    I_F : TYPE, optional
        Current applied to QTF.
        Can be a float array of currents. Must be of same length as I_D.
        The default is 0.0.
    I_D : TYPE, optional
        Current applied to QTD.
        Can be a float array of currents. Must be of same length as I_F.
        The default is 0.0.
    h : Float array of length 4, optional
        Horizontal one turn transfer matrix trace coefficients.
        The default is np.array([0.79716154, 2.007355  , 1.58717351, -1.82842749]).
    v : Float array of length 4, optional
        Vertical one turn transfer matrix trace coefficients.
        The default is np.array([0.99538634, 1.17390209, 2.26760524, -1.4836992 ]).
    superperiods : Int, optional
        Number of superperiods in the machine.
        The default is 10.
    h_period : Int, optional
        arccos period of horizontal phase advance.
        The default is 0.
    h_sign : Int, optional
        arccos sign of horizontal phase advance, should be +/- 1.
        The default is 1.
    v_period : Int, optional
        arccos period of vertical phase advance.
        The default is 0.
    v_sign : Int, optional
        arccos sign of vertical phase advance, should be +/- 1.
        The default is 1.
    Gcal : Float array of length 2, optional
        QTF, QTD calibrations.
        The default is np.array([1.997e-3, 1.997e-3]).
    Brho : Float, optional
        Magnetic rigidity.
        The default is 1.23.
    pn : Float, optional
        Normalised momentum.
        Can be a float array of momenta. Must be of same length as I_F, I_D.
        The default is 1.0.

    Returns
    -------
    Qh : Float, or array of floats with length same as I_F, I_D
        Horizontal tune according to PGH controls.
    Qv : Float, or array of floats with length same as I_F, I_D
        Vertical tune according to to PGH controls.

    '''
    
    # Calculate the trim quad strengths
    Gcal_F, Gcal_D = Gcal
    k_F = current_to_strength(I=I_F, Gcal=Gcal_F, Brho=Brho, pn=pn)
    k_D = current_to_strength(I=I_D, Gcal=Gcal_D, Brho=Brho, pn=pn)

    # Calculate the trace of the one turn transfer matrix
    phih, phiv = strength_to_trace(k_F=k_F, k_D=k_D, h=h, v=v)

    # Calculate the tune
    Qh = trace_to_tune(phi=phih, superperiods=superperiods, arccos_sign=h_sign, arccos_period=h_period)
    Qv = trace_to_tune(phi=phiv, superperiods=superperiods, arccos_sign=v_sign, arccos_period=v_period)

    # Return the tune
    return Qh, Qv


def tune_to_trim_quad_current_di(Qh=4.331, Qv=3.731,
                                 baseQh=4.331, baseQv=3.731, pn=1.0,
                                 z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    '''
    Calculates the trim quad currents required to obtain the given tune under
    Di's tune control system. The default values for the base tunes and
    normalised momentum are those at t = 0.0 ms.

    Parameters
    ----------
    Qh : Float, optional
        Horizontal tune required.
        Can be a float array of tunes. Must be of same length as Qv.
        The default is 4.331.
    Qv : Float, optional
        Vertical tune required.
        Can be a float array of tunes. Must be of same length as Qh
        The default is 3.731.
    baseQh : Float, optional
        Base horizontal tune.
        Can be a float array of tunes. Must be of same length as Qh, Qv.
        The default is 4.331.
    baseQv : Float, optional
        Base vertical tune.
        Can be a float array of tunes. Must be of same length as Qh, Qv.
        The default is 3.731.
    pn : Float, optional
        Normalised momentum.
        Can be a float array of momenta. Must be of same length as Qh, Qv.
        The default is 1.0.
    z : Float array of length 4, optional
        Coefficients of the tune control system.
        The default is np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]).

    Returns
    -------
    I_F : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTF to obtain Qh, Qv according to Di controls.
    I_D : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTD to obtain Qh, Qv according to Di controls.

    '''

    # Get the control coefficients
    z1, z2, z3, z4 = z

    # Calculate the change in tune required
    dQh = Qh - baseQh
    dQv = Qv - baseQv

    # Calculate the currents required
    I_F = pn * (z1 * dQv - z3 * dQh) / (z1 * z4 - z2 * z3)
    I_D = pn * (z4 * dQh - z2 * dQv) / (z1 * z4 - z2 * z3)

    # Return the currents
    return -I_F, -I_D


def trim_quad_current_to_tune_di(I_F=0.0, I_D=0.0,
                                 baseQh=4.331, baseQv=3.731, pn=1.0,
                                 z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    '''
    For given trim quad currents, calculates the tunes reported by
    Di's tune control system. The default values for the base tunes and
    normalised momentum are those at t = 0.0 ms.

    Parameters
    ----------
    I_F : TYPE, optional
        Current applied to QTF.
        Can be a float array of currents. Must be of same length as I_D.
        The default is 0.0.
    I_D : TYPE, optional
        Current applied to QTD.
        Can be a float array of currents. Must be of same length as I_F.
        The default is 0.0.
    baseQh : Float, optional
        Base horizontal tune.
        Can be a float array of tunes. Must be of same length as I_F, I_D.
        The default is 4.331.
    baseQv : Float, optional
        Base vertical tune.
        Can be a float array of tunes. Must be of same length as I_F, I_D.
        The default is 3.731.
    pn : Float, optional
        Normalised momentum.
        Can be a float array of momenta. Must be of same length as I_F, I_D.
        The default is 1.0.
    z : Float array of length 4, optional
        Coefficients of the tune control system.
        The default is np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]).

    Returns
    -------
    Qh : Float, or array of floats with length same as I_F, I_D
        Horizontal tune according to to Di controls.
    Qv : Float, or array of floats with length same as I_F, I_D
        Vertical tune according to to Di controls.

    '''
    
    # Get the control coefficients
    z1, z2, z3, z4 = z

    # Calculate the change in tune
    dQh = (I_F * z2 + z1 * I_D) / pn
    dQv = (I_F * z4 + z3 * I_D) / pn

    # Calculate the tune
    Qh = baseQh + dQh
    Qv = baseQv + dQv
    
    # Return the tune
    return Qh, Qv


def di_to_pgh(Qh_di=4.331, Qv_di=3.371,
              baseQh=4.331, baseQv=3.731, pn=1.0,
              z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]),
              h=np.array([ 0.79319992,  1.99412645,  1.57813958, -1.83046275]),
              v=np.array([ 0.99954791,  1.18077384,  2.28017058, -1.4796457 ]),
              superperiods=10, h_period=0, h_sign=1, v_period=0, v_sign=1,
              Gcal=np.array([1.997e-3, 1.997e-3]), Brho=1.23):
    '''
    This function calculates the new tunes required under PGH controls that 
    return the same currents as the function input tunes would return under
    Di controls.
    This function is intended for finding the new tunes to input when replacing
    Di controls with PGH controls.

    Parameters
    ----------
    Qh_Di : Float, optional
           Horizontal tune in Di controls.
           Can be a float array of tunes. Must be of same length as Qv.
           The default is 4.331.
    Qv_Di : Float, optional
           Vertical tune in Di controls.
           Can be a float array of tunes. Must be of same length as Qh
           The default is 3.731.

    Other Parameters are as per tune_to_trim_quad_current and inverse functions.

    Returns
    -------
    Qh_pgh : Float, or array of floats with length same as Qh
            Horizontal tune for PGH controls to give same currents
            that Qh_di would under Di controls.
    Qv_pgh : Float, or array of floats with length same as Qv
            Vertical tune for PGH controls to give same currents
            that Qv_di would under Di controls.
    '''

    # Get currents Di's controls output
    I_F, I_D = tune_to_trim_quad_current_di(Qh=Qh_di, Qv=Qv_di,
                                            baseQh=baseQh, baseQv=baseQv,
                                            pn=pn, z=z)

    # Get tunes that PGH's controls would give the same currents
    Qh_pgh, Qv_pgh = trim_quad_current_to_tune_pgh(I_F=I_F, I_D=I_D, h=h, v=v,
                                                   superperiods=superperiods,
                                                   h_period=h_period, h_sign=h_sign,
                                                   v_period=v_period, v_sign=v_sign,
                                                   Gcal=Gcal, Brho=Brho, pn=pn)

    # Return the tunes
    return Qh_pgh, Qv_pgh


def tune_to_trace(Qy, superperiods=10, arccos_sign=1, arccos_period=0):
    phiy = 2 * np.cos(-arccos_sign * 2 * np.pi * (Qy / superperiods - arccos_period))
    return phiy


def trace_to_coeffs(phih, phiv,
                    h=np.array([ 0.79319992,  1.99412645,  1.57813958, -1.83046275]),
                    v=np.array([ 0.99954791,  1.18077384,  2.28017058, -1.4796457 ])):
    # Get one turn transfer matrix trace coefficients
    h1, h2, h3, h4 = h
    v1, v2, v3, v4 = v

    # Calculate the coefficients of the quadratic equation that gives k_QTF
    a = (-h1 * v3 - h3 * v1)
    b = (phih * v1 - phiv * h1) + ((h1 * v4 - h4 * v1) + (h3 * v2 - h2 * v3))
    c = (-phih * v2 - phiv * h2) + (h2 * v4 + h4 * v2)

    # Return the coefficients
    return a, b, c


def coeffs_to_strengths(a, b, c, phih, h):
    # Get one turn transfer matrix trace coefficients
    h1, h2, h3, h4 = h

    # Solve the simultaneous equations for k_QTF, k_QTD
    k_F1 = (-b + (b**2 - 4 * a * c)**0.5) / (2 * a)
    k_F2 = (-b - (b**2 - 4 * a * c)**0.5) / (2 * a)
    k_D1 = (phih - h3 * k_F1 - h4) / (h1 * k_F1 + h2)
    k_D2 = (phih - h3 * k_F2 - h4) / (h1 * k_F2 + h2)

    # Find the solutions with the minimum trim quad strengths
    k_F, k_D = find_min_k(k_F1, k_F2, k_D1, k_D2)

    # Return trim quad strengths
    return k_F, k_D


def find_min_k(k_F1, k_F2, k_D1, k_D2):
    k_F, k_D = np.where((abs(k_F1) + abs(k_D1)) < (abs(k_F2) + abs(k_D2)),
                        np.array([k_F1, k_D1]), np.array([k_F2, k_D2]))   
    return k_F, k_D


def strength_to_current(k, Gcal=1.997e-3, Brho=1.23, pn=1.0):
    I = k * Brho / Gcal * pn
    return I


def current_to_strength(I, Gcal=1.997e-3, Brho=1.23, pn=1.0):
    k = I * Gcal / Brho / pn
    return k


def strength_to_trace(k_F, k_D,
                      h=np.array([ 0.79319992,  1.99412645,  1.57813958, -1.83046275]),
                      v=np.array([ 0.99954791,  1.18077384,  2.28017058, -1.4796457 ])):
    # Get one turn transfer matrix trace coefficients
    h1, h2, h3, h4 = h
    v1, v2, v3, v4 = v
    
    # Calculate the trace of the one turn transfer matrix
    phih = h1 * k_F * k_D + h2 * k_D + h3 * k_F + h4
    phiv = v1 * k_F * k_D - v2 * k_D - v3 * k_F + v4

    # Return the trace
    return phih, phiv


def trace_to_tune(phi, superperiods=10, arccos_sign=1, arccos_period=0):
    Qy = superperiods / (2 * np.pi) * (2 * np.pi * arccos_period + arccos_sign * np.arccos(phi / 2))
    return Qy

