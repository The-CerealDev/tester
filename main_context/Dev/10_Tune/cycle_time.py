"""
cycle_time.py

ISIS RCS cycle-time and beam-parameter utilities.

This module provides:
    - standalone relativistic conversion functions
    - an immutable RCSState dataclass
    - an RCSRamp class for converting cycle time into beam state

Default ramp model:
    Idealised sinusoidal main-magnet ramp from injection energy to top energy
    over a 10 ms acceleration cycle at 50 Hz.
"""

from dataclasses import dataclass
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

PROTON_REST_MASS_MEV = 938.27208816
SPEED_OF_LIGHT = 299792458.0
GEV_PER_MEV = 1.0e-3

DEFAULT_INJECTION_ENERGY_MEV = 70.0
DEFAULT_TOP_ENERGY_MEV = 800.0
DEFAULT_CYCLE_START_MS = 0.0
DEFAULT_CYCLE_END_MS = 10.0
DEFAULT_RAMP_FREQUENCY_HZ = 50.0


# ----------------------------------------------------------------------
# Relativistic helper functions
# ----------------------------------------------------------------------

def lorentz_gamma(total_energy_MeV, rest_mass_MeV=PROTON_REST_MASS_MEV):
    return total_energy_MeV / rest_mass_MeV


def lorentz_gamma_from_beta(beta):
    return 1.0 / np.sqrt(1.0 - beta**2)


def lorentz_beta(gamma):
    return np.sqrt(1.0 - 1.0 / gamma**2)


def total_energy_from_gamma_MeV(gamma, rest_mass_MeV=PROTON_REST_MASS_MEV):
    return gamma * rest_mass_MeV


def total_energy_from_kinetic_MeV(kinetic_energy_MeV, rest_mass_MeV=PROTON_REST_MASS_MEV):
    return kinetic_energy_MeV + rest_mass_MeV


def kinetic_energy_from_total_MeV(total_energy_MeV, rest_mass_MeV=PROTON_REST_MASS_MEV):
    return total_energy_MeV - rest_mass_MeV


def momentum_from_total_energy_MeV_c(total_energy_MeV, rest_mass_MeV=PROTON_REST_MASS_MEV):
    return np.sqrt(total_energy_MeV**2 - rest_mass_MeV**2)


def momentum_from_kinetic_energy_MeV_c(kinetic_energy_MeV, rest_mass_MeV=PROTON_REST_MASS_MEV):
    total_energy_MeV = total_energy_from_kinetic_MeV(
        kinetic_energy_MeV,
        rest_mass_MeV=rest_mass_MeV,
    )

    return momentum_from_total_energy_MeV_c(
        total_energy_MeV,
        rest_mass_MeV=rest_mass_MeV,
    )


def relativistic_momentum_MeV_c(gamma, rest_mass_MeV=PROTON_REST_MASS_MEV):
    return gamma * rest_mass_MeV * lorentz_beta(gamma)


def brho_from_momentum_MeV_c(momentum_MeV_c):
    """
    Convert momentum to magnetic rigidity.

    B rho [T m] = p [GeV/c] / 0.299792458
    """

    momentum_GeV_c = momentum_MeV_c * GEV_PER_MEV

    return momentum_GeV_c / 0.299792458


def momentum_from_brho_MeV_c(brho_Tm):
    """
    Convert magnetic rigidity to momentum.

    p [GeV/c] = B rho [T m] * 0.299792458
    """

    momentum_GeV_c = brho_Tm * 0.299792458

    return momentum_GeV_c / GEV_PER_MEV


# ----------------------------------------------------------------------
# State container
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class RCSState:
    """
    Beam state at a single cycle time.
    """

    cycle_time_ms: float
    kinetic_energy_MeV: float
    total_energy_MeV: float
    gamma: float
    beta: float
    momentum_MeV_c: float
    brho_Tm: float
    normalised_momentum: float

    def kinetic_energy_GeV(self):
        return self.kinetic_energy_MeV * GEV_PER_MEV

    def total_energy_GeV(self):
        return self.total_energy_MeV * GEV_PER_MEV

    def momentum_GeV_c(self):
        return self.momentum_MeV_c * GEV_PER_MEV

    def summary_dict(self):
        return {
            "cycle_time_ms": self.cycle_time_ms,
            "kinetic_energy_MeV": self.kinetic_energy_MeV,
            "kinetic_energy_GeV": self.kinetic_energy_GeV(),
            "total_energy_MeV": self.total_energy_MeV,
            "total_energy_GeV": self.total_energy_GeV(),
            "gamma": self.gamma,
            "beta": self.beta,
            "momentum_MeV_c": self.momentum_MeV_c,
            "momentum_GeV_c": self.momentum_GeV_c(),
            "brho_Tm": self.brho_Tm,
            "normalised_momentum": self.normalised_momentum,
        }

    def print_beam(self):
        print(f"M_proton = {PROTON_REST_MASS_MEV:.8g} MeV")
        print(f"Cycle time = {self.cycle_time_ms:.8g} ms")
        print(f"Energy = {self.kinetic_energy_GeV():.8g} GeV")
        print(f"Energy = {self.kinetic_energy_MeV:.8g} MeV")
        print(f"Total energy = {self.total_energy_GeV():.8g} GeV")
        print(f"Gamma = {self.gamma:.8g}")
        print(f"Beta = {self.beta:.8g}")
        print(f"Momentum = {self.momentum_GeV_c():.8g} GeV/c")
        print(f"Rigidity = {self.brho_Tm:.8g} Tm")
        print(f"Normalised momentum = {self.normalised_momentum:.8g}")

    def __str__(self):
        return (
            f"RCSState("
            f"t={self.cycle_time_ms:.4g} ms, "
            f"E_kin={self.kinetic_energy_MeV:.6g} MeV, "
            f"p={self.momentum_GeV_c():.6g} GeV/c, "
            f"Brho={self.brho_Tm:.6g} Tm, "
            f"pn={self.normalised_momentum:.6g}"
            f")"
        )


# ----------------------------------------------------------------------
# Ramp model
# ----------------------------------------------------------------------
class RCSRamp:
    """
    ISIS RCS beam-energy and momentum ramp.

    Model:
        Sinusoidal AC main-magnet ramp over the acceleration cycle.

    Momentum follows:
        p_n(t) = 0.5 * (p_r + 1 - (p_r - 1) * cos(pi * f))

    where:
        f   = (t - t_start) / (t_end - t_start)
        p_r = p_top / p_injection

    Parameters
    ----------
    top_energy_MeV : float
        Final kinetic energy at the end of the acceleration ramp.

    injection_energy_MeV : float, optional
        Injection kinetic energy.

    cycle_start_ms : float, optional
        Start of the modelled ramp.

    cycle_end_ms : float, optional
        End of the modelled ramp.

    ramp_frequency_Hz : float, optional
        Main magnet ramp frequency. Default is 50 Hz.

    rest_mass_MeV : float, optional
        Proton rest mass energy.
    """

    def __init__(
        self,
        top_energy_MeV=DEFAULT_TOP_ENERGY_MEV,
        injection_energy_MeV=DEFAULT_INJECTION_ENERGY_MEV,
        cycle_start_ms=DEFAULT_CYCLE_START_MS,
        cycle_end_ms=DEFAULT_CYCLE_END_MS,
        ramp_frequency_Hz=DEFAULT_RAMP_FREQUENCY_HZ,
        rest_mass_MeV=PROTON_REST_MASS_MEV,
    ):
        self.top_energy_MeV = float(top_energy_MeV)
        self.injection_energy_MeV = float(injection_energy_MeV)
        self.cycle_start_ms = float(cycle_start_ms)
        self.cycle_end_ms = float(cycle_end_ms)
        self.ramp_frequency_Hz = float(ramp_frequency_Hz)
        self.rest_mass_MeV = float(rest_mass_MeV)

        if self.cycle_end_ms <= self.cycle_start_ms:
            raise ValueError("cycle_end_ms must be greater than cycle_start_ms.")

        if self.top_energy_MeV <= self.injection_energy_MeV:
            raise ValueError("top_energy_MeV must be greater than injection_energy_MeV.")

        self.injection_momentum_MeV_c = momentum_from_kinetic_energy_MeV_c(
            self.injection_energy_MeV,
            rest_mass_MeV=self.rest_mass_MeV,
        )

        self.top_momentum_MeV_c = momentum_from_kinetic_energy_MeV_c(
            self.top_energy_MeV,
            rest_mass_MeV=self.rest_mass_MeV,
        )

        self.momentum_ratio = self.top_momentum_MeV_c / self.injection_momentum_MeV_c

    def _as_array(self, value):
        return np.asarray(value, dtype=float)

    def _clip_cycle_time(self, cycle_time_ms):
        return np.clip(
            self._as_array(cycle_time_ms),
            self.cycle_start_ms,
            self.cycle_end_ms,
        )

    def ramp_fraction_at(self, cycle_time_ms):
        """
        Return ramp fraction between 0 and 1.
        """

        t = self._clip_cycle_time(cycle_time_ms)

        return (t - self.cycle_start_ms) / (
            self.cycle_end_ms - self.cycle_start_ms
        )

    def normalised_momentum_at(self, cycle_time_ms):
        """
        Return momentum normalised to injection momentum.
        """

        f = self.ramp_fraction_at(cycle_time_ms)
        pr = self.momentum_ratio

        return 0.5 * (pr + 1.0 - (pr - 1.0) * np.cos(np.pi * f))

    def momentum_at(self, cycle_time_ms):
        """
        Return beam momentum in MeV/c.
        """

        return self.injection_momentum_MeV_c * self.normalised_momentum_at(
            cycle_time_ms
        )

    def total_energy_at(self, cycle_time_ms):
        """
        Return total beam energy in MeV.
        """

        momentum = self.momentum_at(cycle_time_ms)

        return np.sqrt(momentum**2 + self.rest_mass_MeV**2)

    def energy_at(self, cycle_time_ms):
        """
        Return kinetic beam energy in MeV.
        """

        return kinetic_energy_from_total_MeV(
            self.total_energy_at(cycle_time_ms),
            rest_mass_MeV=self.rest_mass_MeV,
        )

    def kinetic_energy_at(self, cycle_time_ms):
        """
        Alias for energy_at().
        """

        return self.energy_at(cycle_time_ms)

    def gamma_at(self, cycle_time_ms):
        """
        Return relativistic gamma.
        """

        return lorentz_gamma(
            self.total_energy_at(cycle_time_ms),
            rest_mass_MeV=self.rest_mass_MeV,
        )

    def beta_at(self, cycle_time_ms):
        """
        Return relativistic beta = v / c.
        """

        return lorentz_beta(self.gamma_at(cycle_time_ms))

    def brho_at(self, cycle_time_ms):
        """
        Return magnetic rigidity in T m.
        """

        return brho_from_momentum_MeV_c(self.momentum_at(cycle_time_ms))

    def state_at(self, cycle_time_ms):
        """
        Return RCSState at a single cycle time.

        This method is intended for scalar cycle times.
        """

        cycle_time_ms = float(cycle_time_ms)

        kinetic_energy_MeV = float(self.energy_at(cycle_time_ms))
        total_energy_MeV = float(self.total_energy_at(cycle_time_ms))
        gamma = float(self.gamma_at(cycle_time_ms))
        beta = float(self.beta_at(cycle_time_ms))
        momentum_MeV_c = float(self.momentum_at(cycle_time_ms))
        brho_Tm = float(self.brho_at(cycle_time_ms))
        normalised_momentum = float(self.normalised_momentum_at(cycle_time_ms))

        return RCSState(
            cycle_time_ms=cycle_time_ms,
            kinetic_energy_MeV=kinetic_energy_MeV,
            total_energy_MeV=total_energy_MeV,
            gamma=gamma,
            beta=beta,
            momentum_MeV_c=momentum_MeV_c,
            brho_Tm=brho_Tm,
            normalised_momentum=normalised_momentum,
        )

    def states_at(self, cycle_time_array_ms):
        """
        Return a list of RCSState objects for multiple cycle times.
        """

        return [self.state_at(t) for t in cycle_time_array_ms]

    def dataframe_at(self, cycle_time_array_ms):
        """
        Return beam parameters as a pandas DataFrame.
        """

        states = self.states_at(cycle_time_array_ms)

        return pd.DataFrame([state.summary_dict() for state in states])

    def full_cycle_dataframe(self, intervals=20):
        """
        Return beam parameters across the full acceleration cycle.
        """

        time_array = np.linspace(
            self.cycle_start_ms,
            self.cycle_end_ms,
            int(intervals) + 1,
        )

        return self.dataframe_at(time_array)
