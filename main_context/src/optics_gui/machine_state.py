"""
machine_state.py

Machine-state object for the ISIS RCS optics GUI backend.

This layer defines the complete machine state at one cycle time and provides
helpers to update tune, harmonic, scaling, corrector and error-table settings.
"""

from collections import OrderedDict
from dataclasses import asdict, dataclass, field
from datetime import datetime

from .machine_state_defaults import (
    DEFAULT_BASE_QX,
    DEFAULT_BASE_QY,
    DEFAULT_HARMONICS,
    DEFAULT_KQTD,
    DEFAULT_KQTF,
    DEFAULT_QX,
    DEFAULT_QY,
    DEFAULT_TUNE_METHOD,
    MAIN_MAGNET_SCALING_PRESETS,
)
from .tune_control import calculate_kqtf_kqtd_di
from .correctors import (
    build_corrector_state,
    default_hd_corrector_currents_A,
    default_hd_corrector_kicks_rad,
    default_vd_corrector_currents_A,
    default_vd_corrector_kicks_rad,
)


@dataclass
class MachineState:
    """
    Complete machine state at one ISIS RCS cycle time.
    """

    timestamp: str
    cycle_time_ms: float
    beam_state: object

    main_magnet_mode: str = "rcs_bare"
    main_magnet_scaling: OrderedDict = field(default_factory=OrderedDict)

    tune_method: str = DEFAULT_TUNE_METHOD

    requested_qx: float = DEFAULT_QX
    requested_qy: float = DEFAULT_QY

    base_qx: float = DEFAULT_BASE_QX
    base_qy: float = DEFAULT_BASE_QY

    kqtd: float = DEFAULT_KQTD
    kqtf: float = DEFAULT_KQTF

    iqtd_A: float = 0.0
    iqtf_A: float = 0.0

    harmonic_tunes: OrderedDict = field(default_factory=lambda: OrderedDict(DEFAULT_HARMONICS))

    error_table_path: str = None

    hd_corrector_kicks_rad: OrderedDict = field(default_factory=default_hd_corrector_kicks_rad)
    vd_corrector_kicks_rad: OrderedDict = field(default_factory=default_vd_corrector_kicks_rad)

    hd_corrector_currents_A: OrderedDict = field(default_factory=default_hd_corrector_currents_A)
    vd_corrector_currents_A: OrderedDict = field(default_factory=default_vd_corrector_currents_A)

    metadata: dict = field(default_factory=dict)

    @classmethod
    def from_defaults(
        cls,
        beam_state,
        main_magnet_mode="rcs_bare",
        requested_qx=DEFAULT_QX,
        requested_qy=DEFAULT_QY,
        base_qx=DEFAULT_BASE_QX,
        base_qy=DEFAULT_BASE_QY,
        tune_method=DEFAULT_TUNE_METHOD,
        harmonic_tunes=None,
        error_table_path=None,
        hd_corrector_kicks_rad=None,
        vd_corrector_kicks_rad=None,
        hd_corrector_currents_A=None,
        vd_corrector_currents_A=None,
        corrector_prefer="kicks",
        metadata=None,
        calculate_tune=True,
    ):
        """
        Construct a MachineState from standard defaults plus user overrides.
        """

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if main_magnet_mode not in MAIN_MAGNET_SCALING_PRESETS:
            valid_modes = list(MAIN_MAGNET_SCALING_PRESETS.keys())
            raise ValueError(
                f"Unknown main_magnet_mode '{main_magnet_mode}'. "
                f"Valid modes are {valid_modes}."
            )

        main_magnet_scaling = OrderedDict(MAIN_MAGNET_SCALING_PRESETS[main_magnet_mode])

        if harmonic_tunes is None:
            harmonic_tunes = OrderedDict(DEFAULT_HARMONICS)
        else:
            complete_harmonics = OrderedDict(DEFAULT_HARMONICS)
            complete_harmonics.update(harmonic_tunes)
            harmonic_tunes = complete_harmonics

        corrector_state = build_corrector_state(
            beam_state=beam_state,
            hd_corrector_kicks_rad=hd_corrector_kicks_rad,
            vd_corrector_kicks_rad=vd_corrector_kicks_rad,
            hd_corrector_currents_A=hd_corrector_currents_A,
            vd_corrector_currents_A=vd_corrector_currents_A,
            prefer=corrector_prefer,
        )

        state = cls(
            timestamp=timestamp,
            cycle_time_ms=float(beam_state.cycle_time_ms),
            beam_state=beam_state,
            main_magnet_mode=main_magnet_mode,
            main_magnet_scaling=main_magnet_scaling,
            tune_method=tune_method,
            requested_qx=float(requested_qx),
            requested_qy=float(requested_qy),
            base_qx=float(base_qx),
            base_qy=float(base_qy),
            harmonic_tunes=harmonic_tunes,
            error_table_path=error_table_path,
            hd_corrector_kicks_rad=corrector_state["hd_corrector_kicks_rad"],
            vd_corrector_kicks_rad=corrector_state["vd_corrector_kicks_rad"],
            hd_corrector_currents_A=corrector_state["hd_corrector_currents_A"],
            vd_corrector_currents_A=corrector_state["vd_corrector_currents_A"],
            metadata={} if metadata is None else dict(metadata),
        )

        if calculate_tune:
            state.calculate_trim_quad_strengths()

        return state

    def calculate_trim_quad_strengths(self):
        """
        Calculate kqtd/kqtf from the requested tunes.
        """

        if self.tune_method == "di_wright":
            result = calculate_kqtf_kqtd_di(
                qx=self.requested_qx,
                qy=self.requested_qy,
                beam_state=self.beam_state,
                base_qx=self.base_qx,
                base_qy=self.base_qy,
            )

            self.kqtf = result["kqtf"]
            self.kqtd = result["kqtd"]
            self.iqtf_A = result["iqtf_A"]
            self.iqtd_A = result["iqtd_A"]

        elif self.tune_method == "madx_match":
            # MAD-X tune matching is handled by MadxModel.
            pass

        elif self.tune_method in ("manual", "none"):
            pass

        else:
            raise ValueError(f"Unknown tune_method: {self.tune_method}")

        return self.kqtf, self.kqtd

    def update_harmonics(self, **kwargs):
        """
        Update harmonic tune variables.
        """

        for key, value in kwargs.items():
            if key not in self.harmonic_tunes:
                raise KeyError(f"Unknown harmonic tune variable: {key}")

            self.harmonic_tunes[key] = float(value)

        return self.harmonic_tunes

    def update_correctors(
        self,
        hd_corrector_kicks_rad=None,
        vd_corrector_kicks_rad=None,
        hd_corrector_currents_A=None,
        vd_corrector_currents_A=None,
        prefer="kicks",
    ):
        """
        Update corrector kicks/currents and keep both representations available.
        """

        corrector_state = build_corrector_state(
            beam_state=self.beam_state,
            hd_corrector_kicks_rad=(
                self.hd_corrector_kicks_rad
                if hd_corrector_kicks_rad is None
                else hd_corrector_kicks_rad
            ),
            vd_corrector_kicks_rad=(
                self.vd_corrector_kicks_rad
                if vd_corrector_kicks_rad is None
                else vd_corrector_kicks_rad
            ),
            hd_corrector_currents_A=(
                self.hd_corrector_currents_A
                if hd_corrector_currents_A is None
                else hd_corrector_currents_A
            ),
            vd_corrector_currents_A=(
                self.vd_corrector_currents_A
                if vd_corrector_currents_A is None
                else vd_corrector_currents_A
            ),
            prefer=prefer,
        )

        self.hd_corrector_kicks_rad = corrector_state["hd_corrector_kicks_rad"]
        self.vd_corrector_kicks_rad = corrector_state["vd_corrector_kicks_rad"]
        self.hd_corrector_currents_A = corrector_state["hd_corrector_currents_A"]
        self.vd_corrector_currents_A = corrector_state["vd_corrector_currents_A"]

        return corrector_state

    def update_main_magnet_mode(self, main_magnet_mode):
        """
        Switch between predefined main-magnet scaling modes.
        """

        if main_magnet_mode not in MAIN_MAGNET_SCALING_PRESETS:
            valid_modes = list(MAIN_MAGNET_SCALING_PRESETS.keys())
            raise ValueError(
                f"Unknown main_magnet_mode '{main_magnet_mode}'. "
                f"Valid modes are {valid_modes}."
            )

        self.main_magnet_mode = main_magnet_mode
        self.main_magnet_scaling = OrderedDict(
            MAIN_MAGNET_SCALING_PRESETS[main_magnet_mode]
        )

        return self.main_magnet_scaling

    def update_error_table(self, error_table_path):
        """
        Store the MAD-X error-table path associated with this machine state.
        """

        self.error_table_path = error_table_path
        return self.error_table_path

    def beam_summary_dict(self):
        """
        Return beam_state as a dictionary.
        """

        if hasattr(self.beam_state, "summary_dict"):
            return self.beam_state.summary_dict()

        try:
            return asdict(self.beam_state)
        except TypeError:
            return dict(self.beam_state.__dict__)

    def summary_dict(self):
        """
        Return a serialisable summary of the complete machine state.
        """

        return {
            "timestamp": self.timestamp,
            "cycle_time_ms": self.cycle_time_ms,
            "beam_state": self.beam_summary_dict(),
            "main_magnet_mode": self.main_magnet_mode,
            "main_magnet_scaling": dict(self.main_magnet_scaling),
            "tune_method": self.tune_method,
            "requested_qx": self.requested_qx,
            "requested_qy": self.requested_qy,
            "base_qx": self.base_qx,
            "base_qy": self.base_qy,
            "kqtd": self.kqtd,
            "kqtf": self.kqtf,
            "iqtd_A": self.iqtd_A,
            "iqtf_A": self.iqtf_A,
            "harmonic_tunes": dict(self.harmonic_tunes),
            "error_table_path": self.error_table_path,
            "hd_corrector_kicks_rad": dict(self.hd_corrector_kicks_rad),
            "vd_corrector_kicks_rad": dict(self.vd_corrector_kicks_rad),
            "hd_corrector_currents_A": dict(self.hd_corrector_currents_A),
            "vd_corrector_currents_A": dict(self.vd_corrector_currents_A),
            "metadata": dict(self.metadata),
        }
