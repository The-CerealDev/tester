# ISIS RCS Optics GUI Backend

This repository is developing the backend/model layer for an ISIS RCS optics GUI. The current aim is a read-only modelling and virtual-diagnostic workflow: build a machine state at a selected cycle time, run MAD-X through `cpymad`, compare model branches, match/correct closed orbits, evaluate beam envelopes, calculate aperture margins, analyse tune working points, and later support measurement-comparison workflows.

This is not a controls-writing layer. The MAD-X lattice files under `Dev/Lattice_Files/` are fixed inputs and should not be modified by the GUI/backend state layer. Operational state is applied through generated machine-state override files or direct `MachineState` application.

The old standalone `strength_scaling.py` concept has been superseded by the current machine-state approach: `cycle_time.py` produces the beam state; `MachineState` stores all GUI/backend settings for one cycle time; `tune_control.py` calculates supported trim-quadrupole settings; `correctors.py` converts controls-style currents and MAD-X kicks; `machine_state_writer.py` writes executable MAD-X override/log files; and `MadxModel` applies either a written machine-state file or a `MachineState` object directly.

## Current Architecture

The active implementation currently lives in development folders rather than a packaged root module. `Dev/03_MADX_Model/` also contains local copies of the foundation modules so its notebook can run from that folder; keep those copies consistent with `Dev/02_Machine_State/` until a package layout replaces the development folders. `Dev/04_Error_Table_Utilities/` contains the standalone error-table utility layer and local copies of the foundation/MAD-X modules for its validation notebook. `Dev/05_Orbit_Branches/` contains the orbit branch workflow layer. `Dev/06_Tune_Matching/` contains the tune matching, tune programme and harmonic tune plotting layer. `Dev/07_Orbit_Correction/` contains the read-only COCU-style closed-orbit correction layer. `Dev/08_Envelope/` contains the MAD-X-backed beam-envelope layer. `Dev/09_Aperture/` contains the MAD-X/source-spreadsheet aperture-margin layer. `Dev/10_Tune/` contains the tune/resonance helper layer. `Dev/11_Snapshot/` contains the current whole-system snapshot integration layer. `Dev/04_MADX_Model_Revisit/` contains revisit copies and generated outputs; do not treat it as the active edit target unless explicitly requested. The current documented validation path is:

```text
FOUNDATION LAYERS
├── Dev/01_Cycle_Time/cycle_time.py
├── Dev/02_Machine_State/cycle_time.py
├── Dev/02_Machine_State/machine_state_defaults.py
├── Dev/02_Machine_State/tune_control.py
├── Dev/02_Machine_State/correctors.py
├── Dev/02_Machine_State/machine_state.py
└── Dev/02_Machine_State/machine_state_writer.py

MAD-X EXECUTION LAYER
└── Dev/03_MADX_Model/madx_model.py

ERROR-TABLE UTILITY LAYER
└── Dev/04_Error_Table_Utilities/errors.py

ORBIT BRANCH WORKFLOW LAYER
└── Dev/05_Orbit_Branches/orbit_branch.py

TUNE MATCHING WORKFLOW LAYER
├── Dev/06_Tune_Matching/tune_matching.py
└── Dev/06_Tune_Matching/tune_plots.py

ORBIT CORRECTION WORKFLOW LAYER
└── Dev/07_Orbit_Correction/orbit_correction.py

ENVELOPE WORKFLOW LAYER
└── Dev/08_Envelope/envelope.py

APERTURE WORKFLOW LAYER
└── Dev/09_Aperture/aperture.py

TUNE / RESONANCE HELPER LAYER
└── Dev/10_Tune/tune.py

SNAPSHOT INTEGRATION LAYER
└── Dev/11_Snapshot/snapshot.py

VALIDATION NOTEBOOKS
├── Dev/01_Cycle_Time/01_Test_Cycle_Time.ipynb
├── Dev/02_Machine_State/machine_state_tests.ipynb
├── Dev/03_MADX_Model/madx_model_tests.ipynb
├── Dev/04_Error_Table_Utilities/errors_tests.ipynb
├── Dev/05_Orbit_Branches/orbit_branch_tests.ipynb
├── Dev/06_Tune_Matching/tune_matching_tests.ipynb
├── Dev/07_Orbit_Correction/orbit_correction_tests.ipynb
├── Dev/08_Envelope/envelope_tests.ipynb
├── Dev/09_Aperture/aperture_tests.ipynb
├── Dev/10_Tune/tune_tests.ipynb
└── Dev/11_Snapshot/snapshot_tests.ipynb

INPUT DATA
├── Dev/Lattice_Files/
└── Dev/Error_Tables/

PLANNED NEXT LAYERS
├── io/
├── plotting/
└── utils/
```

## Implemented Layers

### `cycle_time.py`

Bottom-layer beam and ramp model. It converts ISIS RCS cycle time into kinetic and total energy, momentum, magnetic rigidity, relativistic beta/gamma, and normalised momentum. It provides the reusable `RCSState` dataclass and `RCSRamp` class used by higher layers.

Major functionality includes proton/rest-energy constants, injection/top energy, cycle start/end, ramp frequency, relativistic helper functions, `energy_at`, `total_energy_at`, `gamma_at`, `beta_at`, `momentum_at`, `brho_at`, `normalised_momentum_at`, `state_at`, and dataframe/summary helpers where present. Validation is in `Dev/01_Cycle_Time/01_Test_Cycle_Time.ipynb` and cycle-time boundary checks in `Dev/02_Machine_State/machine_state_tests.ipynb`.

### `machine_state_defaults.py`

Central defaults and operational variable names for the machine-state layer. It defines main magnet scaling presets (`RCS_BARE_SCALING`, `SRM_BARE_SCALING`, `MAIN_MAGNET_SCALING_PRESETS`), tune defaults (`DEFAULT_QX`, `DEFAULT_QY`, `DEFAULT_BASE_QX`, `DEFAULT_BASE_QY`, `DEFAULT_TUNE_METHOD`, `DEFAULT_DI_TUNE_COEFFICIENTS`, `DEFAULT_TQ_GCAL`), `DEFAULT_HARMONICS`, operational corrector names (`HD_CORRECTOR_NAMES`, `VD_CORRECTOR_NAMES`), and `CORRECTOR_CALIBRATION`.

### `tune_control.py`

Pure tune-control helper layer. It converts requested tunes into QTF/QTD currents and strengths for the supported Di Wright control law, including forward/reverse consistency checks. Key functions are `tune_to_trim_quad_current_di`, `trim_quad_current_to_tune_di`, `current_to_strength`, `strength_to_current`, and `calculate_kqtf_kqtd_di`.

MAD-X matching is handled by `MadxModel`, not by this pure helper module.

### `correctors.py`

Operational HD/VD corrector helper layer. It stores and converts both controls-facing currents and MAD-X-facing kicks. Key functions include `current_to_kick_rad`, `kick_rad_to_current`, `currents_to_kicks_rad`, `kicks_rad_to_currents_A`, and `build_corrector_state`.

MAD-X uses kicks, but the GUI and future correction-suggestion tools should retain currents as well.

### `machine_state.py`

Complete one-cycle-time machine-state object for the optics GUI/backend. `MachineState` stores timestamp, cycle time, `RCSState`, main magnet mode/scaling, requested and base tunes, tune method, `kqtd`, `kqtf`, `iqtd_A`, `iqtf_A`, harmonic tune settings, error-table path, HD/VD corrector kicks, HD/VD corrector currents, and metadata.

Key methods include `MachineState.from_defaults`, `calculate_trim_quad_strengths`, `update_harmonics`, `update_correctors`, `update_main_magnet_mode`, `update_error_table`, `beam_summary_dict`, and `summary_dict`. `MachineState` stores the request/state and does not run MAD-X. For `tune_method="madx_match"`, the actual MAD-X matching is performed by `MadxModel`.

### `machine_state_writer.py`

Writes timestamped MAD-X-readable machine-state override/log files. Key functions are `make_machine_state_filename`, `write_machine_state_file`, `write_assignment`, `write_comment_block`, and `format_madx_value`.

The generated `.strength` file is both an executable MAD-X override and a human-readable GUI/backend snapshot. It writes beam momentum and Brho, main-magnet scaling variables, trim-quadrupole strengths, harmonic tune variables, and HD/VD corrector kicks as MAD-X assignments. It retains trim-quadrupole currents and corrector currents as comments/log values, and writes a JSON summary sidecar.

### `madx_model.py`

Thin `cpymad`/MAD-X execution layer in `Dev/03_MADX_Model/`. `MadxModel` starts/resets MAD-X, loads base ISIS lattice files, loads an optional aperture file, applies a generated machine-state file, applies a `MachineState` object directly, applies error tables, runs TWISS, returns normalised pandas DataFrames for TWISS and summary, performs MAD-X tune matching from a `MachineState`, and tracks metadata/loaded files.

Important design rules:

- `MadxModel` deliberately does not load `2023.strength`.
- `MadxModel` deliberately does not load `ISIS.injected_beam`.
- Beam state is supplied by the generated machine-state file or direct `MachineState` application.
- Base load order is `ISIS.elements`, `ISIS.sequence`, `ISIS.strength`, optional aperture file, optional generated machine-state override file.
- After applying errors, TWISS should not silently clear those errors by reissuing an inappropriate `USE` sequence call.

### `errors.py`

Standalone MAD-X error-table utility layer in `Dev/04_Error_Table_Utilities/`. It reads, writes, validates, combines, filters, flips, maps and summarises MAD-X/TFS error tables using pandas DataFrames and paths. It also provides a thin `apply_error_table(...)` wrapper that delegates to `MadxModel`; MAD-X execution remains owned by `MadxModel`.

The validation notebook checks pure DataFrame/file behaviour, SETERR-compatible TFS round trips, mapping against TWISS names, missing/invalid failure modes, and direct MAD-X application. It includes single-element orbit-response tests for QD, QF, QDS and DIP errors across the six MAD-X alignment columns (`dx`, `dy`, `ds`, `dphi`, `dtheta`, `dpsi`), with millimetre-level orbit assertions where this lattice physically responds.

### `orbit_branch.py`

GUI-ready orbit branch workflow layer in `Dev/05_Orbit_Branches/`. It defines `OrbitBranchConfig`, `OrbitBranchResult` and `OrbitBranch` for running independent default, misaligned and manually corrected/corrector branches through separate `MadxModel` instances.

Key helpers include `extract_orbit_df`, `compare_orbit_results`, `summarise_orbit_df` and `summarise_orbit_difference`. Outputs are pandas DataFrames plus compact metadata dictionaries so notebooks and the future Streamlit GUI can plot, tabulate and label branch results without running MAD-X directly. This layer does not solve orbit correction; orbit-correction workflows remain a later layer.

### `tune_matching.py` and `tune_plots.py`

GUI-ready tune matching and harmonic tune workflow layer in `Dev/06_Tune_Matching/`. It provides `TuneWorkflowResult`, harmonic input normalisation, explicit harmonic-support checks, MAD-X tune matching wrappers, Di Wright actual-tune evaluation, full-cycle tune-programme evaluation and single-point tune evaluation with TWISS output.

The result objects and tables expose requested/set tunes, actual MAD-X tunes, matched tunes where applicable, QTF/QTD currents, QTF/QTD MAD-X K values, harmonic settings, summaries, metadata and warnings for GUI display. Plot helpers generate resonance working-point plots, harmonic trim-quad expected/programmed comparisons, beta-function variation plots with harmonic expectation waves and trim-current bar plots. Plotting functions consume data/results and do not run MAD-X.

### `orbit_correction.py`

Read-only COCU-style closed-orbit correction workflow layer in `Dev/07_Orbit_Correction/`. It provides BPM measurement normalisation, BPM enable/disable handling, HD/VD corrector selection, MAD-X measured-orbit fitting, MAD-X `CORRECT` steering correction, corrector-current conversion and plotting helpers for correction suggestions.

The operational workflow is `correct_bpm_measurements_with_madx_correct(...)`: fit a viable measured orbit from BPM measurements using non-steering lattice/error kick knobs, run MAD-X `CORRECT` using selected HD/VD steering correctors, read the MAD-X correction table, and convert proposed kicks to controls-facing currents using ISIS corrector calibrations. It is read-only and does not write to the machine.

Dummy or synthetic physics outputs are not acceptable for this layer or future GUI-facing model layers unless explicitly requested as toy examples. If a workflow claims to model optics, orbit, tune, matching, correction, TWISS, or error response, it must run the relevant cpymad/MAD-X path or fail clearly rather than fabricating stand-in results.

### `envelope.py`

MAD-X-backed beam-envelope workflow layer in `Dev/08_Envelope/`. It consumes real TWISS DataFrames from `MadxModel` or established workflow results and combines them with explicit beam assumptions to produce GUI-ready sigma and n-sigma envelope columns.

The layer supports geometric and normalised RMS emittance inputs, with defaults of 300 pi mm mrad in both planes, arbitrary positive sigma scale with 3 sigma as the default, and user-controlled momentum spread with 0.2% as the default. Normalised emittance conversion uses the selected cycle-time beam state's beta/gamma. It also provides plotting helpers for horizontal/vertical envelopes, RMS beam size, and nominal-versus-large-harmonic envelope comparisons; plotting helpers consume result tables and do not run MAD-X. Plus and minus envelope bounds for the same result must use identical line width, style and colour.

### `aperture.py`

MAD-X/source-spreadsheet aperture-margin workflow layer in `Dev/09_Aperture/`. It extracts real MAD-X APERTURE tables through `MadxModel.run_aperture`, loads the copied ISIS 2024 source aperture spreadsheet, aligns aperture data to envelope/TWISS rows, computes horizontal and vertical clearance margins, summarises limiting locations and plots aperture, closed orbit, envelope and x/y margins together. Combined aperture/envelope/margin plots put positive and negative margin lines on the same millimetre axis as aperture and envelope, avoiding a secondary y-axis. The current aperture notebook starts from the repository-local simplified lattice in `Dev/Lattice_Files/00_Simplified_Lattice` and uses a 1 sigma envelope overlay for the ISIS-style aperture comparison while keeping sigma scale as an explicit user/GUI input.

MAD-X can report zero aperture values for drift rows that cannot carry aperture metadata. The aperture layer excludes those rows from margin calculations and plots the source spreadsheet data alongside MAD-X extracted aperture data so the GUI can distinguish model-table limitations from the intended physical aperture curve.

### `tune.py`

Tune/resonance helper layer in `Dev/10_Tune/`. It extracts compact tune summaries from MAD-X summary tables and workflow results, normalises set-vs-predicted tune programme tables, builds working-point tables, generates resonance-line segments, evaluates resonance proximity and prepares tune-diagram input DataFrames.

Set/requested tunes remain distinct from predicted/actual MAD-X tunes. The plotting helpers consume repository DataFrames/results directly, including sliced DataFrames, and do not require callers to unpack arrays just to plot full-cycle or multi-point working points.

### `snapshot.py`

Whole-system integration and GUI contract layer in `Dev/11_Snapshot/`. It builds one or more complete machine snapshots by orchestrating lower-layer functionality: beam state, `MachineState`, MAD-X TWISS/summary, tune summary, orbit tables, optional envelope, optional aperture margins, optional orbit branches, metadata and warnings.

It supports copying and modifying snapshot configs, such as cloning a nominal snapshot and adding an error table; ordered snapshot series for nominal-vs-error or multi-time comparisons; and full-cycle working-point resonance diagram inputs through the existing tune/resonance DataFrame path. It does not duplicate lower-layer physics logic and does not own plotting.

## Validation and Test Status

The notebooks are currently the foundation test gate before implementing the next physics/workflow layers. Before adding a new layer, run the three validation notebooks and ensure any changed behaviour is intentional.

- `Dev/01_Cycle_Time/01_Test_Cycle_Time.ipynb`: validates ramp behaviour, beam-state construction, relativistic quantities, and cycle-time consistency.
- `Dev/02_Machine_State/machine_state_tests.ipynb`: validates defaults, RCS/SRM scaling presets, harmonic tune defaults and updates, Di Wright tune-control calculations, trim-quadrupole current/strength round trips, corrector current/kick round trips, machine-state update methods, summaries, writer output, JSON sidecar output, cycle-time boundary behaviour, and selected negative tests.
- `Dev/03_MADX_Model/madx_model_tests.ipynb`: validates lattice loading through `MadxModel`, written and direct machine-state application, expected loaded files, exclusion of `2023.strength` and `ISIS.injected_beam`, TWISS and summary DataFrames, error-table application and persistence, default/misaligned/corrected branch comparisons, MAD-X tune matching from `MachineState`, aperture plotting support where present, and missing-file failure tests.
- `Dev/04_Error_Table_Utilities/errors_tests.ipynb`: validates error-table normalisation, validation, read/write round trips, combination, filtering, sign flipping, TWISS-name mapping, `MadxModel.apply_error_table(...)` integration, wrapper delegation, negative failure modes, and single-element orbit responses for QD, QF, QDS and DIP errors.
- `Dev/05_Orbit_Branches/orbit_branch_tests.ipynb`: validates default, misaligned, corrector and combined-error branches, orbit extraction, branch comparison metrics, independent output directories and selected failure modes.
- `Dev/06_Tune_Matching/tune_matching_tests.ipynb`: validates harmonic input normalisation, explicit harmonic lattice support, MAD-X tune matching, Di Wright actual-tune evaluation, full-cycle set-vs-actual tune programmes, trim-quadrupole current/K exposure, harmonic expected-vs-programmed TWISS checks, combined harmonic cases, resonance plots, beta variation plots and harmonic trim-quad plots.
- `Dev/07_Orbit_Correction/orbit_correction_tests.ipynb`: validates BPM/corrector selection, MAD-X measured-orbit fitting from three dummy BPM measurements per plane, MAD-X `CORRECT` with selected steering dipoles, current conversion through ISIS calibrations, measured/corrected orbit plots, corrector-current plots and assertions that correction reduces full-monitor RMS.
- `Dev/08_Envelope/envelope_tests.ipynb`: validates MAD-X-backed envelope evaluation from nominal and large-harmonic `MachineState` TWISS tables, default 300 pi mm mrad RMS emittances, 0.2% momentum spread, arbitrary sigma scale, geometric and normalised RMS modes, envelope plots, harmonic comparison plots and selected failure modes.
- `Dev/09_Aperture/aperture_tests.ipynb`: validates real MAD-X APERTURE extraction, ISIS source aperture spreadsheet loading, zero/invalid MAD-X aperture filtering for drift rows, envelope/aperture alignment, clearance margins, limiting-location summaries, aperture/envelope/margin plots, harmonic and closed-orbit-distortion cases, and selected failure modes.
- `Dev/10_Tune/tune_tests.ipynb`: validates tune summary extraction, copied tune programme values, set-vs-predicted tune programme tables, resonance line/proximity generation, DataFrame-backed tune diagram plotting, harmonic trim-quad plots and selected failure modes.
- `Dev/11_Snapshot/snapshot_tests.ipynb`: validates whole-system integration through snapshots: single full snapshots, copied/modified configs, nominal-vs-error series, branch comparisons, full-cycle working-point resonance inputs, existing plotting helpers consuming snapshot/series outputs directly, and no external repository data dependencies.

## Current Development Status

| Layer / file | Status |
| --- | --- |
| `Dev/01_Cycle_Time/cycle_time.py` | Implemented and notebook-tested |
| `Dev/02_Machine_State/machine_state_defaults.py` | Implemented and notebook-tested |
| `Dev/02_Machine_State/tune_control.py` | Implemented and notebook-tested |
| `Dev/02_Machine_State/correctors.py` | Implemented and notebook-tested |
| `Dev/02_Machine_State/machine_state.py` | Implemented and notebook-tested |
| `Dev/02_Machine_State/machine_state_writer.py` | Implemented and notebook-tested |
| `Dev/03_MADX_Model/madx_model.py` | Implemented and notebook-tested |
| `Dev/04_Error_Table_Utilities/errors.py` | Implemented and notebook-tested |
| `Dev/05_Orbit_Branches/orbit_branch.py` | Implemented and notebook-tested |
| `Dev/06_Tune_Matching/tune_matching.py` | Implemented and notebook-tested |
| `Dev/06_Tune_Matching/tune_plots.py` | Implemented and notebook-tested |
| `Dev/07_Orbit_Correction/orbit_correction.py` | Implemented and notebook-tested |
| `Dev/08_Envelope/envelope.py` | Implemented and notebook-tested |
| `Dev/09_Aperture/aperture.py` | Implemented and notebook-tested |
| `Dev/10_Tune/tune.py` | Implemented and notebook-tested |
| `Dev/11_Snapshot/snapshot.py` | Implemented and notebook-tested |
| `io/` | Planned next |
| `plotting/` | Planned |
| `utils/` | Planned |
| GUI layer | Later, likely Streamlit |

## Planned Next Layer: `io/`

`io/` should provide repository-local configuration, calibration, machine data and measurement data loading once the current snapshot interface has stabilised. It should not reintroduce personal absolute paths or sibling-repository runtime dependencies.

## Later Planned Layers

- `plotting/`: orbit, branch comparison, envelope, aperture margin and tune diagrams. Plotting should consume DataFrames/results rather than running MAD-X directly.
- `utils/`: config hashing, DataFrame comparison helpers, path/log helpers and other small shared utilities.

## Deferred Cleanup Notes

- Older copied tune plotting layers in `Dev/06_Tune_Matching/`, `Dev/07_Orbit_Correction/`, `Dev/08_Envelope/` and `Dev/09_Aperture/` predate the current plotting-interface rule. If those layers are edited again, update their tune plotting helpers to consume repository DataFrame/result outputs directly, including sliced DataFrames, rather than requiring callers to unpack arrays.
- Older validation notebooks in `Dev/06_Tune_Matching/`, `Dev/07_Orbit_Correction/` and `Dev/09_Aperture/` still contain historical markdown provenance references to `isis_2024`. They do not currently represent executable cross-repository reads, but should be rewritten to describe copied repo-local data and internal provenance when those notebooks are next refreshed.

## Development Principles

- Keep layers small and testable.
- Do not modify lattice files as part of normal GUI/backend operation.
- Use generated machine-state override files for reproducible MAD-X runs.
- Preserve both controls-facing values and MAD-X-facing values where useful.
- Keep `cpymad`/MAD-X execution inside `MadxModel` and later workflow wrappers.
- Do not substitute dummy, synthetic, fake or hand-built model outputs for cpymad/MAD-X-backed behaviours unless explicitly requested as toy examples.
- Prefer DataFrame-based interfaces between layers.
- Avoid importing large old helper modules into new clean layers unless a specific function is deliberately migrated.
- Add or update validation notebooks before relying on a new layer.
- Keep the code usable from Jupyter notebooks.

## Example Current Workflow

Run this from `Dev/03_MADX_Model/`, where the relative lattice paths used by the current notebooks resolve correctly:

```python
from cycle_time import RCSRamp
from machine_state import MachineState
from machine_state_writer import write_machine_state_file
from madx_model import MadxModel

ramp = RCSRamp(top_energy_MeV=800.0)
beam_state = ramp.state_at(3.0)

machine_state = MachineState.from_defaults(
    beam_state=beam_state,
    main_magnet_mode="rcs_bare",
    requested_qx=4.331,
    requested_qy=3.731,
    tune_method="di_wright",
)

state_file = write_machine_state_file(
    machine_state,
    output_dir="./machine_states",
)

model = MadxModel(
    lattice_folder="../Lattice_Files/00_Simplified_Lattice/",
    sequence_name="synchrotron",
    machine_state_file=state_file,
    output_dir="./madx_runs/default",
)

model.load_lattice(use_sequence=True)
twiss_df = model.run_twiss()
summary_df = model.get_summary_df()
```
