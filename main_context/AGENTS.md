# Repository Guidelines

## Project Structure & Module Organisation

This repository is the backend/model layer for an ISIS RCS optics GUI. It is currently organised as development folders rather than a packaged Python module.

- `Dev/01_Cycle_Time/`: bottom-layer beam/ramp model and `01_Test_Cycle_Time.ipynb`.
- `Dev/02_Machine_State/`: current machine-state foundation layer: defaults, tune control, correctors, `MachineState`, writer, and `machine_state_tests.ipynb`.
- `Dev/03_MADX_Model/`: current MAD-X execution layer and validation gate: `MadxModel`, local copies of foundation modules, and `madx_model_tests.ipynb`.
- `Dev/04_Error_Table_Utilities/`: current standalone error-table utility layer: `errors.py`, local copies of foundation/MAD-X modules, and `errors_tests.ipynb`.
- `Dev/05_Orbit_Branches/`: current orbit branch workflow layer: `OrbitBranch`, local copies of foundation/MAD-X/error modules, and `orbit_branch_tests.ipynb`.
- `Dev/06_Tune_Matching/`: current tune matching and harmonic tune workflow layer: tune matching wrappers, tune programme evaluation, tune plotting helpers, local copies of lower-layer modules, and `tune_matching_tests.ipynb`.
- `Dev/07_Orbit_Correction/`: current read-only closed-orbit correction workflow layer: BPM/corrector selection, measured-orbit fitting, MAD-X `CORRECT`, corrector-current suggestions, local copies of lower-layer modules, and `orbit_correction_tests.ipynb`.
- `Dev/08_Envelope/`: current beam-envelope workflow layer: MAD-X TWISS consumption, geometric/normalised RMS emittance handling, n-sigma envelope columns, harmonic envelope comparison plots, local copies of lower-layer modules, and `envelope_tests.ipynb`.
- `Dev/09_Aperture/`: current aperture-margin workflow layer: MAD-X APERTURE extraction, ISIS source aperture spreadsheet loading, envelope/aperture alignment, clearance-margin summaries, local copies of lower-layer modules, copied source aperture data, and `aperture_tests.ipynb`.
- `Dev/10_Tune/`: current tune/resonance helper layer: tune summary extraction, tune programme normalisation, resonance-line/proximity tables, DataFrame-backed tune plotting helpers, local copies of lower-layer modules, and `tune_tests.ipynb`.
- `Dev/11_Snapshot/`: current whole-system snapshot integration layer: single and multi-snapshot orchestration, copy/modify snapshot configs, full-cycle working-point resonance inputs, lower-layer integration tests, local copies of lower-layer modules, and `snapshot_tests.ipynb`.
- `Dev/04_MADX_Model_Revisit/`: revisit/reference copies and generated outputs. Do not treat this as the default edit target unless the task explicitly names it.
- `Dev/Lattice_Files/`: fixed MAD-X lattice inputs.
- `Dev/Error_Tables/`: fixed survey/error-table inputs.

Generated outputs live under `machine_states/`, `machine_states_test/`, and `madx_runs/`.

## Architecture Rules

The old standalone `strength_scaling.py` concept is superseded by the machine-state approach.

- `cycle_time.py` produces the beam state.
- `MachineState` stores GUI/backend state for one cycle time and must not run MAD-X.
- `tune_control.py` stays a pure helper for supported tune-control calculations.
- `correctors.py` converts controls-style currents and MAD-X kicks.
- `machine_state_writer.py` writes executable MAD-X override/log snapshots.
- `MadxModel` owns `cpymad`/MAD-X execution, direct `MachineState` application, machine-state file application, error-table application, TWISS, and MAD-X tune matching.
- `errors.py` owns pure MAD-X error-table DataFrame/path utilities and may delegate application to `MadxModel`; it must not own MAD-X execution.
- `orbit_branch.py` owns branch-level orchestration, orbit DataFrame extraction, branch comparison and compact orbit summaries for notebooks and the future GUI. It must use isolated `MadxModel` instances and must not solve orbit correction.
- `tune_matching.py` owns tune workflow orchestration, Di Wright actual-tune evaluation, MAD-X tune-matching wrappers, full-cycle tune-programme tables and single-point tune evaluation. It must keep trim quad currents and MAD-X K values available for GUI display.
- `tune_plots.py` owns resonance, harmonic trim-quad, beta-variation and trim-current plotting helpers. It must consume DataFrames/results and must not run MAD-X.
- `orbit_correction.py` owns read-only COCU-style orbit-correction workflows. It must fit measured BPM data through MAD-X/cpymad, run MAD-X `CORRECT` for selected HD/VD steering correctors, preserve BPM/corrector enable masks, and convert proposed kicks to currents via ISIS corrector calibrations. It must not write to the machine.
- `envelope.py` owns beam-envelope calculations from real MAD-X TWISS tables and explicit user beam assumptions. It may use `MadxModel` workflow wrappers to obtain TWISS, but plotting helpers must only consume DataFrames/results. It must support geometric and normalised RMS emittance, arbitrary positive sigma scale, momentum-spread contribution, and GUI-ready nominal/harmonic envelope comparison outputs. Envelope plots must render the plus and minus envelope bounds for the same result with identical line width, style and colour.
- `aperture.py` owns aperture-table normalisation, source-spreadsheet aperture loading, envelope/aperture alignment, clearance-margin calculation, limiting-location summaries and plotting. It may use `MadxModel.run_aperture` to obtain real MAD-X APERTURE tables, but plotting helpers must only consume DataFrames/results. It must handle zero MAD-X aperture rows from drifts by excluding them from margin calculations and must include the ISIS source aperture spreadsheet data in plots and tests. Aperture notebooks and wrappers should start from repository-local `Dev/Lattice_Files/00_Simplified_Lattice` unless a task explicitly compares another lattice.
- `tune.py` owns tune summary extraction, tune programme normalisation, resonance-line generation, resonance proximity and tune-diagram input tables. It must keep set tunes distinct from predicted/actual MAD-X tunes and preserve DataFrame-ready outputs for plotting and GUI use.
- `snapshot.py` owns whole-system orchestration for one or more machine snapshots. It may build `MachineState` objects and run MAD-X through `MadxModel` or existing workflow wrappers, but it must not duplicate lower-layer physics logic or plotting logic. It must support single snapshots, copied/modified snapshot configs, snapshot series, nominal-vs-error comparisons, full-cycle working-point resonance inputs, GUI-ready table accessors, metadata and warnings.

Do not implement dummy, synthetic, fake, placeholder, or manually fabricated physics/model outputs unless the user explicitly asks for a toy example. If a method represents machine optics, orbit, tune, matching, correction, TWISS, error response, or another MAD-X-backed behaviour, implement the complete cpymad/MAD-X call path or leave the method unimplemented with a clear error. The future GUI depends on real model behaviour; synthetic stand-ins make it useless.

Keep repository data and source dependencies internal. Do not make code, notebooks, tests, examples, or validation paths read from sibling repositories, personal checkout paths, or other absolute external repo paths. If values or static data from another repository are needed, copy them into this repository or into the relevant notebook/file and use only repo-local references, with a concise provenance note if useful.

Plotting helpers should consume the DataFrame/result formats produced by this repository's workflow and helper functions directly, including sliced DataFrames where practical. Do not require callers to manually unpack function outputs into separate arrays just to plot them; keep any column selection or normalisation inside the plotting helper or a shared repo-local adapter.

Do not modify lattice or source error-table files as part of normal backend work. `MadxModel` deliberately does not load `2023.strength` or `ISIS.injected_beam`; beam state is supplied through the generated machine-state override file or direct `MachineState` application.

## Build, Test, and Development Commands

There is no package manifest or Makefile. Run scripts and notebooks from their own development folder so local imports resolve.

Syntax check examples:

```bash
python -m py_compile Dev/02_Machine_State/*.py
python -m py_compile Dev/03_MADX_Model/*.py
python -m py_compile Dev/06_Tune_Matching/*.py
python -m py_compile Dev/07_Orbit_Correction/*.py
python -m py_compile Dev/08_Envelope/*.py
python -m py_compile Dev/09_Aperture/*.py
python -m py_compile Dev/10_Tune/*.py
python -m py_compile Dev/11_Snapshot/*.py
```

Validation notebooks:

```bash
jupyter nbconvert --to notebook --execute --inplace Dev/01_Cycle_Time/01_Test_Cycle_Time.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/02_Machine_State/machine_state_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/03_MADX_Model/madx_model_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/04_Error_Table_Utilities/errors_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/05_Orbit_Branches/orbit_branch_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/06_Tune_Matching/tune_matching_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/07_Orbit_Correction/orbit_correction_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/08_Envelope/envelope_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/09_Aperture/aperture_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/10_Tune/tune_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace Dev/11_Snapshot/snapshot_tests.ipynb
```

If a task requires verifying notebooks without changing generated repository outputs, copy the relevant `Dev/` folders to `/tmp` and execute the copies there.

Expected dependencies include `numpy`, `pandas`, `matplotlib`, `cpymad`, Jupyter, and a working MAD-X installation. Add Streamlit commands only once a GUI entry point exists.

## Coding Style & Naming Conventions

Use standard Python style with 4-space indentation and `snake_case` names. Keep dataclasses for structured state objects. Use `OrderedDict` where magnet, harmonic, or corrector ordering matters.

Keep `cpymad` calls behind `MadxModel` or future workflow wrappers. Do not place MAD-X execution in `MachineState`, `tune_control.py`, plotting helpers, or future Streamlit page callbacks.

## Testing Guidelines

The current foundation test gate is notebook-based:

- `01_Test_Cycle_Time.ipynb` validates the ramp and beam-state layer.
- `machine_state_tests.ipynb` validates defaults, tune control, correctors, `MachineState`, writer output, JSON sidecars, and selected failure modes.
- `madx_model_tests.ipynb` validates lattice loading, generated and direct machine-state application, loaded-file assertions, TWISS outputs, error tables, MAD-X tune matching, and failure modes.
- `errors_tests.ipynb` validates error-table normalisation, validation, read/write round trips, combination, filtering, sign flipping, TWISS-name mapping, `MadxModel` integration, wrapper delegation, failure modes, and single-element orbit responses for QD, QF, QDS and DIP errors.
- `orbit_branch_tests.ipynb` validates default, misaligned, corrector and combined-error branches, orbit extraction, comparison metrics, isolated branch outputs and selected failure modes.
- `tune_matching_tests.ipynb` validates harmonic input normalisation, explicit harmonic lattice support, MAD-X tune matching, Di Wright actual-tune evaluation, full-cycle set-vs-actual tune programmes, trim quad current/K exposure, harmonic expected-vs-programmed TWISS checks, combined harmonic cases, resonance plots, beta variation plots and harmonic trim-quad plots.
- `orbit_correction_tests.ipynb` validates read-only COCU-style orbit correction: BPM and corrector enable masks, dummy BPM measurements fitted into viable measured orbits by MAD-X `MATCH` using non-steering fit knobs, MAD-X `CORRECT` using selected HD/VD steering correctors, current conversion via ISIS calibrations, measured/corrected orbit plots, corrector-current plots, and final assertions that `CORRECT` reduces full-monitor RMS.
- `envelope_tests.ipynb` validates MAD-X-backed envelope evaluation from `MachineState` TWISS tables, default 300 pi mm mrad RMS emittances, 0.2% momentum spread, arbitrary sigma scale, geometric and normalised RMS modes, large-harmonic envelope comparisons, plotting helpers, and selected failure modes.
- `aperture_tests.ipynb` validates real MAD-X APERTURE extraction, ISIS source aperture spreadsheet loading, zero/invalid MAD-X aperture filtering for drift rows, envelope/aperture alignment, clearance margins, limiting-location summaries, aperture/envelope/margin plots, harmonic and closed-orbit-distortion cases, and selected failure modes.
- `tune_tests.ipynb` validates tune summary extraction, set-vs-predicted tune programme tables, resonance line/proximity generation, DataFrame-backed tune diagram plotting, harmonic trim-quad plots and selected failure modes.
- `snapshot_tests.ipynb` validates whole-system integration: single snapshots, copied/modified configs, nominal-vs-error series, branch comparison, full-cycle working-point resonance inputs, lower-layer table accessors, and existing plotting helpers consuming snapshot/series outputs directly.

Before implementing a new layer, run the three validation notebooks or state clearly why a notebook was not run. For Python-only utilities, add focused assertions in the relevant notebook until a formal test suite exists.

## Planned Layer Order

Next planned implementation target: `io/`.

Later planned layers include `plotting/`, `utils/`, and finally the Streamlit GUI layer. Do not claim planned files are implemented until they exist in the repo.

## Commit & Pull Request Guidelines

Recent commits use short bracketed scopes such as `<add>[Dev] ...`, `<update>[Dev/04] ...`, and `[Dev/02]`. Follow that pattern with a concise imperative summary, for example:

```text
<update>[Dev/02-03] add foundation notebook assertions
```

Pull requests should state the affected `Dev/` folders, list validation notebooks or commands run, and call out generated artifacts intentionally changed.

## Agent-Specific Instructions

Only edit the files requested by the user. Do not edit Python modules when the task is documentation-only. Do not edit notebooks when the task is source-only.

Do not remove or rewrite generated MAD-X or machine-state outputs unless explicitly requested. Do not edit `Dev/Lattice_Files/` or `Dev/Error_Tables/` in normal backend work. Avoid broad refactors across older development iterations; keep changes scoped to the current layer and copy consistency changes only when requested.

Avoid committing `.ipynb_checkpoints/`, `__pycache__/`, generated logs, or timestamped machine-state outputs unless the user explicitly wants those artifacts refreshed.
