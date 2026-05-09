# dlsimtools

Automation framework for DL_MONTE and DL_POLY molecular simulations.

## Installation

```bash
pip install dlsimtools
```

## Overview

`dlsimtools` provides a high-level automation layer for running and analysing DL_MONTE (Monte Carlo) and DL_POLY (Molecular Dynamics) simulations, including tools for HPC job submission, input file manipulation, and output analysis.

### Modules

**Simulation control**
- `tmmc_master_control.Controller` — top-level controller for running full TMMC (Transition Matrix Monte Carlo) pipelines, including range-seeking, bias optimisation, and production runs. Supports local execution and HPC job submission (Slurm/Archer2/Isambard3) via the `sched` parameter.
- `MonteCore` — low-level automation of DL_MONTE runs: editing CONTROL/CONFIG files, launching serial/parallel jobs, and monitoring convergence.
- `LSMC` — range-seeking and optimised bin calculation for TMMC simulations.
- `PolyCore` — equivalent automation utilities for DL_POLY MD runs.

**Input/output**
- `InputConverter` — converts between simulation input formats.
- `MonteCon` — reads and manipulates DL_MONTE CONFIG files.
- `FieldTools` — handles FIELD file editing.
- `CifTools` — reads CIF crystallographic files and prepares simulation inputs.
- `PolyDataCore`, `PolyOutput` — DL_POLY output parsing and data extraction.
- `STATIS` — parses DL_POLY STATIS files.
- `TrajAnalysis` — reads and analyses XYZ trajectory files.

**Thermodynamics & sampling**
- `SwitchBias` — switch bias calculation for TMMC.
- `LSMC` — lattice-switch Monte Carlo utilities.
- `CoexistenceMP` — coexistence point calculation.
- `chempotfind` — chemical potential search.
- `MetaSurf`, `MetaUtil` — metadynamics surface utilities.

**Utilities**
- `GeneralUtil` — common file editing and unit conversion utilities.
- `GeneralOptimizer` — general-purpose optimisation routines.
- `OptimizerScript` — scripting interface for optimisation workflows.
- `HPCworker` — writes and submits job scripts for Slurm, Archer2, and Isambard3.

## Quick start

```python
from dlsimtools.tmmc_master_control import Controller

pc = Controller(
    dlm_exec="/path/to/DLMONTE-SRL.X",
    dlm_exec_par="/path/to/DLMONTE-PRL.X",
    nw=16,
    nodes=1,
    sched="archer2",
)

# Run full TMMC pipeline at 300K (requires CONTROL, CONFIG, FIELD in working directory)
pc.master_surf_tmmc(temp=300)

# Or submit jobs for multiple temperatures to HPC
pc.tmmc_looper(temps=[300, 350, 400], timec=24)
```
