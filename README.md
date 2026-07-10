# RGC Skim Framework

This repository contains a CLAS12 RGC HIPO skim framework for producing organized ROOT and filtered HIPO skim files from RGC `sidisdvcs` trains.

The framework supports inclusive electron, single hadron SIDIS, and dihadron `e pi+ pi-` skim modes. It also supports target map matching, run selection, period aware output folders, and preservation of `HEL::scaler` information for FCup and helicity normalization studies.

---

## Purpose

The goal of this framework is to provide a reproducible and organized workflow for RGC data skimming.

It is designed to:

- read one HIPO file or a folder of HIPO files
- select events for inclusive, SIDIS, or dihadron analysis
- write compact ROOT skim files
- optionally write HIPO skim files
- preserve scaler information needed for FCup normalized asymmetries
- organize outputs by dataset, channel, cut configuration, run period, and job tag
- allow YAML based configuration for safer and repeatable production

---

## Supported Skim Modes

| Mode | Physics channel | Description |
|---|---|---|
| `inclusive` | `e` only | Selects inclusive electron events and stores DIS, QE, and resonance kinematics |
| `sidis` | `e + h` | Selects electron plus one hadron channel such as `pi+`, `pi-`, `K+`, `K-`, or proton |
| `dihadron` | `e + pi+ + pi-` | Selects events with an electron and at least one `pi+ pi-` pair |

---

## Repository Layout

```text
rgcskim/
├── Makefile
├── README.md
├── setup_env.sh
├── src/
│   └── rgcskim.cc
├── scripts/
│   ├── run_all_hipo.sh
│   └── run_skim.py
├── configs/
│   ├── skim.example.yaml
│   ├── runlists/
│   │   └── README.md
│   └── target_maps/
│       ├── README.md
│       └── polarization_vs_run_results_preliminary_v1__1_(polarization_vs_run_results_pre).csv
```

### Important files

| File | Purpose |
|---|---|
| `src/rgcskim.cc` | Main C++ skim code |
| `scripts/run_all_hipo.sh` | Bash runner for one HIPO file or a full folder |
| `scripts/run_skim.py` | YAML wrapper that builds and runs the bash command |
| `configs/skim.example.yaml` | Example YAML configuration tracked by Git |
| `configs/skim.yaml` | Local working YAML configuration, ignored by Git |
| `configs/target_maps/` | Target and polarization CSV files |
| `configs/runlists/` | Optional run list text files |
| `setup_env.sh` | Environment setup script |
| `Makefile` | Build instructions |

---

## Clone the Repository

Clone the repository from GitHub:

```bash
git clone https://github.com/utsav-sth/rgcskim.git
cd rgcskim
```

This creates a local folder named:

```text
rgcskim/
```

To clone into a custom folder name, use:

```bash
git clone https://github.com/utsav-sth/rgcskim.git your_folder_name
cd your_folder_name
```

---

## Setup

Start from the repository directory:

```bash
cd /path/to/rgcskim
source setup_env.sh
make clean
make
```

A successful build creates the executable:

```text
./rgcskim
```

Useful syntax checks:

```bash
bash -n scripts/run_all_hipo.sh
python3 -m py_compile scripts/run_skim.py
```

---

## Recommended Workflow

The recommended workflow is to use the YAML wrapper.

### Step 1: Create a local YAML config

```bash
cp configs/skim.example.yaml configs/skim.yaml
```

Edit the local config as per your need, see YAML Configuration guide below:

```bash
vi configs/skim.yaml
```

### Step 2: Run a dry run first

```bash
./scripts/run_skim.py configs/skim.yaml --dry-run
```

The dry run prints:

- resolved skim mode
- input path
- output format
- target map settings
- run selection
- predicted output folders
- full `run_all_hipo.sh` command

No skim is executed during a dry run.

### Step 3: Run the skim

```bash
./scripts/run_skim.py configs/skim.yaml
```


---

## YAML Configuration Guide

The main local configuration file is:

```text
configs/skim.yaml
```

A typical inclusive QE test configuration is:

```yaml
runner:
  script: scripts/run_all_hipo.sh

input:
  dataset: sidisdvcs
  path: /path/to/sidisdvcs/
  period: auto

skim:
  mode: inclusive
  region: qe
  pid:

cuts:
  electrontree: 0
  detpidcut: 1

target:
  map: configs/target_maps/polarization_vs_run_results_preliminary_v1__1_(polarization_vs_run_results_pre).csv
  polsource: offline
  targetfilter: matched

run_selection:
  mode: test_one_run
  run:
  start:
  end:
  runs: []
  runlist:

output:
  format: root
  jobtag: test_one_run
```

---

## Input Section

```yaml
input:
  dataset: sidisdvcs
  path: /path/to/sidisdvcs/
  period: auto
```

| Field | Meaning |
|---|---|
| `dataset` | Dataset tag used in output names |
| `path` | HIPO input file or folder containing HIPO files |
| `period` | Output period folder mode |

Allowed period values:

| Value | Meaning |
|---|---|
| `none` | Do not use a period subfolder |
| `auto` | Infer `summer22`, `fall22`, or `spring23` from the input path |
| `summer22` | Force `summer22` output folder |
| `fall22` | Force `fall22` output folder |
| `spring23` | Force `spring23` output folder |

For normal RGC production paths, use:

```yaml
period: auto
```

Example RGC production paths:

```text
/lustre24/expphy/cache/clas12/rg-c/production/summer22/pass1/10.5gev/ND3/dst/train/sidisdvcs/
/lustre24/expphy/cache/clas12/rg-c/production/fall22/pass1/ND3/dst/train/sidisdvcs/
/lustre24/expphy/cache/clas12/rg-c/production/spring23/pass1/ND3/dst/train/sidisdvcs/
```

---

## Skim Section

```yaml
skim:
  mode: inclusive
  region: qe
  pid:
```

### Inclusive mode

```yaml
skim:
  mode: inclusive
  region: qe
  pid:
```

Allowed inclusive regions:

| Region | Meaning |
|---|---|
| `all` | Keep all accepted inclusive electron events |
| `dis` | Keep DIS candidates |
| `res` | Keep resonance candidates |
| `qe` | Keep broad quasi elastic candidates |

### SIDIS mode

Example for `e pi+`:

```yaml
skim:
  mode: sidis
  region:
  pid: 211
```

Common hadron PID values:

| PID | Channel |
|---|---|
| `211` | `e pi+` |
| `-211` | `e pi-` |
| `321` | `e K+` |
| `-321` | `e K-` |
| `2212` | `e proton` |

### Dihadron mode

```yaml
skim:
  mode: dihadron
  region:
  pid:
```

The current dihadron mode is fixed to:

```text
e pi+ pi-
```

---

## Cuts Section

```yaml
cuts:
  electrontree: 0
  detpidcut: 1
```

| Field | Meaning |
|---|---|
| `electrontree` | Writes an extra diagnostic electron tree in SIDIS mode |
| `detpidcut` | Applies loose detector and PID quality cuts |

Allowed values:

```text
0 = off
1 = on
```

For normal analysis skims, use:

```yaml
detpidcut: 1
```

For inclusive and dihadron modes, `electrontree` is forced to `0`.

---

## Event Selection and Cut Summary

This section summarizes the main event selection logic used by the skim. The exact implementation is in `src/rgcskim.cc`.

### Common electron candidate selection

| Cut category | Requirement | Applied when |
|---|---|---|
| Electron PID | `pid == 11` | All modes |
| Inclusive electron status | `-4000 < status <= -2000` and `REC::Particle` row 0 is used when inclusive row 0 mode is enabled | Inclusive mode |
| SIDIS and dihadron electron status | `status / 1000 == -2` | SIDIS and dihadron modes |
| Electron choice | Row 0 for inclusive mode, highest momentum accepted electron for SIDIS and dihadron modes | Mode dependent |

### Loose electron detector and PID cuts

These cuts are applied only when:

```yaml
detpidcut: 1
```

| Cut category | Requirement |
|---|---|
| Charge | `charge == -1` |
| HTCC match | `has_htcc == 1` |
| HTCC photoelectrons | `htcc_nphe > 2.0` |
| Calorimeter match | `has_cal == 1` |
| PCAL energy | `pcal_e > 0.0` |
| Total calorimeter energy | `cal_e > 0.0` |
| Sampling fraction | `sampling_fraction > 0.10` |
| PCAL local coordinates | `lu > 0.0`, `lv > 0.0`, `lw > 0.0` |

### Inclusive region selection

Inclusive kinematics are calculated from the selected electron using the beam energy defined in `src/rgcskim.cc`. The skim stores both nucleon level and deuteron level quantities.

| Quantity | Definition |
|---|---|
| `Q2` | `-q^2` |
| `nu` | Virtual photon energy |
| `W` | Nucleon level invariant mass, used for DIS and resonance separation |
| `Wd` | Deuteron system invariant mass, stored separately |
| `xB` | Nucleon Bjorken x |
| `xD` | Deuteron x |
| `y` | `nu / beam_energy` |

Allowed inclusive regions:

| Region | Requirement |
|---|---|
| `all` | Keep all accepted inclusive electron events |
| `dis` | `W > 2.0` |
| `res` | `0.0 < W < 2.0` |
| `qe` | Broad quasi elastic skim flag with `Q2 > 0.2`, `0.5 < xB < 2.1`, and `0.0 < y < 1.0` |

The `qe` region is intentionally broad at the skim level. Final physics cuts should be applied later in the analysis stage.

### SIDIS selection

SIDIS mode requires at least one accepted electron and at least one selected hadron.

| Selection item | Requirement |
|---|---|
| Electron | `pid == 11` and `status / 1000 == -2` |
| Selected hadron | `pid == selected_hadron_pid` and `status / 1000 == 2` |
| Event topology | At least one accepted electron and at least one selected hadron |
| Output rows | One row is written for each selected hadron paired with the best electron |

When `detpidcut: 1` is used, the hadron must also pass the loose sanity cuts below.

| Hadron cut category | Requirement |
|---|---|
| Beta | `beta > 0.0` |
| Charge consistency | Positive PID requires positive charge, negative PID requires negative charge |

No SIDIS kinematic cuts such as `Q2`, `W`, `x`, `z`, `pT`, or missing mass are applied in the skim.

### Dihadron selection

Dihadron mode currently selects the `e pi+ pi-` channel.

| Selection item | Requirement |
|---|---|
| Electron | `pid == 11` and `status / 1000 == -2` |
| `pi+` | `pid == 211` and `status / 1000 == 2` |
| `pi-` | `pid == -211` and `status / 1000 == 2` |
| Event topology | At least one accepted electron, one `pi+`, and one `pi-` |
| Output rows | One row is written for each accepted `pi+ pi-` pair with the best electron |

When `detpidcut: 1` is used, the electron and pions must also pass the loose detector and PID sanity checks.

### Target and scaler handling

| Item | Behavior |
|---|---|
| Target map | If `targetfilter: matched`, only runs found in the target map are processed |
| Polarization source | `offline` or `online` target polarization can be selected |
| Scaler tree | `HEL::scaler` information is stored for FCup and helicity normalization |
| HIPO output | HIPO output preserves selected physics events and scaler events |


---

## Target Map Section

```yaml
target:
  map: configs/target_maps/polarization_vs_run_results_preliminary_v1__1_(polarization_vs_run_results_pre).csv
  polsource: offline
  targetfilter: matched
```

| Field | Meaning |
|---|---|
| `map` | Target and polarization CSV file |
| `polsource` | Which target polarization value to use for generic branches |
| `targetfilter` | Whether to process all runs or only runs found in the target map |

Allowed `polsource` values:

```text
offline
online
```

Recommended:

```yaml
polsource: offline
```

Allowed `targetfilter` values:

| Value | Meaning |
|---|---|
| `all` | Process all HIPO files |
| `matched` | Process only runs found in the target map |

For ND3 target polarization analysis, use:

```yaml
targetfilter: matched
```

The tracked shared target map is:

```text
configs/target_maps/polarization_vs_run_results_preliminary_v1__1_(polarization_vs_run_results_pre).csv
```

Expected CSV columns:

```text
run,start_time,stop_time,species,cell,charge_avg_online,charge_avg_offline,charge_avg_offline_err,run_dose(Pe/cm2)
```

---

## Run Selection Section

```yaml
run_selection:
  mode: test_one_run
  run:
  start:
  end:
  runs: []
  runlist:
```

Allowed modes:

| Mode | Meaning |
|---|---|
| `all` | Process all HIPO files in the input path |
| `test_one_run` | Process the first HIPO file that passes run and target filters |
| `single` | Process one run |
| `range` | Process an inclusive run range |
| `list` | Process runs listed in the YAML file |
| `file` | Process runs from a text file |

### Test one run

```yaml
run_selection:
  mode: test_one_run
```

This is recommended for checking a new configuration.

### Single run

```yaml
run_selection:
  mode: single
  run: 16270
```

### Run range

```yaml
run_selection:
  mode: range
  start: 16270
  end: 16320
```

### Run list in YAML

```yaml
run_selection:
  mode: list
  runs: [16270, 16272, 16275]
```

### Run list from file

```yaml
run_selection:
  mode: file
  runlist: configs/runlists/test_runs.txt
```

The run list file should contain one run number per line:

```text
16270
16272
16275
```

Blank lines and lines beginning with `#` are ignored.

---

## Output Section

```yaml
output:
  format: root
  jobtag: test_one_run
```

Allowed output formats:

| Format | Meaning |
|---|---|
| `root` | Write ROOT skim only |
| `hipo` | Write filtered HIPO skim only |
| `both` | Write ROOT and filtered HIPO skim |

Use `jobtag` to separate tests from production output.

Example:

```yaml
output:
  format: root
  jobtag: test_one_run
```

For large production jobs, choose a meaningful tag:

```yaml
output:
  format: root
  jobtag: full_summer22_qe
```

---

## Output Folder Structure

Output folders are built from:

```text
dataset + channel + option tag + period + job tag
```

Example inclusive QE output with target matching:

```text
rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1tf1/summer22/test_one_run/
hipofiles/skim_sidisdvcs_inclusive_qe_eT0pid1tf1/summer22/test_one_run/
logs/skim_sidisdvcs_inclusive_qe_eT0pid1tf1/summer22/test_one_run/
```

Example output filename:

```text
sidisdvcs_inclusive_qe_eT0pid1tf1_016270.root
sidisdvcs_inclusive_qe_eT0pid1tf1_016270.hipo
sidisdvcs_inclusive_qe_eT0pid1tf1_016270.log
```

### Option tag

The option tag has the form:

```text
eT<electrontree>pid<detpidcut>tf<targetfilter>
```

Example:

```text
eT0pid1tf1
```

means:

```text
electrontree = 0
detpidcut    = 1
targetfilter = matched
```

---

## Standard Checks Before Production

Before running a larger skim, run:

```bash
make clean
make
bash -n scripts/run_all_hipo.sh
python3 -m py_compile scripts/run_skim.py
./scripts/run_skim.py configs/skim.yaml --dry-run
```

Then run one test file first:

```yaml
run_selection:
  mode: test_one_run
```

After the test succeeds, switch to full production:

```yaml
run_selection:
  mode: all
```

or use a controlled run list.

---

## Recommended Production Procedure

1. Update `configs/skim.yaml`.
2. Run a dry run.
3. Check the resolved command and predicted output folders.
4. Run `test_one_run`.
5. Check the log file.
6. Check the output ROOT or HIPO file.
7. Switch to the intended run selection.
8. Run the production skim.
9. Keep the output folders outside Git.