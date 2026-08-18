# Inclusive QE Tensor Analysis

This directory contains the inclusive quasi-elastic (QE) tensor-analysis workflow for RGC Spring 2023.

The workflow does four things:

1. scans many Spring 2023 QE runs,
2. ranks possible run pairs using charge/yield/cut-flow diagnostics,
3. applies RGC QADB-based run-quality penalties,
4. runs the final two-run QE tensor-asymmetry analysis for selected pairs.

Directory location in the repository:

```text
rgc-workflow/analysis/inclusive/qe
```

---

## 1. Main goal

The purpose of this analysis is to find stable two-run pairs for the inclusive QE tensor-asymmetry analysis.

A useful pair must have more than a large tensor-polarization difference. The pair should also have compatible charge-normalized yields, compatible cut-flow survival, consistent trigger behavior, and acceptable QADB status.

The main pair-selection logic is:

```text
good pair =
large enough DeltaQ
+ Y1/Y2 close to 1
+ reasonable R2Q screening value
+ compatible raw/final yield behavior
+ compatible cut survival
+ compatible trigger behavior
+ clean or acceptable QADB status
+ stable behavior from kin0 to kin1
```

---

## 2. Current recommended pairs

The current clean late-block pairs are:

| Role | Pair | Reason |
|---|---|---|
| Primary | 17780 / 17792 | highest clean ranking and strong tensor lever arm |
| Stability check | 17780 / 17793 | most stable from kin0 to kin1 |
| Cross-check | 17781 / 17792 | clean independent nearby pair |
| Backup check | 17781 / 17793 | optional consistency check |

Current integrated results:

| Pair | Role | DeltaQ | kin0 Y1/Y2 | kin0 R2Q | kin1 Y1/Y2 | kin1 R2Q |
|---|---|---:|---:|---:|---:|---:|
| 17780 / 17792 | primary | 0.13583 | 0.990019 | -0.0732 | 0.979191 | -0.1520 |
| 17780 / 17793 | stability | 0.12515 | 0.992243 | -0.0617 | 0.991672 | -0.0663 |
| 17781 / 17792 | cross-check | 0.131006 | 0.985108 | -0.1130 | 0.981143 | -0.1429 |
| 17781 / 17793 | backup | 0.120326 | 0.987321 | -0.1047 | 0.993648 | -0.0526 |

Use `kin1` as the nominal narrow QE result. Use `kin0` as the broad detector/PID stability control.

---

## 3. Physics quantities

For each run, the tensor polarization is computed from the vector polarization:

```text
Q = 2 - sqrt(4 - 3*P^2)
```

For a run pair:

```text
DeltaQ = |Q1 - Q2|
```

The charge-normalized selected yield is:

```text
Y = N_selected / FCup_gated_charge
```

The integrated two-run tensor-asymmetry screening quantity is:

```text
R2Q = (Y1 - Y2) / (Q1*Y2 - Q2*Y1)
```

`R2Q` is a diagnostic screening quantity for run-pair comparison. The final interpretation should also inspect the binned `R2Q_vs_xB` plots.

---

## 4. Cut-profile convention

The input sample is the broad inclusive QE `kin0` skim. The analysis then applies one of two analysis-level cut profiles.

| Profile | Definition | Purpose |
|---|---|---|
| kin0 | detector/PID baseline only, then xB counting range | broad stability/control profile |
| kin1 | detector/PID baseline plus tight QE kinematic cuts | nominal narrow QE profile |

The `kin1` cuts are the relevant ones for the narrow QE analysis. The `kin0` profile is used to check whether the same run pair is stable before applying tight QE kinematics.

Current `kin1` cuts:

| Variable | Cut |
|---|---|
| theta_e | 7.80 < theta_e < 8.20 deg |
| vz | -5.758 < vz < 1.5165 cm |
| pe | pe > 2.0 GeV |
| Q2 | 1.9433 < Q2 < 2.0574 GeV^2 |
| W | 0.0 < W < 1.073 GeV |
| xB counting range | 0.85 < xB < 1.15 |

The detector/PID baseline cuts are applied in both `kin0` and `kin1`.

---

## 5. Directory contents

Expected directory structure:

```text
analysis/inclusive/qe/
  README.md
  .gitignore
  Makefile
  setup_env.sh

  qe_analysis.cc
  qe_pair_scan.cc

  run_qe_analysis.sh
  run_pair_scan.sh
  run_selected_pairs.sh

  configs/
    pair_scan_kin0.yaml
    pair_scan_kin1.yaml
    qe_analysis_kin0.yaml
    qe_analysis_kin1.yaml

  tools/
    qadb/
      qadb_pair_postprocess.py
      plot_qadb_pair_diagnostics.py
      compare_pair_scan_cut_profiles.py

  do-compare_pair_scan_cut_profiles.mac
  do-run_pair_scan.mac
  do-run_plot_qadb_pair_diagnostics.mac
  do-run_qadb_postprocess.mac
  do-run_qe_analysis_kin0.mac
  do-run_qe_analysis_kin1.mac
```

Generated files go under `outputs/`. They should not be committed to Git.

---

## 6. What each code/script does

| File | Purpose |
|---|---|
| `qe_pair_scan.cc` | Scans many runs and creates candidate run pairs. Computes run-level and pair-level diagnostics, including the base pair-quality score. |
| `qe_analysis.cc` | Runs the final two-run QE analysis for one selected pair. Produces ROOT output, CSV diagnostics, and binned plots. |
| `run_pair_scan.sh` | Main driver for many-run pair scans. Reads YAML, runs `qe_pair_scan`, applies QADB postprocessing, and makes diagnostic plots. |
| `run_qe_analysis.sh` | Main driver for one selected run pair. Reads YAML and runs `qe_analysis`. |
| `run_selected_pairs.sh` | Convenience script to run the selected final pairs in both `kin0` and `kin1`. |
| `tools/qadb/qadb_pair_postprocess.py` | Adds QADB run-quality information to the pair-scan CSV and computes QADB-adjusted pair scores. |
| `tools/qadb/plot_qadb_pair_diagnostics.py` | Makes diagnostic plots from the QADB-merged pair CSV. |
| `tools/qadb/compare_pair_scan_cut_profiles.py` | Compares pair ranking and pair stability between `kin0` and `kin1`. |
| `configs/*.yaml` | Controls input paths, run numbers, cuts, output names, and QADB/plot options. |
| `Makefile` | Builds `qe_analysis` and `qe_pair_scan`. |
| `setup_env.sh` | Loads the CLAS12/ROOT environment on ifarm. |
| `do-*.mac` | Convenience csh/tcsh macros for running common commands. |

---

## 7. Build

From this directory:

```bash
cd /home/utsav/workhallc/rgc/rgc-workflow/analysis/inclusive/qe
source setup_env.sh
make clean
make
```

Expected build products:

```text
qe_analysis
qe_pair_scan
```

These are local executables and should not be committed.

---

## 8. Input skim path

The YAML files currently point to the Spring 2023 broad QE skim:

```text
../../../skim/rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin0tf1/spring23/spring23_qe_full
```

Check the configured input paths with:

```bash
grep -R "skim_dir\|indir\|input_dir\|rootfiles" -n configs
```

Expected examples:

```text
configs/qe_analysis_kin0.yaml:skim_dir: ../../../skim/rootfiles/...
configs/qe_analysis_kin1.yaml:skim_dir: ../../../skim/rootfiles/...
configs/pair_scan_kin0.yaml:indir: ../../../skim/rootfiles/...
configs/pair_scan_kin1.yaml:indir: ../../../skim/rootfiles/...
```

---

## 9. Full workflow

Run all commands below from:

```bash
cd /home/utsav/workhallc/rgc/rgc-workflow/analysis/inclusive/qe
```

### Step 1: Load environment and build

```bash
source setup_env.sh
make clean
make
```

### Step 2: Optional clean start

Use this only when you want to archive old outputs before a fresh run.

```bash
mkdir -p archived_outputs

if [ -d outputs ]; then
  mv outputs archived_outputs/outputs_$(date +%Y%m%d_%H%M%S)
fi

mkdir -p outputs
```

### Step 3: Dry-run the pair-scan configs

```bash
./run_pair_scan.sh --config configs/pair_scan_kin0.yaml --dry-run
./run_pair_scan.sh --config configs/pair_scan_kin1.yaml --dry-run
```

Confirm that the dry run prints the correct input directory:

```text
../../../skim/rootfiles/skim_sidisdvcs_inclusive_qe_eT0pid1kin0tf1/spring23/spring23_qe_full
```

### Step 4: Run the full many-run pair scans

```bash
./run_pair_scan.sh --build --config configs/pair_scan_kin0.yaml
./run_pair_scan.sh --config configs/pair_scan_kin1.yaml
```

This creates outputs for both cut profiles:

```text
outputs/pair_scan/kin0/
outputs/pair_scan/kin1/
```

Expected files:

```text
outputs/pair_scan/kin0/qe_pair_scan_spring23_kin0.root
outputs/pair_scan/kin0/qe_pair_scan_spring23_kin0_run_summary.csv
outputs/pair_scan/kin0/qe_pair_scan_spring23_kin0_pair_summary.csv
outputs/pair_scan/kin0/qe_pair_scan_spring23_kin0_pair_summary_qadb.csv
outputs/pair_scan/kin0/qadb_spring23_run_metadata.csv

outputs/pair_scan/kin1/qe_pair_scan_spring23_kin1.root
outputs/pair_scan/kin1/qe_pair_scan_spring23_kin1_run_summary.csv
outputs/pair_scan/kin1/qe_pair_scan_spring23_kin1_pair_summary.csv
outputs/pair_scan/kin1/qe_pair_scan_spring23_kin1_pair_summary_qadb.csv
outputs/pair_scan/kin1/qadb_spring23_run_metadata.csv
```

Expected diagnostic plots:

```text
outputs/pair_scan/kin0/diagnostics/spring23_kin0_yield_ratio_vs_abs_deltaQ.png
outputs/pair_scan/kin0/diagnostics/spring23_kin0_R2Q_vs_abs_deltaQ.png
outputs/pair_scan/kin0/diagnostics/spring23_kin0_raw_vs_final_yield_ratio.png
outputs/pair_scan/kin0/diagnostics/spring23_kin0_theta_survival_vs_final_yield_ratio.png
outputs/pair_scan/kin0/diagnostics/spring23_kin0_score_vs_abs_deltaQ.png
outputs/pair_scan/kin0/diagnostics/spring23_kin0_top50_qadb_score.png
outputs/pair_scan/kin0/diagnostics/spring23_kin0_qadb_pair_category_counts.png

outputs/pair_scan/kin1/diagnostics/spring23_kin1_yield_ratio_vs_abs_deltaQ.png
outputs/pair_scan/kin1/diagnostics/spring23_kin1_R2Q_vs_abs_deltaQ.png
outputs/pair_scan/kin1/diagnostics/spring23_kin1_raw_vs_final_yield_ratio.png
outputs/pair_scan/kin1/diagnostics/spring23_kin1_theta_survival_vs_final_yield_ratio.png
outputs/pair_scan/kin1/diagnostics/spring23_kin1_score_vs_abs_deltaQ.png
outputs/pair_scan/kin1/diagnostics/spring23_kin1_top50_qadb_score.png
outputs/pair_scan/kin1/diagnostics/spring23_kin1_qadb_pair_category_counts.png
```

### Step 5: Compare `kin0` and `kin1`

```bash
python3 tools/qadb/compare_pair_scan_cut_profiles.py \
  --kin0-csv outputs/pair_scan/kin0/qe_pair_scan_spring23_kin0_pair_summary_qadb.csv \
  --kin1-csv outputs/pair_scan/kin1/qe_pair_scan_spring23_kin1_pair_summary_qadb.csv \
  --outdir outputs/pair_scan/comparison \
  --prefix spring23_kin0_vs_kin1 \
  --pair-class pos_neg \
  --top-n 50
```

Expected outputs:

```text
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_pair_comparison.csv
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_score_comparison.png
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_score_comparison.pdf
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_yield_ratio_comparison.png
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_yield_ratio_comparison.pdf
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_R2Q_screen_comparison.png
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_R2Q_screen_comparison.pdf
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_rank_shift_topN.png
outputs/pair_scan/comparison/spring23_kin0_vs_kin1_rank_shift_topN.pdf
```

### Step 6: Run selected final pairs

```bash
./run_selected_pairs.sh
```

This runs the current selected pairs in both `kin0` and `kin1`.

Expected output directories:

```text
outputs/pairs/kin0/17780_17792/
outputs/pairs/kin1/17780_17792/

outputs/pairs/kin0/17780_17793/
outputs/pairs/kin1/17780_17793/

outputs/pairs/kin0/17781_17792/
outputs/pairs/kin1/17781_17792/

outputs/pairs/kin0/17781_17793/
outputs/pairs/kin1/17781_17793/
```

Expected files in each pair/profile directory:

```text
qe_Azz_approx_RUN1_RUN2_PROFILE.root
qe_Azz_approx_RUN1_RUN2_PROFILE_run_diagnostics.csv
qe_Azz_approx_RUN1_RUN2_PROFILE_pair_diagnostics.csv
qe_Azz_approx_RUN1_RUN2_PROFILE_R2Q_vs_xB.png
qe_Azz_approx_RUN1_RUN2_PROFILE_R2Q_vs_xB.pdf
qe_Azz_approx_RUN1_RUN2_PROFILE_Azz_approx_vs_xB.png
qe_Azz_approx_RUN1_RUN2_PROFILE_Azz_approx_vs_xB.pdf
qe_Azz_approx_RUN1_RUN2_PROFILE_Azz_times2_diagnostic_vs_xB.png
qe_Azz_approx_RUN1_RUN2_PROFILE_Azz_times2_diagnostic_vs_xB.pdf
```

---

## 10. Quick output checks

Check the pair-scan outputs:

```bash
find outputs/pair_scan -maxdepth 3 -type f | sort | head -50
```

Check the selected-pair diagnostic CSVs:

```bash
find outputs/pairs -name "*_pair_diagnostics.csv" | sort
```

Check the selected-pair binned `R2Q` plots:

```bash
find outputs/pairs -name "*R2Q_vs_xB.png" | sort
```

For the current selected pairs, there should be eight pair diagnostic CSVs and eight `R2Q_vs_xB.png` plots:

```text
kin0 and kin1 for 17780_17792
kin0 and kin1 for 17780_17793
kin0 and kin1 for 17781_17792
kin0 and kin1 for 17781_17793
```

---

## 11. Base pair-quality score

The many-run scanner first computes a base pair-quality score.

This score is only a ranking metric. It is not the physics asymmetry.

The base score rewards:

| Component | Preferred behavior |
|---|---|
| Tensor lever arm | large `DeltaQ` |
| Final yield ratio | `Y1/Y2` close to 1 |
| Raw-yield ratio | raw run behavior close between the two runs |
| Cut-survival ratios | similar survival through detector and kinematic cuts |
| Trigger compatibility | same or consistent dominant trigger behavior |

The tensor-lever-arm factor is:

```text
LQ = min(|DeltaQ| / 0.15, 2.0)
```

The value `0.15` is a chosen reference scale for ranking. It is not a physics constant.

Example for pair 17780 / 17792:

| Quantity | Value |
|---|---:|
| Q1 | 0.188357 |
| Q2 | 0.0525271 |
| DeltaQ = abs(Q1 - Q2) | 0.13583 |
| LQ = DeltaQ / 0.15 | 0.9055 |

The simplified base-score structure is:

```text
base pair score =
100
* tensor-lever-arm factor
* final-yield compatibility factor
* raw-yield compatibility factor
* cut-survival compatibility factors
* trigger-compatibility factor
```

Example for the primary `kin1` pair:

| Component | Input | Factor |
|---|---:|---:|
| DeltaQ lever arm | 0.13583 | 0.9055 |
| Y1/Y2 compatibility | 0.979191 | 0.9792 |
| raw-yield compatibility | 1.0316 | 0.9846 |
| cut/trigger compatibility | stable | 0.9956 |
| Base score | 100 times all factors | 86.919 |

A high base score means the pair has good tensor leverage and good run-to-run compatibility before QADB penalties are applied.

---

## 12. QADB-adjusted pair score

After the base pair-quality score is computed, QADB information is applied as a quality penalty.

```text
QADB-adjusted score = base pair-quality score * QADB penalty
```

QADB penalty table:

| QADB pair category | Penalty | Treatment |
|---|---:|---|
| clean | 1.00 | keep full score |
| cnd_warning | 0.90 | minor down-rank |
| qadb_warning | 0.70 | stronger down-rank |
| target_pol_warning | 0.30 | avoid for final |
| helicity_warning | 0.20 | avoid or inspect manually |
| junk_or_tuning | 0.00 | reject |
| empty_or_special | 0.00 | reject |

Example from the `kin1` ranking:

| Pair | Base score | QADB category | Penalty | QADB-adjusted score |
|---|---:|---|---:|---:|
| 17780 / 17792 | 86.919 | clean | 1.00 | 86.919 |
| 17492 / 17545 | 93.3606 | cnd_warning | 0.90 | 84.0245 |
| 17781 / 17792 | 83.8814 | clean | 1.00 | 83.8814 |
| 17780 / 17793 | 81.0141 | clean | 1.00 | 81.0141 |

QADB does not create the pair score. It only modifies the scanner score so known warning or problematic runs do not dominate the final ranking.

---

## 13. Important outputs

### Pair-scan outputs

| Output | Meaning |
|---|---|
| `qe_pair_scan_spring23_kin*_run_summary.csv` | run-level diagnostics |
| `qe_pair_scan_spring23_kin*_pair_summary.csv` | pair-level diagnostics before QADB adjustment |
| `qe_pair_scan_spring23_kin*_pair_summary_qadb.csv` | pair-level diagnostics after QADB classification and adjusted score |
| `qadb_spring23_run_metadata.csv` | QADB metadata used for run classification |
| `diagnostics/*.png` and `diagnostics/*.pdf` | pair-ranking and diagnostic figures |

### Selected-pair outputs

| Output | Meaning |
|---|---|
| `*_run_diagnostics.csv` | diagnostics for each run in the pair |
| `*_pair_diagnostics.csv` | integrated pair result, including yields and R2Q |
| `*_R2Q_vs_xB.png/pdf` | binned R2Q plot |
| `*_Azz_approx_vs_xB.png/pdf` | binned Azz-approx plot |
| `*.root` | ROOT output for the pair analysis |

---

## 14. Git policy

Commit source code, configs, scripts, and documentation.

Commit these:

```text
README.md
.gitignore
Makefile
setup_env.sh
qe_analysis.cc
qe_pair_scan.cc
run_qe_analysis.sh
run_pair_scan.sh
run_selected_pairs.sh
configs/*.yaml
tools/qadb/*.py
do-*.mac
```

Do not commit generated outputs or local build products:

```text
qe_analysis
qe_pair_scan
outputs/
outputs_old/
archived_outputs*/
qadb_cache/
*.root
*.png
*.pdf
*.log
*.out
*.err
```

Before committing, check:

```bash
git status --short
git diff --cached --stat
git diff --cached --name-only | grep -E 'outputs|qadb_cache|\.root$|\.png$|\.pdf$|\.log$|/qe_analysis$|/qe_pair_scan$'
```

If the last command prints nothing, no generated analysis output is staged.

---

## 15. Minimal command summary

From this directory:

```bash
source setup_env.sh
make clean
make

./run_pair_scan.sh --config configs/pair_scan_kin0.yaml --dry-run
./run_pair_scan.sh --config configs/pair_scan_kin1.yaml --dry-run

./run_pair_scan.sh --build --config configs/pair_scan_kin0.yaml
./run_pair_scan.sh --config configs/pair_scan_kin1.yaml

python3 tools/qadb/compare_pair_scan_cut_profiles.py \
  --kin0-csv outputs/pair_scan/kin0/qe_pair_scan_spring23_kin0_pair_summary_qadb.csv \
  --kin1-csv outputs/pair_scan/kin1/qe_pair_scan_spring23_kin1_pair_summary_qadb.csv \
  --outdir outputs/pair_scan/comparison \
  --prefix spring23_kin0_vs_kin1 \
  --pair-class pos_neg \
  --top-n 50

./run_selected_pairs.sh
```
