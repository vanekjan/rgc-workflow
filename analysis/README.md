# Analysis

This folder contains downstream physics analysis code using skimmed ROOT files produced by the skim workflow.

The analysis code is organized by physics channel, not forced into one generic analysis program.

## Layout

```text
analysis/
├── common/
├── configs/
├── inclusive/
│   ├── qe/
│   ├── dis/
│   └── resonance/
├── sidis/
└── dihadron/
```

## Channel organization

### inclusive/

Inclusive electron analyses using the inclusive skim tree.

This includes:

```text
inclusive/qe/          Quasi-elastic region analysis
inclusive/dis/         DIS region analysis
inclusive/resonance/   Resonance-region analysis
```

### sidis/

Single-hadron SIDIS analysis.

This is for channels such as:

```text
e pi+
e pi-
e K+
e K-
e proton
```

### dihadron/

Dihadron analysis, especially electron plus pi+ pi- final states.

### common/

Shared helper code used by more than one analysis channel.

Examples:

```text
ROOT file loading helpers
run-list helpers
target-map helpers
charge and FCup normalization helpers
binning helpers
plot helpers
constants
```

### configs/

Shared analysis configuration files.

Channel-specific configs may also live inside each channel folder later.

## Notes

Generated ROOT files, logs, and plots should normally not be committed unless they are small final reference outputs.

Large analysis outputs should stay outside Git or be ignored.