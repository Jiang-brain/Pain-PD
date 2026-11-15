# Pain Burden, Plasma Proteomics, and Incident Parkinson Disease

Analysis code for the study of pain burden, genetic risk, plasma proteomics, and
incident Parkinson disease (PD).

## Data Availability and Access

This repository contains analysis code only. Individual-level participant data
and proteomic measurements from UK Biobank, the China Health and Retirement
Longitudinal Study (CHARLS), the Parkinson's Progression Markers Initiative
(PPMI), and the Olink-based proteomic datasets used in this study are not
included or redistributed.

Access to these data is controlled by the respective data providers and is
subject to their application procedures, data-use agreements, and participant
confidentiality requirements. Researchers wishing to reproduce the analyses
must apply independently to the relevant data custodians and obtain the
necessary approvals. The authors of this repository cannot provide, transfer,
or authorize access to these datasets.

After obtaining authorized access, users may place the required local input
files under `data/` or set the environment variable `PD_PAIN_DATA_DIR` to
another secure, private data directory. The filenames listed below refer to
locally prepared analysis files and are not distributed with this repository.

## Repository Layout

| Directory | Contents |
| --- | --- |
| `00_setup/` to `03_proteomics/` | Setup, cohort preparation, epidemiology, and protein association analyses |
| `04_proteomics_mediation/` | Independent pain scoring, 1,000-bootstrap score stability, and shared-protein mediation |
| `05_prediction/` | Elastic Net, Cox prediction analyses, and LightGBM |
| `06_figures/` | Figure scripts |
| `figure_data/` | Aggregate enrichment inputs for Figure 4G |

## Running the Code

Use R 4.4.x with `data.table`, `survival`, `rms`, `ggplot2`, and `jsonlite`.
Use Python 3.12+ and install dependencies with
`python -m pip install -r requirements.txt`. Figure 4 and 5 plotting requires Arial.

Run from the repository root. After cohort preparation and association analyses
(`01` to `03`), run the independent pain-score workflow in
[`04_proteomics_mediation/`](04_proteomics_mediation/), then its mediation
workflow and the [prediction workflow](05_prediction/README.md). See each
workflow README for required inputs and commands. The score-stability workflow
is separate and optional for reproducing Figures 4 and 5.

```bash
python 06_figures/02_plot_final_figure4.py
python 06_figures/03_plot_final_figure5.py
```

These commands plot saved results without refitting models. Analysis outputs
default to `results/`; the Figure 4 and 5 scripts accept `--output-dir`.

## Data Access

Individual-level UK Biobank, CHARLS, PPMI, and Olink proteomic data are not
included. Access must be obtained from the respective data custodians under
their approval and data-use requirements. Keep authorized inputs in a private
`data/` directory or set `PD_PAIN_DATA_DIR`. Use
`PD_PAIN_PREDICTION_SOURCE_DIR` for private reference inputs where required.
Do not upload participant-level inputs, intermediate files, or predictions.
