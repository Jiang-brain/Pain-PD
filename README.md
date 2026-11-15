# Pain Burden, Plasma Proteomics, and Incident Parkinson Disease

Analysis code for the study of pain burden, genetic risk, plasma proteomics, and
incident Parkinson disease (PD).

## Analysis Overview

- **Figure 4:** Composite pain-score weights were estimated among participants
  without proteomic data and applied to the proteomic cohort. The analysis
  identified 24 shared proteins. Single-protein and joint mediation used
  inverse odds-ratio weighting with 95% bootstrap confidence intervals.
- **Figure 5:** Elastic Net (L1 and L2 regularization) selected a stable panel
  of 370 proteins for LightGBM prediction. Performance was evaluated using
  out-of-fold predictions. Panel C covers diagnoses after five years, with no
  ten-year cutoff. No PCA was used in the final analyses.

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
