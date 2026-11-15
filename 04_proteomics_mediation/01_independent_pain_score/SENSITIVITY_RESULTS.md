# Non-Proteomic Weight Derivation: Sensitivity Results

## Material Passport

- Mode: requested analysis implementation, execution, and reproducibility checks.
- Status: COMPLETED; all 16 implemented output-integrity checks passed.
- Source: local final PA scripts and data, with no external data transfers.
- Scope: new pain-score weights, pain-protein associations, and shared proteins.
- Main analyses, manuscript, and Figure 4: not modified or overwritten.

## Cohorts and Score Development

| Cohort | Participants | PD events |
| --- | ---: | ---: |
| Original eligible score-development cohort | 457,639 | 3,742 |
| No measured proteins: new weight derivation | 408,583 | 2,995 |
| Measured proteins: application of fixed weights | 49,056 | 747 |

The union of unfiltered and QC-filtered protein IDs contained 53,039 unique
participants. None entered the weight-derivation group. Eligibility for score
development further restricted the proteomic scoring group to 49,056 people.
The subsequent pain-protein analysis retained 48,586 base complete cases, exactly
as in the original PA analysis; protein-specific sample sizes also matched for
every one of the 2,920 proteins.

The 14 pain coefficients were estimated only in the 408,583 non-proteomic
participants. Following the original PA code, negative coefficients were set to
zero and the sum of the resulting coefficients was fixed at 1.19553693636065.
The same fixed weights and denominator were applied to the proteomic group.
No PD outcomes from the proteomic group were used to estimate these new weights.
Original and sensitivity scores correlated at r = 0.9701604 in the 49,056 scored
proteomic participants. This correlation is descriptive, not a validation metric
for future disease prediction.

## Protein Associations and Shared Proteins

| Result | Original PA analysis | Sensitivity |
| --- | ---: | ---: |
| Proteins tested | 2,920 | 2,920 |
| Pain-protein P < 0.05 | 980 | 1,055 |
| Pain-protein BH FDR < 0.05 | 780 | 834 |
| Pain-protein Bonferroni < 0.05 | 359 | 398 |
| Protein-PD BH FDR < 0.05 | 205 | 205 (reused) |
| Both associations FDR < 0.05, any direction | 128 | 135 |
| Both FDR < 0.05 and directionally concordant | 25 | 24 |

The shared-set criterion was `pain_beta * protein_PD_log_HR > 0`, in addition to
FDR < 0.05 for both associations. The original 25-protein set was reproduced
from the original result tables before comparing membership.

Twenty-four of the original 25 proteins (96%) were retained; no new shared
proteins were added. COMP alone no longer met the pain-protein FDR criterion:
its FDR changed from 0.0232955 to 0.1224054. Its pain-protein beta remained
negative (-0.00335514 originally; -0.00239867 in the sensitivity), so exclusion
was due to statistical significance, not a direction reversal. Its unchanged
protein-PD FDR was 0.00113349.

Retained shared proteins: ADGRG2, CA14, CD209, CDA, CHGA, CST5, CTSB, EGFR, FAP,
GPRC5C, HPGDS, IGSF21, ITGA11, ITGA2, ITGAM, ITGAV, ITGB1, ITGB2, MET, MYOC,
NEFL, PEPD, RAB6A, and VAT1.

These results indicate that most of the original shared-protein selection was
retained when pain-score weights were derived outside the proteomic cohort.
They do not demonstrate causality, external validity, or intervention effects.

## Verification and Scope Limits

- Original scores reproduced for all 457,639 participants: maximum absolute
  difference 4.9960e-16.
- Independent recalculation of the new scores directly from raw pain fields:
  maximum absolute difference 5.5511e-16.
- All 2,920 association models returned finite estimates and P values.
- Derivation/protein overlap, frozen denominator, score range, protein universe,
  per-protein sample sizes, FDR/Bonferroni correction, shared-set direction,
  original shared-set reproduction, and unchanged PD-result hash checks passed.
- The source weight model omits genetic PCs and clips negative coefficients;
  the protein-association model includes 20 genetic PCs. These original PA
  specifications were preserved, not harmonised with reorganized helper code.
- The PA protein-association implementation models pain score as the outcome,
  whereas the current manuscript describes protein level as the outcome. This
  source-to-manuscript discrepancy requires author confirmation; it was not
  silently changed for the sensitivity analysis.
- The score source reproduces 3,742 PD events in its eligible cohort, whereas
  the manuscript reports 3,412 in its primary cohort. Confirm differences in
  analytical scope or data version before reporting these cohort counts in the
  manuscript. This sensitivity deliberately reproduces the saved PA scores.
- Mediation, prediction, PPI, enrichment, and the assembled Figure 4 were not
  rerun or revised. A 24-protein analysis must not be described as a rerun of
  the original 25-protein mediation or prediction results.

## Result Files

Outputs are in `results/sensitivity_nonproteomic_weights/`. The principal
tables are `shared_protein_summary.csv`, `shared_proteins_fdr05_concordant.csv`,
`shared_protein_retained_lost_gained.csv`, `pain_score_fixed_beta_comparison.csv`,
`pain_score_protein_association_results.csv`, and `validation_checks.csv`.
Participant-level outputs in this directory must remain private.
