# Revision Summary: Observed S2 Distribution and Outlier Association

## Purpose

Address the reviewer concern that spin contamination may affect the transition-state energetics and the BEP residual outliers. The updated analysis now includes both TS logs and the open-shell reactant/component logs used in `3_Build_DataBase.ipynb`.

The analysis is implemented in:

- Notebook: `output/jupyter-notebook/revision_spin_contamination_outlier_association.ipynb`
- Executed notebook: `output/jupyter-notebook/revision_spin_contamination_outlier_association.executed.ipynb`
- TS spin cache: `output/jupyter-notebook/spin_contamination_ts.csv`
- TS/outlier analysis table: `output/jupyter-notebook/spin_contamination_outlier_association.csv`
- Reactant/component spin table: `output/jupyter-notebook/reactant_component_spin_contamination.csv`
- Reactant/component summary: `output/jupyter-notebook/reactant_component_spin_summary.csv`
- Reaction-level reactant spin table: `output/jupyter-notebook/reaction_reactant_spin_contamination.csv`
- Reaction-level reactant spin correlations: `output/jupyter-notebook/reaction_reactant_spin_correlations.csv`
- High-TS-contamination case table: `output/jupyter-notebook/high_ts_spin_contamination_cases.csv`

Figures:

- `Figure/Revision_SpinContamination_Distribution.png`
- `Figure/Revision_SpinContamination_Outlier_Association.png`
- `Figure/Revision_ReactantComponent_SpinContamination_Distribution.png`
- `Figure/Revision_ReactantSpin_Outlier_Association.png`

## Log Lookup

The TS log lookup follows `3_Build_DataBase.ipynb`: search `E:/work/B_Cl_Nu/Sum/TS_needIRC`, first using the base filename and then the conformer-specific filename.

| item | value |
|---|---:|
| TS reactions | 9237 |
| TS logs found | 9237 |
| TS log coverage | 100.0% |
| normal termination | 9237 / 9237 |
| conformer-specific filenames | 9114 |
| base filenames without conformer suffix | 123 |

For reactants/components, the notebook parses the open-shell logs used by `3_Build_DataBase.ipynb`:

- `B_radical`: `Data/GS_OPT/B_single` and `Data/GS_SPE/B_single`
- `BN_complex_r`: `Data/GS_OPT/B_N_r` and `Data/GS_SPE/B_N_r`
- `C_radical_from_Cl_p`: `Data/GS_OPT/Cl_p` and `Data/GS_SPE/Cl_p`

All 9192 open-shell component log rows were found.

## TS S2 Distribution

All TS calculations are doublets, for which the ideal value is `<S^2> = 0.75`. The values below are the observed `<S^2>` values directly, without subtracting 0.75.

| metric | pre-annihilation `<S^2>` | after-annihilation/final `<S^2>` |
|---|---:|---:|
| mean | 0.7643 | 0.750153 |
| median | 0.7642 | 0.750100 |
| 90th percentile | 0.7705 | n/a |
| 95th percentile | 0.7722 | 0.750400 |
| 99th percentile | 0.7766 | n/a |
| maximum | 0.8323 | 0.7526 |

The pre-annihilation TS values span `0.7500-0.8323`, with 95% at or below `0.7722`. Only 16 of 9237 TS structures have `<S^2> >= 0.80`, corresponding to 0.17%.

### High-TS-Contamination Tail

The 16 TS structures with pre-annihilation `<S^2> >= 0.80` show a clear substrate pattern:

| feature | observation |
|---|---:|
| high-S2 TS rows | 16 |
| unique chloride substrates among high-S2 TS rows | 1 |
| chloride substrate | `Cl_Index = 488`, `CN1CN(C)C(Cl)=C1Cl` |
| total rows with `Cl_Index = 488` | 213 |
| high-S2 rows among `Cl_Index = 488` | 16 |
| high-S2 fraction within `Cl_Index = 488` | 7.51% |
| high-S2 fraction in full TS dataset | 0.17% |
| most common LB among high-S2 rows | `N_Index = 143`, `N#CP(C#N)C#N`, 6 / 16 |
| most common B radical among high-S2 rows | `B_Index = 417`, `Bc1c(F)c(F)c(F)c(F)c1F`, 4 / 16 |

Thus, the high-S2 tail is not broadly distributed across all components. It is completely concentrated in reactions involving `Cl_Index = 488`. The B radical and Lewis base identities vary, although `N_Index = 143` and `B_Index = 417` are enriched within this small tail. Half of the 16 high-S2 TS rows are also BEP 2-sigma outliers, so a high TS `<S^2>` value is not equivalent to the BEP outlier label.

## Reactant/Component S2 Distribution

The reaction energy in `Data/TS/Borane_all.csv` is constructed from `B_N_r + Cl_r`. `Cl_r` is closed-shell, so the reaction-level open-shell reactant is `BN_complex_r`.

Observed pre-annihilation `<S^2>` by open-shell species:

| species class | stage | unique species | median | 95th percentile | maximum | count at or above 0.80 |
|---|---|---:|---:|---:|---:|---:|
| B radical | GS_OPT | 55 | 0.7519 | 0.7668 | 0.7782 | 0 |
| B radical | GS_SPE | 55 | 0.7527 | 0.7735 | 0.7844 | 0 |
| B-LB complex_r | GS_OPT | 4363 | 0.7623 | 0.7773 | 0.8096 | 6 |
| B-LB complex_r | GS_SPE | 4363 | 0.7691 | 0.7948 | 0.8494 | 107 |
| C radical from Cl_p | GS_OPT | 178 | 0.7556 | 0.7794 | 0.7884 | 0 |
| C radical from Cl_p | GS_SPE | 178 | 0.7567 | 0.7939 | 0.8056 | 3 |

For reaction-level association, the `BN_complex_r / GS_OPT` value was mapped to all 9237 reactions:

| metric | value |
|---|---:|
| mapped reaction rows | 9237 / 9237 |
| median BN-reactant `<S^2>` | 0.7649 |
| 95th percentile BN-reactant `<S^2>` | 0.7797 |
| maximum BN-reactant `<S^2>` | 0.8096 |
| reaction rows at or above 0.80 | 43 |

## BEP Outlier Definition

The global BEP fit was recomputed from the released final TS table:

`DeltaG_act = 0.460 * DeltaG_rxn + 30.257`

Outliers were defined as reactions with absolute BEP residual at least two residual standard deviations.

| item | value |
|---|---:|
| residual standard deviation | 4.372 kcal/mol |
| 2-sigma cutoff | 8.744 kcal/mol |
| BEP outliers | 445 |
| BEP outlier fraction | 4.82% |

## Association With BEP Outliers

Observed TS `<S^2>`:

| group | n | median TS `<S^2>` | 95th percentile | count at or above 0.80 | fraction at or above 0.80 |
|---|---:|---:|---:|---:|---:|
| non-outlier | 8792 | 0.7643 | 0.7721 | 8 | 0.091% |
| outlier | 445 | 0.7622 | 0.7753 | 8 | 1.80% |

Observed `<S^2>` of the actual open-shell reactant (`BN_complex_r / GS_OPT`):

| group | n | median BN-reactant `<S^2>` | 95th percentile | count at or above 0.80 | fraction at or above 0.80 |
|---|---:|---:|---:|---:|---:|
| non-outlier | 8792 | 0.7649 | 0.7797 | 42 | 0.478% |
| outlier | 445 | 0.7645 | 0.7788 | 1 | 0.225% |

Correlation with absolute BEP residual:

| metric | Pearson | Spearman |
|---|---:|---:|
| TS `<S^2>` vs absolute BEP residual | -0.0019 | -0.0037 |
| BN-reactant `<S^2>` vs absolute BEP residual | 0.0502 | 0.0771 |
| max(TS, BN-reactant) `<S^2>` vs absolute BEP residual | 0.0297 | 0.0484 |

Interpretation:

- TS `<S^2>` values are tightly concentrated near the ideal doublet value of 0.75 and show essentially no monotonic relationship with BEP residual magnitude.
- The actual open-shell reactant, `BN_complex_r`, shows a similarly narrow GS_OPT `<S^2>` distribution.
- The BEP outlier group does not show elevated BN-reactant `<S^2>`. Its median and fraction at or above 0.80 are slightly lower than those of the non-outlier group.
- Therefore, elevated `<S^2>` in either the TS or the open-shell reactant is not a systematic source of the BEP residual outliers.

## Recommended Reviewer-Response Text

We added a new analysis of the observed `<S^2>` values for all 9237 transition-state calculations and for the open-shell reactant/component logs used in the database construction. Following the same log-file lookup used in `3_Build_DataBase.ipynb`, all TS Gaussian logs were located and parsed. For these doublet TS calculations, the ideal `<S^2>` value is 0.75. The pre-annihilation values are tightly distributed around this reference, with a median of 0.7642, a 95th percentile of 0.7722, and a range of 0.7500-0.8323. After annihilation, the median is 0.7501 and the maximum is 0.7526.

We also parsed the open-shell reactant/component logs for the B radical, the B-LB radical complex, and the chloride-derived carbon radical. The actual open-shell reactant entering the reaction energy is the `B-LB complex_r`; its GS_OPT `<S^2>` was mapped to all 9237 reactions. This reactant also remains close to the doublet reference, with a median of 0.7649, a 95th percentile of 0.7797, and a maximum of 0.8096.

To test whether elevated `<S^2>` is associated with the BEP residual outliers, we recomputed the global BEP residuals and defined outliers as reactions with absolute residuals greater than two residual standard deviations. The absolute BEP residual and TS `<S^2>` show essentially no correlation (Pearson r = -0.0019, Spearman rho = -0.0037). The corresponding correlation for the BN-reactant `<S^2>` is also weak (Pearson r = 0.0502, Spearman rho = 0.0771). The BEP outlier group does not show elevated BN-reactant `<S^2>` compared with the non-outlier group (median 0.7645 vs 0.7649). Therefore, elevated `<S^2>` in either the TS or the open-shell reactant is not a systematic source of the BEP residual outliers.

## Possible Manuscript/SI Sentence

The observed `<S^2>` values were extracted from the Gaussian logs for all 9237 transition states and the corresponding open-shell B-LB reactants. The pre-annihilation values were tightly distributed near the ideal doublet value of 0.75 for both TSs (median 0.7642, 95th percentile 0.7722) and B-LB reactants (median 0.7649, 95th percentile 0.7797), while the annihilation-corrected TS values were essentially ideal (median 0.7501). No meaningful correlation was observed between absolute BEP residuals and either TS `<S^2>` (Spearman rho = -0.0037) or reactant `<S^2>` (Spearman rho = 0.0771), indicating that elevated `<S^2>` is not a systematic origin of the BEP residual outliers.
