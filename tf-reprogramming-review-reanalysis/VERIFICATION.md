# Verification

## Cold start

Copied the release tree and independently extracted the release ZIP into new temporary paths on a different drive from the historical experiment directories. Used CPython 3.12.14, NumPy 2.4.4, Matplotlib 3.11.2 and Pillow 13.0.0. A temporary Python audit hook blocked reads and directory access under the old experiment root for the quick runner and its three child processes. A deliberate legacy file probe confirmed the guard was active. The full quick command exited successfully. ZIP checks verified the CRCs and SHA-256 of all archived files against the local release tree.

```sh
python -m pip install -r requirements.txt && python reproduce.py
python scripts/figures/verify_figures.py --mirrored-ranks build/figures/MIRRORED_SLOT_RANKS.tsv --output build/figure_verify
```

| Check | Recomputed result | Output |
| --- | --- | --- |
| Figure 4 | 533 eligible pairs, 69 positives overall and 3 positives in the first ten. The first ten averaged 0.011% gated yield. Rank 11 yielded 1.60%. | `build/tables/pairs/FIG4_RECOMPUTED.tsv` |
| Figure 5 | A, B and C contain 210, 90 and 65 measured pairs. The respective positive counts are 37, 7 and 18. Class medians, tie-aware AUROC and fixed-budget counts agree with the manuscript. | `build/tables/pairs/PANEL_METHOD_SUMMARY.tsv` and `BOUNDARY_TIES.tsv` |
| TFcomb eligibility | 533/69 expression-eligible, 230/45 after graph coverage, then 65/18 after requiring both coefficients to be positive. The old common-panel intersection contains 60 slots. | `build/tables/pairs/TFCOMB_ELIGIBILITY.tsv` |
| Paper pair diagnostics | The TFcomb and expression bHLH-only diagnostic ranking files reproduce byte for byte. Budget, partner and sensitivity statistics are regenerated from the frozen inputs. | `build/tables/pairs/diagnostic` and corresponding summary tables |
| Single-TF score recovery | Both restored gzip files decompress to the SHA-256 digests in the original score lock and checked archive. Rebuilt top-ten identities and positions, primary eligibility counts and exact reference ranks agree with the manuscript. Table 6 includes the five CompForce components. | `data/single/SCORE_LOCK.json`, `build/tables/single/TABLE6_RANKS.tsv` and `TABLE6_COUNTS_AND_RANDOM_EXPECTATIONS.tsv` |
| Figure 6 | 18 plotted cells, 10 QC hatches and two hatched zero-hit cells. The latter are heart IMCell singleton and β-TC6 ANANSE. After QC exclusion, liver recovers 6/11 and β-TC6 recovers 3/5. No heart configuration remains. | `build/tables/single/TOP10_QC_EXCLUSION.tsv` and `build/figures/Figure_6.svg` |

Figures 4–6 opened as SVG, PDF and PNG. At 125 dpi with the matching system fonts, the regenerated PDFs rendered pixel for pixel like the three PDF inputs in the final manuscript ZIP. The comparison used the PDF sources recorded in [figure provenance](config/figures/PROVENANCE.md). The repository also supports a bundled Matplotlib font fallback.

## Staged commands

The following staging commands also ran against empty directories. Both produced score-only workspaces. The pair label attachment check failed as expected before a score lock was present. Bash syntax checks passed for the included shell scripts.

```sh
python scripts/pairs/prepare_workspace.py score work/pairs_preflight
python scripts/pairs/prepare_workspace.py labels work/pairs_preflight
python scripts/single/stage_workspace.py init work/single_preflight
```

The pair utility and core scoring scripts were also replayed in their offline scoring sandboxes using frozen derived expression inputs. Their score SHA-256 values matched the original run. The quick run uses frozen scores and derived records. The public downloads, raw RNA quantification, ATAC alignment and peak calling, network reconstruction, full method fitting and Monte Carlo scoring were documented for staged replay. They were not re-run from raw inputs during packaging. See [pair stages](docs/PAIRS.md) and [single-TF stages](docs/SINGLE.md) for external files and tool requirements.
