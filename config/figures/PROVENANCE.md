# Manuscript Figures 4–6

The final author archive is `TF_and_intervention_prioritization__5_.zip` (SHA-256
`110aaddc5302f04d69333247edcd72d9d0e2d7d0051c6996466b7a361b54739d`).
Its `evaluation.tex` includes the three PDF files below. Their SHA-256 hashes
match the PDFs in the original figure directories. The PDF includes identify
the manuscript figures when the archive also has a PNG export from an earlier
revision. The archived manuscript and software environments are external to
this repository.

| Manuscript asset | Original generated PDF | SHA-256 of original PDF | Renderer in this repository |
|---|---|---|---|
| `figures/fig4_budget.pdf` | `figure_section5_2_budget10/bold_typography_larger/SECTION5_2_BUDGET10_PANEL_C_BOLDER.pdf` | `3feb9365c26f21961a0d98deae695bf2c91d3a13442265aef7e317c2be4b05de` | `scripts/figures/render_budget10_panel_c_bolder.py` |
| `figures/fig5_pairs.pdf` | `section4_rank_distribution_20261120_115703/mirrored_bold_typography_larger/PAIR_RANK_DISTRIBUTIONS_MIRRORED_BOLD_LARGER.pdf` | `c86844e02413d18feffa461a68a64b5a3a24c7dceb0655222073937de95e17db` | `scripts/figures/render_mirrored_larger.py` |
| `figures/fig6_single.pdf` | `hammelman_top10_count_rank_qc_20261120_165652/TOP10_HIT_COUNT_AND_RANK.pdf` | `dbb99f1c3bfe2485731505dd3063d4d342e60635af956ea41647225b55bf39ad` | `scripts/figures/plot_top10_count_and_rank.py` |

The original renderers retain their panel dimensions, data selection and plot
geometry. The public adaptations change file paths, output names, input
validation paths, PDF/SVG metadata and system font discovery. The Figure 6
input digest follows the public status table, which omits private run paths.
The output table of hits, ranks and status matches the original figure table
byte for byte. The renderer checks the Figure 4 values, Figure 5 class
medians and AUROC, and Figure 6 hatch glyphs against the frozen figure data.
With the original system font faces, the rebuilt PDFs matched the archived
PDFs pixel for pixel when rendered with Poppler at 125 dpi. A separate run
with the DejaVu fallback passed the numerical and layout checks.

Original renderer file SHA-256 digests before the public path adaptation

| Original file | SHA-256 |
|---|---|
| `plot_section5_2_budget10.py` | `98173588b6abf5095ac113980d12d3176e784a2be1beff4740fff6f679c7c99e` |
| `render_budget10_panel_c_bolder.py` | `01ab7eaf49f58f01dc82334977205752b6a316398722b415051b5eeb25e205de` |
| `render_mirrored_larger.py` | `48fb886518c1504a9cf16eac8a643fcf3dec48f624c31e8adf0cdc266349fe11` |
| `render_rounded.py` | `290dace1f5df2c73f9d91a1d05a664d7104e9afe88e2f5c630ee92a04c824f75` |
| `plot_top10_count_and_rank.py` | `305deaa16653b94284319d6e731326b9f6eda280f72c30efea900ea99ddde61c` |

| Input | Origin |
|---|---|
| `data/figures/SECTION5_2_BUDGET10_TOP12.tsv` | Main Figure 4 source directory; exact archived plot rows |
| `data/pairs/manuscript/budget10_first12.tsv`, `budget10_figure_summary.tsv`, and `mirrored_panel_metrics.tsv` | Author archive `data/pairs`, used for numerical identity checks |
| `data/figures/PLOT_SLOTS.tsv` and `PANEL_METHOD_SUMMARY.tsv` | Frozen rank-distribution reconstruction; plot-ready Figure 5 tables also rebuilt by the pair evaluator |
| `data/figures/TOP10_HIT_COUNT_AND_RANK_DATA.tsv` | Final Figure 6 source directory, with plot data expected after recomputation |
| `config/figures/fig4_source_hashes.tsv` | Five relevant source rows of the frozen Tsunemoto result manifest |

Run after installing the quick reproduction dependencies from the repository root.

```sh
python scripts/figures/render.py --output results/figures
```

Use `--rank-tables` for a directory containing `PLOT_SLOTS.tsv` and
`PANEL_METHOD_SUMMARY.tsv` from a fresh pair evaluation. Each figure is
exported as PDF, editable SVG and a PNG at 300 dpi or higher. Figure 4
used Arial Narrow Bold Italic for gene labels and Figure 5 used Arial Bold.
The scripts use the licensed copies installed on the system if available,
or fall back to Matplotlib's DejaVu family. To select a directory with
licensed font files, set `FIGURE_FONT_DIR` to its path. Font fallback or
rendering engine changes can alter raster hashes while preserving ranks,
axis positions and the plotted values. No fonts are redistributed.

The author archive's numerical verification script `data/pairs/verify_figures.py`
(SHA-256 `83d2e5a7b98cdb42d12b45ee326346c7d8ddad813bfdbfac4cf3d6ebe4dcef51`) is
included as `scripts/figures/verify_figures.py` with configurable paths and
unchanged calculations. The expected compressed plot rows were absent from
the manuscript ZIP. The adapted verifier reads the Figure 5 renderer's
`MIRRORED_SLOT_RANKS.tsv` and the repository's frozen initial pair ranks in
the same schemas. Run it after the standalone renderer to regenerate
the two paired-screen summaries independently.

```sh
python scripts/figures/verify_figures.py --mirrored-ranks results/figures/MIRRORED_SLOT_RANKS.tsv --output results/figure_checks
```
