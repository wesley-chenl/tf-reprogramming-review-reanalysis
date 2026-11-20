# TF reprogramming review reanalysis

```sh
python -m pip install -r requirements.txt && python reproduce.py
```

Run from the repository root in a Python 3.12 environment. The install step needs access to Python packages. Once the environment is ready, `python reproduce.py` runs offline using only the bundled inputs. Use `--output PATH` to place the regenerated tables and figures elsewhere.

## Reproduce the manuscript numbers

| Analysis | Frozen inputs | Quick-run outputs |
| --- | --- | --- |
| §5.2. Measured TF-pair ranking in the Tsunemoto screen | `data/pairs` | Fixed-budget yield and hits, matched-panel comparisons, partner controls and sensitivity tables in `build/tables/pairs` |
| §5.4. Historical single-TF recovery in three Hammelman comparisons | `data/single` | Full eligibility and rank checks, reference matches, CompForce components and QC-exclusion tables in `build/tables/single` |
| Figures 4–6 | `data/figures` and frozen analysis inputs | `Figure_4`, `Figure_5` and `Figure_6` in `build/figures`, each as SVG, PDF and 300 dpi or higher PNG |

The quick run checks bundled input hashes, recomputes the evaluation from frozen ranks, compares against the manuscript tables, and runs the original plotting code with portable paths. Figure 6 uses the final version with ten QC hatches. Font availability can change text layout on different systems. See [verification](VERIFICATION.md) for the cold-start test and [source records](SOURCES.md) for provenance and external licenses.

## Rescore from public inputs

The full workflow uses method-specific environments. Keep scoring workspaces separate from the screen outcomes or historical references. Prepare the data and scoreable sets, run each scoring method, freeze the scores, eligibility and tie order, then join labels or references for evaluation. The quick workflow starts at that last stage.

1. **Tsunemoto pairs.** Follow [the pair workflow](docs/PAIRS.md) for accessions, expression and chromatin inputs, scoring, eligibility, score locking, and evaluation. The paired comparisons use the relevant physical-slot intersection before re-ranking.
2. **Hammelman single TFs.** Follow [the single-TF workflow](docs/SINGLE.md) for the six methods, the CompForce components, and the staged RNA and ATAC inputs. Match reference genes only after freezing rankings. The QC-exclusion comparison removes the prespecified flagged method–target configurations from the shortlist union.

The full re-scoring stages require separate downloads and were not re-run during packaging. This repository contains the portable scoring and evaluation code, frozen parameters, hashes and small derived inputs. See the method-specific workflow for which stages need external files and their exact layout.

## Citation and release

Cite the review and the relevant original datasets and methods. Annotations are specific to the stated experimental panels. Historical reference recovery, screen-positive retrieval, and measured yield each answer a distinct question. The repository is a local release candidate. The authors will supply the repository license, citation metadata, release tag and archive DOI.
