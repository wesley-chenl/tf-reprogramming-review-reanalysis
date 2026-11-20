# Section 5.2 pair analysis

Run `python scripts/pairs/reproduce.py --out build/tables/pairs` from the repository root for the frozen-input replay. It verifies the Figure 4 quantities, reproduces the 210-, 90- and 65-pair tables for Figure 5, and evaluates the manuscript sensitivity analyses and both bHLH controls. The original frozen initial and follow-up slot ledgers are compressed without changing their contents. Rank 1 is the highest priority. The comparison unit is a physical two-factor position, with identity variants and exclusions recorded before scoring. Use `variant=main` for the manuscript comparisons. The experiment's own expression series GSE112381 was excluded from scoring.

The binary endpoint is the published screen call. Gated yield is the printed corrected percent of starting MEFs for positives, and zero by definition for negatives. The raw S1 field retains measured zero and unavailable values separately. Both native and matched budgets use same-panel reference distributions. Fixed-bHLH means weight the uniform mean among each nominated backbone's eligible partners by its number of nominations. Supplemental output keeps frozen integer tie breaks, midranks, and attainable and expected counts for ties crossing the budgets. The exact means, rank statistics and selection tables are derived afresh by the quick script. The original utility evaluator below also reruns the recorded Monte Carlo intervals.

## Source and version records

| Payload | Download and staging record |
| --- | --- |
| Tsunemoto screen | [Source ledger and SHA256](../data/pairs/manuscript/outcome_sources.tsv) for DOI [10.1038/s41586-018-0103-5](https://doi.org/10.1038/s41586-018-0103-5), Figure 1b, the original S1 workbook and identity evidence. The archival program extracts calls and uses the frozen interpretation variants. |
| GSE53951 expression | Joint RMA of GSM1303947/49 (source) and GSM1303943/45 (target). Mouse GPL1261 probe annotation and original R package versions are in [the package table](../config/pairs/utility/EXPRESSION_R_PACKAGES.tsv). R 4.4.3, Bioconductor 3.20, reference BLAS and LAPACK 3.12.x. |
| GSE114809 expression | Author-provided log2 sample matrix. The [expression source ledger](../data/pairs/provenance/SOURCE_EXPRESSION.tsv) records the raw slots, exact file names and SHA256. [The sample sheet](../config/pairs/utility/EXPRESSION_SAMPLE_MANIFEST.tsv) contains all selected samples, exclusions and published titles. |
| GSE43916 expression and GSE101397 ATAC | GSM1074345/46 and GSM1074357. RNA and SOFT files are in the [Taiji metadata source table](../data/pairs/provenance/TAIJI_METADATA_SOURCES.tsv). Selected paired SRA runs SRR5822253/2288/2289/2290/2291, layouts and public SRA/FASTQ locations are in the [run ledger](../data/pairs/provenance/GSE101397_RAW_RUNS.tsv). |
| TFcomb | [MIT source at d80164b](https://github.com/Chen-Li-17/TFcomb/tree/d80164bef403d567206c0127d802f9dce90aaba8), v1.0.1. Original Ridge and native helper adapter use Python 3.11.18, NumPy 2.2.6, pandas 3.0.6, SciPy 1.15.2 and scikit-learn 1.6.1. [Source hashes](../data/pairs/manuscript/software_sources.tsv) and [run receipt](../config/pairs/initial/tfcomb_audit_GSE53951.json). The fitted predictor is a bulk-data adaptation of the published module. |
| TRRUST | Mouse v2 release 2018-04-16, signed evidence by 2015-12-31. The 3,964-edge derived layer and its extracted diagnostic edge columns use CC BY-SA 4.0. [Attribution and terms](../data/pairs/TRRUST_ATTRIBUTION.md), [source ledger](../data/pairs/provenance/SOURCE_TRRUST.tsv). |
| CollecTRI | [Zenodo 8192729](https://zenodo.org/records/8192729) v2.0 and the frozen HCOP 2026-10-01 snapshot. Modern and direction-supported through-2015 mouse one-to-one projections. [Prior protocol](../config/pairs/followups/ADDENDUM_PROTOCOL.yaml). |
| IMCell and SFA | [IMCell v0.1.0](https://zenodo.org/records/20563020), commit `cfa3a869138f9d8f6395ac3509f4ff8a6423a3f9` (MIT), R 4.4.3 and igraph 2.3.4. [SFA](https://github.com/dwgoon/sfa/tree/4ff875e5272ba190ed8c684c97da77471ebb4d21) commit `4ff875e5272ba190ed8c684c97da77471ebb4d21` (BSD 3-Clause). SFA uses the registered sparse 48-term Neumann response with alpha 0.5 and verifies convergence against the independent solve. |
| Taiji and msVIPER | [Taiji](https://github.com/Taiji-pipeline/Taiji/releases/tag/v1.3.0) v1.3.0 and the [author reprogramming example](https://github.com/Wang-lab-UCSD/Taiji-reprogram/tree/0b6337aae2bebdcb5f1db35031ecc38e8b49179a). Bowtie2 2.5.5, SRA Toolkit 3.4.1, samtools and htslib 1.24, MACS3 3.0.4. [Genome and asset checksums](../config/pairs/taiji/STATIC_RESOURCES.tsv). VIPER 1.40.0 RELEASE_3_20 e859341 uses the Columbia license and must be fetched separately. |

The source manifests identify upstream software and data. No third-party code archive is bundled. The original scripts, path adaptations and local sources are itemized in [ORIGINS.tsv](../data/pairs/ORIGINS.tsv).

## Stage 1. Public inputs and outcomes

The commands below use Bash on Linux or WSL and an ordinary UID. Run from the repository root. Set `pair` and `old` to new directories. Set `runtime` to the separately installed native environment where shown. The quick run and the fixture replays use frozen inputs without SRA downloads.

```bash
export PAIR_STAGE="$(pwd)/build/pair-source"
mkdir -p "$PAIR_STAGE"/{validation,inputs_scoring,design_public,sources_private}
python scripts/pairs/utility/prepare_expression.py raw
python scripts/pairs/utility/prepare_expression.py rna
Rscript --vanilla scripts/pairs/utility/prepare_expression_rma.R "$PAIR_STAGE"
python scripts/pairs/utility/prepare_expression.py microarray
cp config/pairs/utility/TRRUST_MIXED_SIGN_REVIEW.tsv "$PAIR_STAGE/inputs_scoring/"
python scripts/pairs/utility/prepare_trrust_prior.py --root "$PAIR_STAGE"
```

Copy the frozen `SCORING_POLICY.json` to `$PAIR_STAGE/design_public/`. To rebuild the positional outcome sheet separately, place the sources named in `prepare_outcomes.py` under `$PAIR_STAGE/sources_private` following the exact downloaded names and checksums from the outcome ledger. Fetch the publisher and NIHMS Figure 1 images, the S1 workbook, PMC XML and frozen mouse Gene Info. Use the included mapping conflict and orthology records as manual identity inputs. With NumPy, Pillow and openpyxl installed, run `python scripts/pairs/utility/prepare_outcomes.py`. This preparation writes `design_public/CLEAN_PAIRS.tsv` for scoring and puts the positional outcomes in `sources_private` for use after locking. Alternatively use the curated frozen outcome TSVs. A changing NCBI Gene/HCOP snapshot can alter mapping and will not reproduce the frozen SHA256.

## Stage 2. Replay the label-blind utilities and core methods

```bash
repo="$(pwd)"
pair="$repo/build/pair-rescore"
old="$repo/build/pair-initial"
python scripts/pairs/prepare_workspace.py score "$pair"
# To use prepared sources instead, append --source "$PAIR_STAGE" to the previous command.
bash "$pair/code/run_scoring.sh" "$pair" "$pair/scores"
python "$pair/code/seal.py" "$pair" SCORE
python scripts/pairs/prepare_workspace.py initial "$old" --source "$pair"
bash "$old/code/run_core.sh" "$old" "$old/scores/core"
```

The first sandbox gives the scorer only two expression contrasts, the locked signed graph and the unlabeled design and policy. The second replays the historical expression and direct scores and computes two-step signed coverage and SFA. Its fixture checks validate exact historical score bytes. Both replays were tested against the archived results with the package's frozen derived inputs. To rerun with a different version of a public snapshot, choose a fresh work directory and record the new input and score hashes.

## Stage 3. Native methods and follow-up scores

Download and unpack the exact author releases to private local directories. Prepare an IMCell helper directory containing upstream `IMCell.R`, `IMCell_setup.R` and the included indexed compatibility layer from `scripts/pairs/followups/imcell_weights/imcell_indexed_cpu.R`. `authors_tfcomb` contains only the upstream `TFcomb/tools/tf_inference.py` and `utils.py`. Verify hashes against the source records. The Linux runtime prefix must contain the appropriate Python and R executable and pinned packages, with packages installed inside the mounted prefix.

```bash
run="$repo/scripts/pairs/run_isolated.sh"
mkdir -p "$old/scores"/{tfcomb,imcell}
for ep in GSE53951 GSE114809; do
    bash "$run" "$old/inputs" "$old/code" "$old/scores/tfcomb" "$runtime" "$authors_tfcomb" \
        "$runtime/bin/python" -B /code/native_modules_tfcomb.py /inputs /out /authors "$ep"
done
bash "$run" "$old/inputs" "$old/code" "$old/scores/imcell" "$runtime" "$authors_imcell" \
    /usr/bin/bash -c '
        python -B /code/imcell_prepare.py /inputs /out &&
        Rscript --vanilla /code/imcell_pair.R /inputs /out /authors 12 score &&
        python -B /code/imcell_expand.py /out
    '
# Transfer from the isolated IMCell output directory to the freezer's expected location.
install -m 644 "$old/scores/imcell/imcell.tsv" "$old/scores/imcell.tsv"
awk 'FNR==1 && NR!=1 {next} {print}' "$old"/scores/tfcomb/tfcomb_GSE*.tsv > "$old/scores/native_modules.tsv"
```

TFcomb uses signed one-, two- and three-hop weights `0.8, 0.15, 0.05` with Ridge alpha 1, and the upstream directing-score PCC of the two retained direct response columns. The adapter records every failed gate. After reveal on the primary contrast, expression retained 533/69, measured outgoing support retained 230/45 and dual positive coefficients retained 65/18 (pairs/positives). The manuscript diagnostic works from fitted coefficients without refitting and restores all 65 full-model scores at 12 decimal places.

For IMCell the six blocks use seed `481991 + block_index`, L'Ecuyer-CMRG and three independent 1,000-trial replicates. To replay the joint/summed single-seed weight sensitivity, copy `mc_block_1.tsv` and `mc_block_2.tsv` into a copy of `$old/inputs`, and `imcell_jobs.tsv`, `imcell_target_program.tsv` and `imcell_expand.tsv` into the weight output directory. Run the following applications in order, with those directories mounted as `/inputs` and `/out` using `run_isolated.sh` and the script folder mounted as `/code`.

```bash
Rscript --vanilla /code/imcell_weights.R /inputs /out /authors 12 initial
python -B /code/imcell_weights_expand.py /inputs /out initial
Rscript --vanilla /code/imcell_weights.R /inputs /out /authors 12 extend
python -B /code/imcell_weights_expand.py /inputs /out final
```

The rule extends both endpoints to 3 × 5,000 when specified by the blinded budget boundary indicators. The original settings and registered constants are saved under `config/pairs`.

Prepare the CollecTRI layers with the author preparation script and a separate public-source directory. Place HCOP at `$audit/sources/annotation/human_mouse_hcop_fifteen_column.txt.gz`. Place `TRRUST_PUBMED_DATES.tsv` under `$PAIR_STAGE/inputs_scoring`. The expression stage provides the frozen Gene Info and `old/inputs` the signed TRRUST cross-reference.

```bash
python scripts/pairs/followups/prior_prepare.py "$prior" "$old" "$PAIR_STAGE" "$audit"
```

For each alternative prior, copy the corresponding prepared edge layer to a blank input directory as `GRAPH.tsv`. Add `old/inputs` tables and `old/scores/core/program.tsv` as `PROGRAM.tsv`. Create a code-only directory with `prior_coverage_score.py`, `prior_viper.R`, `core_score.py` and `historical_score.py`. Use `run_isolated.sh` with those input/code directories to run `python /code/prior_coverage_score.py collectri_modern` or `collectri_through2015` and `Rscript --vanilla /code/prior_viper.R GSE53951 /inputs /out` or `GSE114809`. The two measured-target minima are both five, with `level=2`, `processAll=TRUE`, unit likelihood and `nullmodel=NULL`. Keep each layer's output separate.

For Taiji, choose an independent `PAIR_TAIJI_ROOT` and `PAIR_TAIJI_ASSETS`. Stage the five SRA blocks under `$PAIR_TAIJI_ROOT/taiji/sra/` using the exact selectors in the run ledger. Source run SRR5822253 is one biological unit. Technical target runs SRR5822288/89 and SRR5822290/91 form two biological units. Stage the two GSE43916 files under `metadata/`, the mouse Gene Info file under `expression_raw/`, GENCODE M24 unique gene IDs and symbols as `inputs/annotation/GENES_M24.tsv`, and the blind design under `scoring/common_allowed/`. Copy `ADDENDUM_PROTOCOL.yaml` from `config/pairs/followups`. The analysis assets include GRCm38.p6, GENCODE vM24, the same FASTA-derived Bowtie2 index and Taiji mm10 Index and the mouse cisBP motif file. Match the resource checksums. `TAIJI_BIN` is the versioned executable and `TAIJI_COMPAT_LIB` is its compatibility library folder.

```bash
for sr in SRR5822253 SRR5822288 SRR5822289 SRR5822290 SRR5822291; do
    bash scripts/pairs/followups/taiji_align_run.sh "$sr"
done
bash scripts/pairs/followups/taiji_finish_atac.sh all
python scripts/pairs/followups/taiji_rna.py
bash scripts/pairs/followups/taiji_pagerank.sh
python scripts/pairs/followups/taiji_score.py \
    "$PAIR_TAIJI_ROOT/taiji/TAIJI_NATIVE_GENE_RANKS.tsv" "$PAIR_TAIJI_ASSETS/cisBP_mouse.meme"
```

The registered alignment uses `--very-sensitive-local --seed 0 -I 0 -X 2000` and retains nuclear primary proper pairs with both MAPQ ≥30. Pair-based deduplication is applied within each biological unit. MACS3 uses BAMPE and `-q 0.01 --keep-dup all`. The archived sandbox configs record mm10 and the GeneQuant / NarrowPeak file types. The model requires both members to be increased on the independent RNA contrast and both PageRanks to be positive and finite. Scores multiply the two bilateral PageRank ratios. The expression comparator uses role percentiles of the target log2(TPM+1) minus the mean of source log2(TPM+1).

## Stage 4. Freeze before evaluation

Keep screen labels in a separate directory until each score branch has finished. Record SHA256 of every blind input, config, score and eligibility table. The original utility freeze and evaluator can be run as follows.

```bash
python scripts/pairs/prepare_workspace.py labels "$pair"
# To use rebuilt screen labels instead, append --source "$PAIR_STAGE/sources_private".
python "$pair/code/seal.py" "$pair" LABEL
python "$pair/code/evaluate.py" "$pair" --workers 4
```

For the primary expanded pairs, place the SFA source from the pinned commit in `$old/external/sfa/` for provenance checks, then run `python "$old/code/freeze.py" "$old"` followed by `python "$old/code/evaluate.py" "$old" "$pair"`. The original freeze constructs the expression-plus-IMCell rank integration from the exact component intersection, called before any score evaluation. Freeze the separately scored Taiji and sensitivity outputs before joining them to published positional labels. This publication preserves the completed follow-up ranking artifact. The quick entry recomputes its reported matched evaluations rather than assembling all raw score branches into a new follow-up freeze.

The quick verification is the tested publication path. The public-source prep, native method blocks, legacy Taiji pipeline and complete Monte Carlo replay have additional platform and upstream artifact requirements and were not rerun end to end during publication packaging.
