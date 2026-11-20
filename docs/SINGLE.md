# Section 5.4 single-TF reanalysis

The article reports TRIAGE, VIPER, Taiji, ANANSE, IMCell singleton and CompForce RRA for P0 liver, P1 PCM1-positive cardiomyocyte nuclei and the pancreatic β-TC6 cell line. The native score lock also contains the original DE comparator. The article tables and QC exclusion use the six reported workflows. Mouse Ensembl IDs are matched exactly against each target's documented TF set. The shared universe has 1,374 stable IDs. Hhex is outside that universe. The reference lists derive from Hammelman et al., [Figure 1b](https://www.nature.com/articles/s41592-022-01522-2/figures/1). The lists are historical references for the target and permit no conclusion about single-TF sufficiency.

## Frozen-data analysis

From the repository root, run:

~~~sh
python scripts/single/recompute.py --output reproduced/single
~~~

The two full score gzip files missing from the manuscript archive were recovered from the original Round 2 score lock. Their uncompressed bytes match both the original run outputs and the validated Round 2 light archive. The primary table contains 28,854 rows, including the three legacy DE cells. The CompForce input-scale version contains 4,122 rows. See [the score hashes and archive members](../data/single/SCORE_LOCK.json). The gzip wrappers are new archival encodings of those TSVs. The historical fold and Top-50 summaries can now be recomputed from full scores.

The fast evaluator verifies hashes before opening evaluation references, reorders all eligible scores at full decimal precision and uses lexicographic stable-ID tie breaking. It checks the final manuscript top tens, every displayed reference rank, native CompForce component rank and eligibility count, conditional null expectations and overlaps. It writes Table 6 ranks and counts, exact shortlists, budget unions, common eligibility and the quality-filtered analysis.

Figure 6 reports 18 workflow-target configurations. Its 10 flagged configurations comprise six for cardiomyocytes, two for liver and two for β-TC6. The two flagged zero-hit cells are cardiomyocyte IMCell and β-TC6 ANANSE. The quality exclusion retains no cardiac workflows. It retains four workflows each for liver and β-TC6, recovering six of 11 and three of five documented factors, respectively. Cardiac RNA mapping was 28.615% and 27.912%; pooled MEF ATAC peak overlap was 3.418% of alignment records, and cardiac ATAC peak overlap was 0.030% of alignment records.

## Staged re-scoring from published inputs

Run the raw replay in Linux or WSL on a POSIX volume. The fast entry uses the small files in the repository. This separate entry needs substantial disk for full reads and network outputs. The release verifies workspace initialization without downloads. It preserves the production method runners with path-only changes and separates raw scoring from reference access. Rebuilding the public sequencing reads, installing the archived third-party software and rerunning the native methods remain full-run steps for the reader.

### 1. Prepare the score-only workspace

~~~sh
REPO="$(pwd)"
export HAMMELMAN_WORK="$PWD/work/single"
python "$REPO/scripts/single/stage_workspace.py" init "$HAMMELMAN_WORK"
export HAMMELMAN_BASE="$HAMMELMAN_WORK/base"
export HAMMELMAN_ROUND2="$HAMMELMAN_WORK/round2"
python "$REPO/scripts/single/stage_workspace.py" download "$HAMMELMAN_WORK" --assay reference
python "$REPO/scripts/single/stage_workspace.py" download "$HAMMELMAN_WORK" --assay RNA
python "$REPO/scripts/single/stage_workspace.py" download "$HAMMELMAN_WORK" --assay ATAC
python "$REPO/scripts/single/stage_workspace.py" download "$HAMMELMAN_WORK" --assay software
~~~

The [sample table](../data/single/SAMPLES.tsv) and [public read sources](../data/single/RAW_SOURCE_POINTERS.tsv) specify SRA, GEO and ENCODE accessions, source and target roles, file locations in the replay workspace, byte sizes and checksums. RNA uses source GSM1866625–27 and targets ENCSR096STK, GSM2525582/83 and GSM2666438/40. ATAC uses source GSM1867737/38 and targets ENCSR609OHJ, GSM2525591 and GSM2666437. For cardiac ATAC, the single-end run is SRR5318031. This stage produces a checksum-verified download manifest without historical labels. Large source data and software are obtained separately. The packaged gzip support files contain derived stable-ID mappings and the exact TRIAGE and CollecTRI inputs used by the original scorers.

Use the [base explicit Linux environment](../config/single/base/conda-explicit.txt), [base R session](../config/single/base/R_SESSION.txt), [Round 2 R session](../config/single/ROUND2_R_SESSION.txt), [ANANSE explicit environment](../config/single/ANANSE_CONDA_EXPLICIT.txt) and [ANANSE pip record](../config/single/ANANSE_PIP_FREEZE.txt). Environment paths are selectable. Example:

~~~sh
export HAMMELMAN_BASE_ENV="$HAMMELMAN_WORK/env-base"
export HAMMELMAN_ANANSE_ENV="$HAMMELMAN_WORK/env-ananse"
export HAMMELMAN_R_LIB="$HAMMELMAN_WORK/Rlib"
micromamba create -y -p "$HAMMELMAN_BASE_ENV" --file "$REPO/config/single/base/conda-explicit.txt"
micromamba create -y -p "$HAMMELMAN_ANANSE_ENV" --file "$REPO/config/single/ANANSE_CONDA_EXPLICIT.txt"
export PATH="$HAMMELMAN_BASE_ENV/bin:$PATH"
export R_LIBS_USER="$HAMMELMAN_R_LIB"
mkdir -p "$R_LIBS_USER" "$HAMMELMAN_BASE/environments/sratoolkit"
tar -xzf "$HAMMELMAN_BASE/environments/sratoolkit.tar.gz" -C "$HAMMELMAN_BASE/environments/sratoolkit"
chmod +x "$HAMMELMAN_BASE/environments/bin/taiji"
~~~

The sequencing tools were Salmon 2.11.0, SRA Toolkit 3.4.1, Bowtie2 2.5.5, samtools 1.24, bedtools 2.31.1 and MACS3 3.0.4. R 4.4.3 used Bioconductor VIPER 1.40.0, DESeq2 1.46.0 and tximport 1.34.0. Taiji 1.3.0 may need the compatible libtinfo5 shared library. If needed, extract the library from the Debian package into the work directory, using the location specified by the Taiji wrappers.

The software download verifies and extracts the three [source releases](../config/single/SOURCE_RELEASES.json). Use ANANSE v0.5.1 at `580907581bda42e7bae172bcd346076cb87dcd41`, CompForce v0.1.0 at `c0d64c336cf4fd4506315eb82cc1dac23a92d7dc` (Zenodo 20563029) and IMCell v0.1.0 at `cfa3a869138f9d8f6395ac3509f4ff8a6423a3f9` (Zenodo 20563020). Keep these original sources untouched. The installed working copies used the patches below. Install CellNet 0.1.1 at commit `6049909f2f943aad656186e5bbd1145ec829104f` as a CompForce helper dependency. The [R dependency ledger](../config/single/R_INSTALL_PROVENANCE.json) gives package versions and archive checksums.

~~~sh
# Obtain the pinned CellNet helper and install the working copies.
git clone https://github.com/CahanLab/CellNet.git "$HAMMELMAN_WORK/CellNet-build"
git -C "$HAMMELMAN_WORK/CellNet-build" checkout --detach 6049909f2f943aad656186e5bbd1145ec829104f
PATCHES="$REPO/config/single/patches"
patch "$HAMMELMAN_WORK/CellNet-build/NAMESPACE" < "$PATCHES/CellNet_unused_AnnotationDbi_import.diff"
R CMD INSTALL -l "$R_LIBS_USER" "$HAMMELMAN_WORK/CellNet-build"
"$HAMMELMAN_ANANSE_ENV/bin/python" -m pip install --no-deps --no-build-isolation "$HAMMELMAN_ROUND2/sources/ANANSE"
cp -a "$HAMMELMAN_ROUND2/sources/CompForce" "$HAMMELMAN_WORK/CompForce-build"
cp -a "$HAMMELMAN_ROUND2/sources/IMCell" "$HAMMELMAN_WORK/IMCell-build"
patch "$HAMMELMAN_WORK/CompForce-build/DESCRIPTION" < "$PATCHES/CompForce_optional_dependencies.diff"
patch "$HAMMELMAN_WORK/CompForce-build/R/pm_cellnet.R" < "$PATCHES/CompForce_bounded_exact_pearson.diff"
patch "$HAMMELMAN_WORK/CompForce-build/R/pm_diffexp.R" < "$PATCHES/CompForce_undefined_welch_empty_table.diff"
patch "$HAMMELMAN_WORK/IMCell-build/DESCRIPTION" < "$PATCHES/IMCell_optional_dependencies.diff"
R CMD INSTALL -l "$R_LIBS_USER" "$HAMMELMAN_WORK/CompForce-build" "$HAMMELMAN_WORK/IMCell-build"
~~~

The optional dependency diffs move unused imports to Suggests. The CompForce compatibility diffs preserve the edge-specific Pearson correlation and the Welch result schema, with NA for essentially constant inputs. The original installation also removed CellNet's unused `AnnotationDbi::as.list` import for this matrix-only use. The ANANSE memory adapter is selected by the staged wrapper and preserves native network terms. The IMCell wrapper checks the indexed CPU adapter hash against documented full-result and RNG parity.

### 2. Rebuild RNA and ATAC inputs

~~~sh
python "$HAMMELMAN_BASE/code/quantify_rna.py"
python "$HAMMELMAN_BASE/code/quantify_sra_stream.py"
for target in liver cardiomyocyte betapanc; do
  Rscript --vanilla "$HAMMELMAN_BASE/code/prepare_rna.R" "$HAMMELMAN_BASE" "$target"
done
python "$HAMMELMAN_BASE/code/preprocess_atac.py"
python "$HAMMELMAN_BASE/code/validate_atac_cache.py" --require-all
python "$HAMMELMAN_BASE/code/prepare_atac_condition_peaks.py"
python "$HAMMELMAN_BASE/code/prepare_taiji_genome_index.py" --compact-gc
samtools faidx "$HAMMELMAN_BASE/inputs/annotation/GRCm38.p6.genome.fa"
python "$REPO/scripts/single/stage_workspace.py" freeze-base "$HAMMELMAN_WORK"
~~~

RNA uses GENCODE M24 and a transcriptome-only Salmon k31 index with complete biological reads, library inference, sequence bias and paired-read GC bias. Single-end fragment parameters are mean 200 and SD 80. The fixed target/source Wald contrast uses tximport length-scaled TPM estimated counts, design `~condition` and alpha 0.1. ATAC uses the same mm10 coordinate system, Bowtie2 sensitive-local alignment, seed 0, MAPQ 30, deduplication and removal of mitochondrial reads. MACS3 uses q 0.01 with all remaining duplicates kept. Paired libraries use BAMPE; the single-end cardiac call uses `-f BAM --nomodel --shift -100 --extsize 200`. The pooled source peak object uses the two complete, independently deduplicated source BAMs. The archived flags and exact preprocessing operations are in [the QC evidence](../config/single/qc/ATAC_LAYOUT_AND_ACTUAL_COMMANDS.tsv).

### 3. Score the six workflows and lock full ranks

~~~sh
python "$REPO/scripts/single/stage_workspace.py" prepare-round2 "$HAMMELMAN_WORK"
python "$HAMMELMAN_ROUND2/scripts/prepare_r_inputs.py"
Rscript --vanilla "$HAMMELMAN_ROUND2/scripts/imcell_freeze_jobs.R" "$HAMMELMAN_ROUND2"
"$HAMMELMAN_ANANSE_ENV/bin/python" "$HAMMELMAN_ROUND2/scripts/prepare_ananse.py"
for target in liver cardiomyocyte betapanc; do
  python "$HAMMELMAN_BASE/code/score_A.py" "$target"
  python "$HAMMELMAN_BASE/code/score_R.py" "$target" B
  python "$HAMMELMAN_BASE/code/score_C.py" "$target"
  Rscript --vanilla "$HAMMELMAN_ROUND2/scripts/compforce_run.R" "$HAMMELMAN_ROUND2" "$target"
  Rscript --vanilla "$HAMMELMAN_ROUND2/scripts/imcell_run.R" "$HAMMELMAN_ROUND2" "$target" 6
  for state in source target; do
    "$HAMMELMAN_ANANSE_ENV/bin/python" "$HAMMELMAN_ROUND2/scripts/run_ananse.py" "$target" binding "$state" 4
    "$HAMMELMAN_ANANSE_ENV/bin/python" "$HAMMELMAN_ROUND2/scripts/run_ananse.py" "$target" network "$state" 2
  done
  "$HAMMELMAN_ANANSE_ENV/bin/python" "$HAMMELMAN_ROUND2/scripts/run_ananse.py" "$target" influence task 4
done
python "$REPO/scripts/single/stage_workspace.py" prelock "$HAMMELMAN_WORK"
python "$HAMMELMAN_ROUND2/scripts/rank_round2.py"
~~~

TRIAGE uses the mean biological-replicate `ln(TPM + 1) × RTS`. VIPER takes the finite signed-Wald signature, with no internal score transformation, NES, pleiotropy correction, likelihood 1 and a minimum of 25 input targets. Taiji computes the maximum of the target/source and source/target PageRank ratios. ANANSE uses native mouse motifs, same-region ATAC for both states, mean state TPM, eligible positive target fold changes and original influence scores (`-i 500000 -j 0.05 -c prob`). The seed is 20261118 and `PYTHONHASHSEED=0`.

CompForce uses full rounded length-scaled TPM estimated counts and the official mouse ChipX network mapped to M24. It retains native lists for Welch, DESeq2, CellNet NIS, Mogrify and d’Alessio JSD and combines them by archived RRA. The differential-expression threshold is 0.05, and the Mogrify source expression threshold is 20. The sensitivity re-scales only Welch input to logCP100k. CompForce uses R seed 20261118. IMCell uses one activating seed with `in_degree` edge weights and three 1,000-draw estimates per candidate. Per-job L'Ecuyer-CMRG streams start from R seed 20261118, with one native core per seeded job and six outer workers in this command. These scientific settings and full native eligibility rules are recorded in the [Round 2 configuration](../config/single/FROZEN_CONFIG_ROUND2.yaml).

Keep this workspace free of reference and outcome files until the score lock exists. Historical QC values in the rank runner annotate the original data contract. If replayed sequencing inputs or software change preprocessing outputs, compare the fresh QC metrics to the [published QC inputs](../config/single/qc/) and report the change.

### 4. Attach evidence and evaluate

~~~sh
python "$REPO/scripts/single/stage_workspace.py" attach-references "$HAMMELMAN_WORK"
python "$HAMMELMAN_ROUND2/scripts/evaluate_round2.py"
~~~

The attachment step verifies every fresh lock artifact before copying the documented TF evidence. The evaluator verifies the lock again before exact-ID and motif-family matching. For the paper's six-workflow QC exclusion, compare the fresh ranks with the computed `TOP10_QC_EXCLUSION.tsv` from the frozen entry.

## Source and licensing notes

- Hammelman et al. (2022), [published benchmark](https://doi.org/10.1038/s41592-022-01522-2) and [ReprogrammingRecovery inputs](https://github.com/gifford-lab/ReprogrammingRecovery/tree/047cd83e90a31b527d034f3151d4ffcef551bf99).
- TRIAGE [reference RTS and upstream code](https://github.com/woojunshim/TRIAGE) at `e57492c6eb67f4ac96a7ccb121496f12d3225414`. Upstream software uses GPL-3.0. This repository contains the derived M24 mapping and a benchmark-specific runner.
- VIPER [Bioconductor source](https://bioconductor.org/packages/viper/) uses its own upstream license. The R adapter applies the documented one-column shape fix to the installed function at runtime.
- CollecTRI [Zenodo 8192729](https://zenodo.org/records/8192729), mapped to mouse with the original frozen HCOP rules. [GENCODE mouse M24](https://www.gencodegenes.org/mouse/release_M24.html) supplies the genome and annotation.
- Taiji [source and release](https://github.com/Taiji-pipeline/Taiji/releases/tag/v1.3.0) uses BSD-3-Clause. The bibliography and native resource manifests identify the CIS-BP mouse motifs. No Taiji executable or motif archive ships here.
- ANANSE, CompForce and IMCell use MIT. Their license notices accompany the bundled compatibility diffs under [config/single/licenses](../config/single/licenses/). CellNet uses [Artistic-2.0](https://spdx.org/licenses/Artistic-2.0.html). Its only bundled change is the recorded one-line import removal. All third-party executable packages are obtained separately.
