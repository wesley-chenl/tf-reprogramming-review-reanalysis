"""Generate per-task ANANSE-compatible annotation, expression and frozen regions."""
from collections import Counter, defaultdict
import gzip
from pathlib import Path
import re
import shutil

import numpy as np
import pandas as pd
from ananse.peakpredictor import _get_species, get_motif_factors
from ananse.utils import merge_regionfiles
from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs
from r2lib import ATAC_TGT, BASE, ROOT, TASKS, RNA_SRC, RNA_TGT, input_record, read_json, read_tsv, sha, write_json, write_tsv

out = ROOT / "inputs_manifest/ananse"
out.mkdir(exist_ok=True)
checked = read_json(ROOT / "inputs_manifest/VERIFIED_INPUTS.json")
fa = Path(checked["assembly"])
assert _get_species(str(fa)) == "mouse"
sizes = {v[0]: int(v[1]) for v in (l.split("\t") for l in Path(str(fa) + ".fai").read_text().splitlines())}
genes = read_tsv(BASE / "inputs/annotation/GENES_M24.tsv")
name_count = Counter(r["gene_symbol"] for r in genes)
keep = {g["gene_id"]: g for g in genes if name_count[g["gene_symbol"]] == 1 and g["chr"] in sizes and
        g["strand"] in ("+", "-") and g["gene_symbol"] and 0 <= int(g["start"]) - 1 < int(g["end"]) <= sizes[g["chr"]]}
mapping = {g: r["gene_symbol"] for g, r in keep.items()}
write_tsv(out / "GENE_SYMBOL_MAPPING.tsv", [
    {"gene_id": r["gene_id"], "symbol": r["gene_symbol"], "ANANSE_unique_GENCODE_M24_symbol": r["gene_id"] in keep,
     "reason": "" if r["gene_id"] in keep else "duplicate_or_missing_symbol_or_coordinate"}
    for r in genes])
gtf_gz = BASE / "inputs/annotation/gencode.vM24.annotation.gtf.gz"
exons = defaultdict(list)
gene_attr = re.compile(r'(?:^|;\s*)gene_id "([^"]+)"')
with gzip.open(gtf_gz, "rt") as stream:
    for line in stream:
        if line.startswith("#"):
            continue
        columns = line.rstrip("\n").split("\t")
        if columns[2] != "exon":
            continue
        stable = gene_attr.search(columns[8]).group(1).split(".", 1)[0]
        if stable in keep:
            g = keep[stable]
            assert columns[0] == g["chr"] and columns[6] == g["strand"]
            a, b = int(columns[3]) - 1, int(columns[4])
            assert int(g["start"]) - 1 <= a < b <= int(g["end"])
            exons[stable].append((a, b))
bed = out / "M24_gene_symbols.bed"
no_exon = []
with bed.open("w") as stream:
    for g in sorted(keep.values(), key=lambda r: (r["chr"], int(r["start"]), int(r["end"]), r["gene_id"])):
        start, end = int(g["start"]) - 1, int(g["end"])
        blocks = []
        for a, b in sorted(exons[g["gene_id"]]):
            if blocks and a <= blocks[-1][1]:
                blocks[-1] = (blocks[-1][0], max(blocks[-1][1], b))
            else:
                blocks.append((a, b))
        if not blocks:
            no_exon.append(g["gene_id"])
            blocks = [(start, end)]
        values = [g["chr"], start, end, g["gene_symbol"], 0, g["strand"], start, start, "0,0,0",
                  len(blocks), ",".join(str(b - a) for a, b in blocks), ",".join(str(a - start) for a, b in blocks)]
        stream.write("\t".join(map(str, values)) + "\n")
motif_config = MotifConfig()
motifs_in = Path(motif_config.get_motif_dir()) / "gimme.vertebrate.v5.0.pfm"
motif2 = motifs_in.with_suffix(".motif2factors.txt")
for p in [motifs_in, motif2]:
    shutil.copyfile(p, out / p.name)
motifs = out / motifs_in.name
tfs = defaultdict(set)
for m in read_motifs(str(motifs)):
    for factor in get_motif_factors(m, indirect=True):
        if factor[1:].islower():
            tfs[factor].add(m.id)
write_tsv(out / "MOUSE_MOTIF_TFS.tsv", [
    {"symbol": s, "number_of_native_motifs": len(tfs[s]), "M24_symbol_is_unique": name_count[s] == 1}
    for s in sorted(tfs)])
source_peaks = BASE / "inputs/atac_conditions/source/source_pool_peaks.narrowPeak"
outputs = [gtf_gz, bed, out / motif2.name, motifs, out / "MOUSE_MOTIF_TFS.tsv", out / "GENE_SYMBOL_MAPPING.tsv"]
task_summaries = []
for task in TASKS:
    task_dir = out / task
    task_dir.mkdir(exist_ok=True)
    target_peaks = BASE / "inputs/atac" / ATAC_TGT[task] / "peaks_peaks.narrowPeak"
    combined = Path(merge_regionfiles([str(source_peaks), str(target_peaks)], str(fa), str(task_dir)))
    regions = set()
    widths = Counter()
    for line in combined.read_text().splitlines():
        chrom, start, end, *_ = line.split("\t")
        start, end = int(start), int(end)
        assert 0 <= start < end <= sizes[chrom], line
        widths[end-start] += 1
        regions.add((chrom, start, end))
    assert len(regions) == sum(widths.values()) and len(regions) > 0
    folder = BASE / "inputs/rna" / task
    tpm = pd.read_csv(folder / "TPM.tsv", sep="\t", index_col=0, float_precision="round_trip")
    assert list(tpm.columns) == [*RNA_SRC, *RNA_TGT[task]]
    genes_in_rna = [g for g in tpm.index if g in mapping]
    for state, samples in [("source", list(RNA_SRC)), ("target", list(RNA_TGT[task]))]:
        expr = tpm.loc[genes_in_rna, samples].mean(axis=1).rename(index=mapping)
        assert expr.index.is_unique and np.isfinite(expr).all()
        expr.index.name = "gene"
        dest = task_dir / f"{state}_mean_TPM.tsv"
        expr.to_frame("tpm").to_csv(dest, sep="\t", float_format="%.17g")
        outputs.append(dest)
    de = pd.read_csv(folder / "DESEQ2_TARGET_V_SOURCE.tsv", sep="\t", index_col=0, float_precision="round_trip")
    de = de.loc[[g for g in de.index if g in mapping]].rename(index=mapping)
    dest = task_dir / "DESEQ2_TARGET_V_SOURCE_gene_symbols.tsv"
    de.to_csv(dest, sep="\t", na_rep="NA", float_format="%.17g")
    outputs += [dest, combined]
    task_summaries.append({"target": task, "n_regions": len(regions), "width_counts": dict(widths),
                           "source_peaks_sha256": sha(source_peaks), "target_peaks_sha256": sha(target_peaks),
                           "common_state_region_sha256": sha(combined), "n_tpm_genes": len(genes_in_rna),
                           "n_DE_genes": len(de), "same_regions_for_both_states": True})
    print(task, len(regions), "unique official union regions", flush=True)
write_json(out / "REGION_AND_MOTIF_PROVENANCE.json", {
    "assembly": str(fa), "species_detected_by_native_tool": "mouse",
    "region_rule": "fixed ANANSE 0.5.1 merge_regionfiles / CombineBedFiles width=200",
    "bed12_rule": "one row per unique exact M24 symbol; full gene span and merged annotated exons; native network uses columns 1-6",
    "gene_span_fallback_single_block_n": len(no_exon), "ANANSE_mouse_motif_TF_n": len(tfs),
    "fixed_ANANSE_provided_ATAC_motif_models": "native ATAC+motif models fitted by upstream; custom mouse coordinates, no genomic human REMAP columns",
    "tasks": task_summaries, "inputs": [input_record(p, "prepared round2 ANANSE lineage") for p in outputs]})
