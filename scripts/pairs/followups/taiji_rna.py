"""Prepare independent GSE43916 gene-level TPM and a blind expression-only comparison."""
import csv
import gzip
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

ROOT = Path(os.environ["PAIR_TAIJI_ROOT"]).resolve()
OLD = ROOT
UTILITY = ROOT
BENCH = ROOT
SAMPLES = (
    ("GSM1074345", "MEFs_rtTA_Rep1", "source"),
    ("GSM1074346", "MEFs_rtTA_Rep2", "source"),
    ("GSM1074357", "cortical_neurons", "target"),
)
META = ROOT / "metadata"
OUT = ROOT / "taiji" / "rna"
SCORE = ROOT / "scoring" / "taiji"
DESIGN = ROOT / "scoring" / "common_allowed" / "CLEAN_PAIRS.tsv"
ANNOTATION = BENCH / "inputs" / "annotation" / "GENES_M24.tsv"
SCORE_COLUMNS = ["endpoint", "variant", "pair_id", "method", "score", "score_num",
                 "score_den", "eligible", "exclusion_reason", "details"]

def load(path, delimiter="\t"):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="", encoding="utf-8-sig") as source:
        return list(csv.DictReader(source, delimiter=delimiter))

def emit(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0])
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)

def sha(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()

def site_open(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace")

def assert_sample_metadata():
    sample_data = defaultdict(list)
    sample = None
    path = META / "GSE43916_family.soft.gz"
    with site_open(path) as src:
        for line in src:
            if line.startswith("^SAMPLE = "):
                sample = line.strip().split(" = ", 1)[1]
            if line.startswith(("!Sample_title = ", "!Sample_data_processing = ")):
                sample_data[sample].append(line.strip())
    for sample, title, _ in SAMPLES:
        lines = sample_data[sample]
        assert f"!Sample_title = {title}" in lines
        assert any("RPKM" in x for x in lines)
        assert any("FPKM" in x for x in lines)
    return path

def normalize():
    metadata = assert_sample_metadata()
    matrix = META / "GSE43916_RPKM.tsv.gz"
    rows = load(matrix)
    symbols = [r["GeneSymbol"] for r in rows]
    assert len(symbols) == len(set(symbols)) and all(symbols)
    assert all(r["Transcript_id"] == r["GeneSymbol"] for r in rows)
    selected = [next(key for key in rows[0] if key.startswith(name + "_accepted_hits"))
                for _, name, _ in SAMPLES]
    per_sample = [[float(r[key]) for r in rows] for key in selected]
    assert all(math.isfinite(x) and x >= 0 for col in per_sample for x in col)
    totals = [math.fsum(col) for col in per_sample]
    assert all(total > 0 for total in totals)
    tpm = [[1e6 * x / total for x in col] for total, col in zip(totals, per_sample, strict=True)]
    annotation = load(ANNOTATION)
    native_symbol_to_id = defaultdict(set)
    for row in annotation:
        native_symbol_to_id[row["gene_symbol"]].add(row["gene_id"])
    # GeneQuant requires an annotation gene_name in column 1 and unlogged TPM in column 2.
    # Keep unmatched RPKM rows in the TPM denominator and in the transparent RNA export.
    table = {}
    for symbol, sg1, sg2, tg in zip(symbols, *tpm, strict=True):
        table[symbol] = {
            "gene_symbol": symbol, "annotation_match_count": len(native_symbol_to_id[symbol]),
            **{spec[0]: value for spec, value in zip(SAMPLES, (sg1, sg2, tg), strict=True)},
            "source_mean_tpm": (sg1 + sg2) / 2, "target_tpm": tg,
            "target_minus_source_log2p1_tpm": math.log2(1 + tg) - (math.log2(1 + sg1) + math.log2(1 + sg2))/2,
        }
    emit(OUT / "FULL_GENE_TPM.tsv", list(table.values()))
    manifest = []
    for i, (sample, title, role) in enumerate(SAMPLES):
        destination = OUT / f"{sample}.tsv"
        with destination.open("w", encoding="utf-8", newline="") as file:
            for symbol, row in sorted(table.items()):
                if row["annotation_match_count"] == 1:
                    file.write(f"{symbol}\t{row[sample]:.17g}\n")
        manifest.append({
            "GEO": sample, "study": "GSE43916", "title": title, "rna_condition": role,
            "original_column": selected[i], "original_unit": "linear RefSeq RPKM (SOFT also says FPKM)",
            "original_gene_rows": len(symbols), "original_rpkm_sum": totals[i],
            "converted_full_tpm_sum": math.fsum(tpm[i]),
            "genequant_export_rows": sum(r["annotation_match_count"] == 1 for r in table.values()),
            "export_matched_tpm_sum": math.fsum(r[sample] for r in table.values() if r["annotation_match_count"] == 1),
            "full_tpm_denominator_before_m24_matching": True, "sha256_original": sha(matrix),
            "sha256_geo_soft": sha(metadata), "sha256_annotation_table": sha(ANNOTATION),
            "genequant_path": destination.name, "genequant_sha256": sha(destination),
        })
    emit(OUT / "TAIJI_RNA_INPUTS.tsv", manifest)
    return table

def old_percentiles(node_effect):
    ordered = sorted((effect, gene) for gene, effect in node_effect.items())
    i = 0
    percent = {}
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][0] == ordered[i][0]:
            j += 1
        for _, member in ordered[i:j]:
            percent[member] = Fraction(i+j, 2*len(ordered))
        i = j
    return percent

def measured_probe_symbol(symbol, geneid, names, secondary):
    if geneid in {"", "NA", "unresolved"}:
        return None
    if symbol in names and names[symbol] == geneid:
        return symbol
    aliases = [alias for alias in secondary[geneid] if alias in names and names[alias] == geneid]
    return aliases[0] if len(aliases) == 1 else None

def expression_control(table, design):
    # Only a current name or an unambiguous old NCBI synonym with a common Entrez ID.
    info = load(UTILITY / "expression_raw" / "Mus_musculus.gene_info.gz")
    evidence, aliases, canonical = defaultdict(set), defaultdict(set), defaultdict(set)
    for row in info:
        geneid = "GeneID:" + row["GeneID"]
        canonical[row["Symbol"]].add(geneid)
        evidence[row["Symbol"]].add(geneid)
        for alias in row["Synonyms"].split("|"):
            if alias in table:
                evidence[alias].add(geneid)
                aliases[geneid].add(alias)
    source_names = {
        symbol: next(iter(canonical[symbol] or ids)) for symbol, ids in evidence.items()
        if symbol in table and len(canonical[symbol] or ids) == 1
    }
    members = sorted({(p[f"{s}_symbol"], p[f"{s}_id"]) for p in design for s in ("bhlh", "partner")})
    audit = []
    mapped = {}
    for symbol, geneid in members:
        found = measured_probe_symbol(symbol, geneid, source_names, aliases)
        mapped[(symbol, geneid)] = found
        audit.append({"candidate_symbol": symbol, "candidate_id": geneid, "rna_gene_symbol": found or "",
                      "resolution": ("exact_symbol_entrez" if found == symbol else "unique_entrez_alias") if found else "unmapped_or_ambiguous"})
    emit(OUT / "TAIJI_RNA_CANDIDATE_MAP.tsv", audit)
    by_variant = defaultdict(list)
    for p in design:
        by_variant[p["variant"]].append(p)
    scores = []
    for variant, pairs in sorted(by_variant.items()):
        ranks = {}
        for role in ("bhlh", "partner"):
            ids = {(p[f"{role}_symbol"], p[f"{role}_id"]) for p in pairs}
            ranks[role] = old_percentiles({key: table[mapped[key]]["target_minus_source_log2p1_tpm"]
                                          for key in ids if mapped[key] is not None})
        for p in pairs:
            identifiers = [(p[f"{role}_symbol"], p[f"{role}_id"]) for role in ("bhlh", "partner")]
            bhlh, partner = (ranks[role].get(key) for role, key in zip(("bhlh", "partner"), identifiers, strict=True))
            e = (bhlh+partner)/2 if bhlh is not None and partner is not None else None
            scores.append({
                "endpoint": "GSE43916", "variant": variant, "pair_id": p["pair_id"],
                "method": "expression_GSE43916", "score": "" if e is None else format(float(e), ".17g"),
                "score_num": "" if e is None else e.numerator, "score_den": "" if e is None else e.denominator,
                "eligible": int(e is not None),
                "exclusion_reason": "ambiguous_or_unmeasured_member" if e is None else "",
                "details": json.dumps({"GSE43916_bhlh": mapped[identifiers[0]], "GSE43916_partner": mapped[identifiers[1]],
                    "bhlh_role_percentile": str(bhlh), "partner_role_percentile": str(partner)}, separators=(",", ":")),
            })
    assert len(scores) == 9 * 598
    emit(SCORE / "TAIJI_EXPRESSION_CONTROL.tsv", scores, SCORE_COLUMNS)
    print("GSE43916 expression coverage", Counter((r["variant"], r["eligible"]) for r in scores))

def main():
    assert (ROOT / "ADDENDUM_PROTOCOL.yaml").exists()
    OUT.mkdir(parents=True, exist_ok=True)
    SCORE.mkdir(parents=True, exist_ok=True)
    table = normalize()
    expression_control(table, load(DESIGN))

if __name__ == "__main__":
    main()
