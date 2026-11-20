#!/usr/bin/env python3
"""Prepare only the independent, prespecified mouse expression endpoints."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import os
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
from statistics import median
import subprocess
import tarfile

ROOT = Path(os.environ["PAIR_STAGE"]).resolve()
RAW = ROOT / "expression_raw"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"
PUBLIC = {
    "GSE53951_RAW.tar": (f"{GEO}/series/GSE53nnn/GSE53951/suppl/GSE53951_RAW.tar", "GSE53951"),
    "GSE53951_family.soft.gz": (f"{GEO}/series/GSE53nnn/GSE53951/soft/GSE53951_family.soft.gz", "GSE53951"),
    "GSE114809_fpkm.samples.log2.txt.gz": (f"{GEO}/series/GSE114nnn/GSE114809/suppl/GSE114809_fpkm.samples.log2.txt.gz", "GSE114809"),
    "GSE114809_family.soft.gz": (f"{GEO}/series/GSE114nnn/GSE114809/soft/GSE114809_family.soft.gz", "GSE114809"),
    "Mus_musculus.gene_info.gz": ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz", "NCBI Gene 10090"),
    "gene_history.gz": ("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz", "NCBI Gene 10090"),
}
SAMPLES = {
    "GSE53951": [
        ("GSM1303947", "MEFs_mock_sample1", "source"),
        ("GSM1303949", "MEFs_mock_sample2", "source"),
        ("GSM1303943", "neurons_mock_sample1", "target"),
        ("GSM1303945", "neurons_mock_sample2", "target"),
    ],
    "GSE114809": [
        ("GSM3150587", "PMEF-1", "source"),
        ("GSM3150588", "PMEF-2", "source"),
        ("GSM3150589", "PMEF-3", "source"),
        ("GSM3150590", "Pneuron-1", "target"),
        ("GSM3150591", "Pneuron-2", "target"),
        ("GSM3150592", "Pneuron-3", "target"),
    ],
}
EFFECT_FIELDS = ("stable_id", "symbol", "d_g", "mean_source", "mean_target", "measurement_rows")

def write_tsv(path, fields, rows):
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as stream:
        out = csv.writer(stream, delimiter="\t", lineterminator="\n")
        out.writerow(fields)
        out.writerows(rows)
    partial.replace(path)

def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def header_field(headers, key):
    entries = re.findall(rf"^{key}:\s*(.*?)\s*$", headers, re.I | re.M)
    return entries[-1] if entries else ""

def verify_file(path, size=""):
    if size and path.stat().st_size != int(size):
        raise ValueError(f"Unexpected Content-Length for {path.name}")
    if path.name.endswith(".gz"):
        with path.open("rb") as fh:
            if fh.read(2) != b"\x1f\x8b":
                raise ValueError(f"Not gzip: {path.name}")
    elif path.suffix == ".tar":
        with tarfile.open(path, "r:") as tar:
            if len(tar.getmembers()) != 10:
                raise ValueError("Expected 10 CEL archive members")

def obtain(entry):
    name, (url, series) = entry
    path = RAW / name
    info = RAW / f"{name}.headers"
    if not path.is_file() or not path.stat().st_size:
        partial = RAW / f"{name}.partial"
        subprocess.run(["curl", "-fsSL", "--retry", "3", "--connect-timeout", "15", "--max-time", "200",
                        "-D", str(info), url, "-o", str(partial)], check=True)
        partial.replace(path)
    headers = info.read_text(encoding="latin1") if info.exists() else ""
    verify_file(path, header_field(headers, "Content-Length"))
    return (str(path.relative_to(ROOT)), url, header_field(headers, "Date") or datetime.fromtimestamp(
        path.stat().st_mtime, timezone.utc).isoformat(), path.stat().st_size, digest(path), series,
        header_field(headers, "Last-Modified"), "verified_public_download", "public GEO/NCBI access")

def soft_samples(series):
    result = {}
    current = None
    with gzip.open(RAW / f"{series}_family.soft.gz", "rt", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if line.startswith("^"):
                kind, _, value = line.strip().partition(" = ")
                current = None
                if kind == "^SAMPLE":
                    current = defaultdict(list)
                    if value in result:
                        raise ValueError(f"Duplicate SOFT sample {value}")
                    result[value] = current
            elif current is not None and line.startswith("!Sample_"):
                key, _, value = line.strip().partition(" = ")
                current[key].append(value)
    return result

def extract_cels():
    target = RAW / "GSE53951_cels"
    target.mkdir(exist_ok=True)
    whitelist = {gsm for gsm, _, _ in SAMPLES["GSE53951"]}
    result = {}
    with tarfile.open(RAW / "GSE53951_RAW.tar", "r:") as tar:
        for member in tar:
            name = Path(member.name).name
            gsm = name.split("_", 1)[0]
            if gsm not in whitelist:
                continue
            assert member.isfile() and gsm not in result and name.endswith(".CEL.gz"), name
            destination = target / (gsm + ".CEL")
            zipped = tar.extractfile(member)
            assert zipped is not None
            with gzip.GzipFile(fileobj=zipped) as cel, destination.open("wb") as output:
                shutil.copyfileobj(cel, output)
            with destination.open("rb") as stream:
                magic = stream.read(5)
            assert magic == b"\x40\x00\x00\x00\x04" or magic.startswith((b";\x01", b"[CEL]")), (name, magic.hex())
            result[gsm] = (name, digest(destination), magic.hex(), str(destination.relative_to(ROOT)))
    assert set(result) == whitelist and len({v[1] for v in result.values()}) == 4
    return result

def manifest_and_samples():
    RAW.mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=min(6, len(PUBLIC))) as pool:
        meta = list(pool.map(obtain, PUBLIC.items()))
    write_tsv(ROOT / "SOURCE_EXPRESSION.tsv",
              ("path", "url", "downloaded_utc", "bytes", "sha256", "accession", "version_last_modified",
               "status", "licence_or_access"), sorted(meta))
    cel = extract_cels()
    with gzip.open(RAW / "GSE114809_fpkm.samples.log2.txt.gz", "rt", encoding="utf-8") as fh:
        head = fh.readline().rstrip().split("\t")
    rows = []
    for series, whitelist_rows in SAMPLES.items():
        table = soft_samples(series)
        whitelist = {gsm: (title, role) for gsm, title, role in whitelist_rows}
        assert set(whitelist) <= set(table)
        if series == "GSE53951":
            assert len(table) == 10
        else:
            assert len(table) == 12
        for gsm, sample in table.items():
            title, = sample["!Sample_title"]
            species, = sample["!Sample_organism_ch1"]
            platform, = sample["!Sample_platform_id"]
            source, = sample["!Sample_source_name_ch1"]
            characteristics = sample["!Sample_characteristics_ch1"]
            assert species == "Mus musculus" and platform == ("GPL1261" if series == "GSE53951" else "GPL17021")
            cel_rows = cel.get(gsm, ("", "", "", ""))
            role = whitelist.get(gsm, ("", "exclude"))[1]
            notes = "donor matching and technical replication not asserted"
            expression_col = ""
            if role != "exclude":
                assert title == whitelist[gsm][0]
                if series == "GSE53951":
                    assert "treated with: mock" in characteristics and "E14.5" in " ".join(characteristics)
                    assert ("fibroblast" if role == "source" else "cortical/hippocampal") in source
                    assert next(iter(sample["!Sample_supplementary_file"])).endswith(cel_rows[0])
                else:
                    assert ("embryonic fibroblast" if role == "source" else "primary cortical neuron") in source
                    expression_col, = [label for label in head if label.lower() == title.lower()]
                    if expression_col != title:
                        notes += "; case-only title to submitted column"
            else:
                notes = "IFN-treated or L929" if series == "GSE53951" else "induced neuron; neither read nor transformed for scoring"
            rows.append([series, gsm, role, title, species, platform, source, " | ".join(characteristics), expression_col,
                         *cel_rows, f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}", notes])
    write_tsv(ROOT / "validation" / "EXPRESSION_SAMPLE_MANIFEST.tsv",
              ("series", "gsm", "use", "title", "organism", "platform", "source_name", "characteristics",
               "author_expression_column", "raw_tar_member", "uncompressed_cel_sha256", "cel_magic_hex",
               "cel_local_path", "sample_url", "replication_status"), rows)

def gene_crosswalk():
    official = {}
    alias = defaultdict(set)
    with gzip.open(RAW / "Mus_musculus.gene_info.gz", "rt", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["#tax_id"] != "10090":
                continue
            gid, symbol = row["GeneID"], row["Symbol"]
            assert gid not in official
            official[gid] = symbol
            if row["Synonyms"] not in {"", "-"}:
                for alt in row["Synonyms"].split("|"):
                    alias[alt].add(gid)
    by_symbol = defaultdict(set)
    for gid, symbol in official.items():
        by_symbol[symbol].add(gid)
    return official, by_symbol, alias

def retired_ids():
    with gzip.open(RAW / "gene_history.gz", "rt", encoding="utf-8") as stream:
        return {row["Discontinued_GeneID"]: row["GeneID"] for row in csv.DictReader(stream, delimiter="\t") if row["#tax_id"] == "10090"}

def unique_current_id(current, retired, gid):
    seen = set()
    while gid not in current and gid in retired and gid not in seen:
        seen.add(gid)
        gid = retired[gid]
    return gid if gid in current else None

def rna_effects(official, by_symbol, alias):
    rows_by_gid = defaultdict(list)
    audit = []
    stats = Counter()
    source = RAW / "GSE114809_fpkm.samples.log2.txt.gz"
    prior = {"source": [x[1].lower() for x in SAMPLES["GSE114809"] if x[2] == "source"],
             "target": [x[1].lower() for x in SAMPLES["GSE114809"] if x[2] == "target"]}
    with gzip.open(source, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        assert len(header) == 13 and header[0] == "test_id" and len(set(h.lower() for h in header)) == 13
        ix = {k: [header.index(next(h for h in header if h.lower() == n)) for n in prior[k]] for k in prior}
        assert len(set(ix["source"] + ix["target"])) == 6
        for line, row in enumerate(reader, 2):
            stats["input_rows"] += 1
            if len(row) != len(header):
                raise ValueError(f"RNA sheet line {line} has unexpected number of columns")
            name = row[0].strip()
            group = by_symbol.get(name)
            method = "official_exact_symbol"
            if not group:
                group = alias.get(name)
                method = "unique_ncbi_synonym" if group else "unmapped_submitted_id"
            if group and len(group) != 1:
                method = "ambiguous_submitted_id"
            gene = next(iter(group)) if group and len(group) == 1 else ""
            vals = None
            try:
                vals = {k: [float(row[i]) for i in ix[k]] for k in ix}
            except ValueError:
                pass
            if vals is None or not all(math.isfinite(v) for group_values in vals.values() for v in group_values):
                method = "non_finite"
            stats[method] += 1
            if gene and method not in {"non_finite", "ambiguous_submitted_id"}:
                rows_by_gid[gene].append((line, name, method, vals))
            else:
                audit.append((line, name, gene, method, "excluded"))
    effects = []
    finite_min, finite_max = math.inf, -math.inf
    for gene, measurements in rows_by_gid.items():
        exact = [row for row in measurements if row[2] == "official_exact_symbol"]
        use = (exact or measurements)[0] if len(exact or measurements) == 1 else None
        for item in measurements:
            decision = "retained" if item is use else ("prefer_single_official_record" if use else "duplicate_gene_ambiguous")
            audit.append((*item[:2], gene, item[2], decision))
            stats[decision] += 1
        if use is None:
            continue
        vals = use[-1]
        finite_min = min(finite_min, *vals["source"], *vals["target"])
        finite_max = max(finite_max, *vals["source"], *vals["target"])
        avg_s = math.fsum(vals["source"]) / 3
        avg_t = math.fsum(vals["target"]) / 3
        effects.append((gene, official[gene], avg_t - avg_s, avg_s, avg_t, 1))
    write_tsv(ROOT / "validation" / "GSE114809_GENE_ID_AUDIT.tsv",
              ("source_line", "submitted_test_id", "stable_id", "mapping_basis", "decision"), sorted(audit))
    stats["retained_measured_zero_effect"] = sum(row[2] == 0 for row in effects)
    write_json(ROOT / "validation" / "GSE114809_EXPRESSION_QC.json", {
        "sample_columns": {k: [header[i] for i in ix[k]] for k in ix},
        "raw_value_units": "files described as log2 FPKM; pseudocount not independently specified in GEO SOFT",
        "submitted_range_retained_selected_values": [finite_min, finite_max],
        "mapping_and_qc_counts": dict(stats), "score_genes": len(effects),
    })
    if len(effects) < 1000:
        raise ValueError("Unexpected gene coverage for GSE114809")
    write_tsv(ROOT / "inputs_scoring" / "GSE114809_EFFECT.tsv", EFFECT_FIELDS,
              sorted(effects, key=lambda row: int(row[0])))
    print(f"GSE114809: {len(effects)} stable mouse genes", flush=True)

def microarray_effects(official, retired):
    rma = RAW / "GPL1261_RMA_PROBES.tsv"
    grouping = defaultdict(list)
    audit = []
    counts = Counter()
    cols = [s[0] for s in SAMPLES["GSE53951"]]
    with rma.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        assert reader.fieldnames == ["probeset_id", "bioc_entrez_ids", *cols], reader.fieldnames
        for row in reader:
            probe = row["probeset_id"]
            ids = row["bioc_entrez_ids"].split("|") if row["bioc_entrez_ids"] else []
            reason = "included"
            stable = ""
            if probe.upper().startswith("AFFX"):
                reason = "control_probe"
            elif len(ids) > 1:
                reason = "bioconductor_multiple_gene_ids"
            elif not ids:
                reason = "bioconductor_no_gene_id"
            else:
                stable = unique_current_id(official, retired, ids[0]) or ""
                if not stable:
                    reason = "retired_or_not_current_ncbi_mouse_id"
            if stable:
                vals = tuple(float(row[c]) for c in cols)
                assert all(math.isfinite(v) for v in vals)
                grouping[stable].append(vals)
            counts[reason] += 1
            audit.append((probe, row["bioc_entrez_ids"], stable, official.get(stable, ""), reason))
    effects = []
    for gene, values in grouping.items():
        medians = tuple(median(v[i] for v in values) for i in range(4))
        av_s = (medians[0] + medians[1]) / 2
        av_t = (medians[2] + medians[3]) / 2
        effects.append((gene, official[gene], av_t - av_s, av_s, av_t, len(values)))
    stats = {"probeset_mapping_counts": dict(counts),
             "score_genes": len(effects),
             "same_stable_id_multiple_probe_groups": sum(r[-1] > 1 for r in effects),
             "retained_measured_zero_effect": sum(r[2] == 0 for r in effects)}
    write_tsv(ROOT / "validation" / "GPL1261_PROBE_MAPPING_AUDIT.tsv",
              ("probeset_id", "bioc_entrez_ids", "current_stable_id", "current_ncbi_symbol", "decision"), sorted(audit))
    write_json(ROOT / "validation" / "GSE53951_EXPRESSION_QC.json", stats)
    if len(effects) < 1000:
        raise ValueError("Unexpected gene coverage for GSE53951")
    write_tsv(ROOT / "inputs_scoring" / "GSE53951_EFFECT.tsv", EFFECT_FIELDS,
              sorted(effects, key=lambda row: int(row[0])))
    print(f"GSE53951: {len(effects)} stable mouse genes", flush=True)

def main():
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("stage", choices=["raw", "rna", "microarray"])
    args = command.parse_args()
    if args.stage == "raw":
        manifest_and_samples()
    else:
        genes, symbols, aliases = gene_crosswalk()
        if args.stage == "rna":
            rna_effects(genes, symbols, aliases)
        else:
            microarray_effects(genes, retired_ids())

if __name__ == "__main__":
    main()
