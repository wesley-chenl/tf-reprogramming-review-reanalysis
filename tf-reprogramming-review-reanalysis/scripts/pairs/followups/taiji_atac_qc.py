"""ATAC molecular fragment QC and whole pair filtering before replicate-level markdup."""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
import itertools
import json
from pathlib import Path
import sys

RADIUS = 2000
FAIL_TAGS = 4 | 8 | 256 | 512 | 1024 | 2048
MITO = {"chrM", "M", "MT", "chrMT"}

def tsv(path, records, cols=None):
    fields = list(records[0]) if cols is None else cols
    with path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

def filtered_pairs(path):
    counts = Counter()
    group = []
    for line in sys.stdin:
        if line.startswith("@"):
            sys.stdout.write(line)
            continue
        fields = line.rstrip("\n").split("\t", 9)
        flag = int(fields[1])
        if flag & (256 | 2048):
            counts["non_primary_alignment_records"] += 1
            continue
        if group and fields[0] != group[0][0][0]:
            screen(group, counts)
            group = []
        group.append((fields, flag, line))
    if group:
        screen(group, counts)
    tsv(path, [{"measure": key, "counts": count} for key, count in sorted(counts.items())])

def screen(group, counts):
    counts["raw_qname_groups"] += 1
    if len(group) != 2:
        counts["orphan_or_duplicate_primary_qname"] += 1
        return
    (ma, f1, a), (mb, f2, b) = group
    if not (f1 & 1 and f2 & 1 and (f1 & 64) != (f2 & 64) and (f1 & 128) != (f2 & 128)):
        counts["non_pair_flags"] += 1
        return
    if ma[2] in MITO and mb[2] in MITO and not ((f1 | f2) & 4):
        counts["mitochondrial_pairs"] += 1
        return
    if (f1 & 2) == 0 or (f2 & 2) == 0 or ((f1 | f2) & FAIL_TAGS) or ma[2] != mb[2] or int(ma[8]) == 0:
        counts["not_nuclear_proper_pair_or_failed_flag"] += 1
        return
    if min(int(ma[4]), int(mb[4])) < 30:
        counts["nuclear_proper_pair_below_mapq30_either_mate"] += 1
        return
    counts["nuclear_proper_pair_both_mapq30_pre_dedup"] += 1
    sys.stdout.write(a)
    sys.stdout.write(b)

def merged_intervals(path):
    regions = defaultdict(list)
    with path.open() as stream:
        for line in stream:
            fields = line.split()
            a, b = map(int, fields[1:3])
            if a < 0 or a >= b:
                raise ValueError(f"Invalid peak bounds: {line!r}")
            regions[fields[0]].append((a, b))
    result = {}
    for chrom, values in regions.items():
        merged = []
        for a, b in sorted(values):
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        result[chrom] = tuple(map(list, zip(*merged)))
    return regions, result

def gene_tss(path):
    points = defaultdict(set)
    with path.open() as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            col = line.split("\t", 8)
            if col[2] == "gene" and col[6] in "+-":
                pos = int(col[3] if col[6] == "+" else col[4]) - 1
                points[(col[0], col[6])].add(pos)
    return {key: sorted(values) for key, values in points.items()}

def overlaps(start, end, ranges):
    if not ranges:
        return False
    loc = bisect_left(ranges[0], end) - 1
    return loc >= 0 and ranges[1][loc] > start

def profile_insert(profile, tss, chrom, end):
    for strand in "+-":
        points = tss.get((chrom, strand), ())
        lo = bisect_left(points, end - RADIUS)
        hi = bisect_right(points, end + RADIUS)
        sgn = 1 if strand == "+" else -1
        for pos in points[lo:hi]:
            profile[(end - pos) * sgn + RADIUS] += 1

def metrics(gtf, peaks, prefix):
    raw_peaks, merged = merged_intervals(peaks)
    tss = gene_tss(gtf)
    length_histogram = Counter()
    profile = [0] * (2 * RADIUS + 1)
    frag_n = ali_n = frip_n = 0
    for line in sys.stdin:
        fields = line.split("\t", 9)
        flag, chrom, pos, tlen = int(fields[1]), fields[2], int(fields[3])-1, int(fields[8])
        if not (flag & 2) or (flag & FAIL_TAGS) or chrom in MITO or int(fields[4]) < 30:
            raise ValueError(f"Invalid retained nuclear fragment record: {fields[:9]!r}")
        ali_n += 1
        if tlen < 1:
            continue
        frag_n += 1
        length_histogram[tlen] += 1
        frip_n += overlaps(pos, pos+tlen, merged.get(chrom))
        profile_insert(profile, tss, chrom, pos+4)
        profile_insert(profile, tss, chrom, pos+tlen-5)
    if ali_n != frag_n * 2 or frag_n == 0:
        raise ValueError(f"Expected two retained primary alignments per fragment, observed {ali_n}/{frag_n}")
    flank = (sum(profile[:101]) + sum(profile[-101:])) / 202
    center = sum(profile[RADIUS-50:RADIUS+51]) / 101
    summary = {
        "state_or_run": prefix.name,
        "retained_nuclear_proper_fragments": frag_n,
        "retained_alignment_records": ali_n,
        "peaks": sum(map(len, raw_peaks.values())),
        "frip_fragments_numerator": frip_n,
        "frip_fragments_denominator": frag_n,
        "frip_fragments": frip_n/frag_n,
        "tss_central_mean_over_flank_mean": center/flank if flank else None,
        "tss_central_mean_insertions": center,
        "tss_flank_mean_insertions": flank,
        "mononucleosome_length_band_180_247bp": sum(n for length,n in length_histogram.items() if 180<=length<=247),
        "subnucleosomal_length_lt120bp": sum(n for length,n in length_histogram.items() if length<120),
        "fragments_250_450_bp": sum(n for length,n in length_histogram.items() if 250<=length<=450),
        "bam_qc_reason": "Completed valid paired fragment QC; screen-independent",
        "peak_q": 0.01,
        "peak_genome": "mm10 shared Bowtie2 GRCm38.p6 annotation resources",
    }
    tsv(prefix.with_suffix(".fragment_lengths.tsv"), [{"length_bp": bp, "fragments": n} for bp,n in sorted(length_histogram.items())])
    tsv(prefix.with_suffix(".tss_profile.tsv"), [
        {"distance_from_tss_bp": i-RADIUS, "insertions": n, "flank_normalized_insertions": n/flank if flank else ""}
        for i,n in enumerate(profile)
    ])
    tmp = prefix.with_suffix(".metrics.partial")
    tmp.write_text(json.dumps(summary, indent=2) + "\n")
    tmp.replace(prefix.with_suffix(".metrics.json"))

def summarize(receipts, dest):
    rows = []
    for f in sorted(receipts.glob("*.metrics.json")):
        rec = json.loads(f.read_text())
        label = rec["state_or_run"]
        funnel = receipts / f"{label}.peak_filter.tsv"
        markdup = receipts / f"{label}.markdup.txt"
        source = {}
        if funnel.exists():
            with funnel.open(newline="") as inp:
                source = {entry["measure"]: int(entry["counts"]) for entry in csv.DictReader(inp, delimiter="\t")}
        mark = {}
        if markdup.exists():
            for line in markdup.read_text().splitlines():
                if ": " in line:
                    key, val = line.split(": ", 1)
                    if val.isdigit():
                        mark[key] = int(val)
        rec.update(raw_qname_groups=source.get("raw_qname_groups"), mitochondrial_pairs=source.get("mitochondrial_pairs",0) if source else None,
                   pct_mitochondrial_all_qname_groups=(100*source.get("mitochondrial_pairs", 0)/source["raw_qname_groups"]) if source else None,
                   both_mates_mapq30_pre_dedup=source.get("nuclear_proper_pair_both_mapq30_pre_dedup"), duplicate_records_removed=mark.get("DUPLICATE TOTAL"),
                   duplicate_rate_eligible_pre_dedup=(mark.get("DUPLICATE TOTAL",0)/mark["READ"]) if mark.get("READ") else None)
        rows.append(rec)
    tsv(dest / "TAIJI_RUN_QC.tsv", rows)
    b1, reg1 = merged_intervals(receipts / "cortex_B1_peaks.narrowPeak")
    b2, reg2 = merged_intervals(receipts / "cortex_B2_peaks.narrowPeak")
    concordance = {
        "B1_open_intervals": sum(map(len,b1.values())),
        "B2_open_intervals": sum(map(len,b2.values())),
        "B1_peaks_with_B2_overlap": sum(overlaps(a,b,reg2.get(c)) for c,rs in b1.items() for a,b in rs),
        "B2_peaks_with_B1_overlap": sum(overlaps(a,b,reg1.get(c)) for c,rs in b2.items() for a,b in rs),
    }
    tsv(dest / "TAIJI_TARGET_BIOREPLICATE_PEAK_OVERLAP.tsv", [concordance])

def main():
    step, *args = sys.argv[1:]
    if step == "filter":
        filtered_pairs(Path(args[0]))
    elif step == "metrics":
        metrics(*map(Path, args))
    elif step == "summarize":
        summarize(*map(Path, args))
    else:
        raise ValueError(step)

if __name__ == "__main__":
    main()
