"""Rerank the frozen full score exports, then evaluate the historical TF references."""
import argparse
import csv
import gzip
import hashlib
import io
import itertools
import json
import math
from collections import defaultdict
from decimal import Decimal
from pathlib import Path
from statistics import median


DATA = Path(__file__).resolve().parents[2] / "data/single"
TARGETS = ("liver", "cardiomyocyte", "betapanc")
PRIMARY = ("TRIAGE", "VIPER", "Taiji", "ANANSE", "IMCell", "CompForce RRA")
COMPONENTS = {
    "DiffExp_ttest": "padj",
    "DiffExp_deseq": "padj",
    "CellNet": "medVal",
    "Mogrify": "sum_rank",
    "dAlessio": "specificity_score",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def truth(value):
    return str(value).lower() in {"true", "1"}


def sha(blob):
    return hashlib.sha256(blob).hexdigest()


def grouped_records(records, *fields):
    groups = defaultdict(list)
    for row in records:
        groups[tuple(row[f] for f in fields)].append(row)
    return groups


def blind_lists(name, data, lock, universe):
    compressed = (data / name).read_bytes()
    original = gzip.decompress(compressed)
    record = lock["scores"][name]
    require(sha(compressed) == record["gzip_sha256"], f"Compressed score checksum mismatch: {name}")
    require(sha(original) == record["uncompressed_sha256"], f"Original score checksum mismatch: {name}")
    rows = list(csv.DictReader(io.StringIO(original.decode("utf-8-sig")), delimiter="\t"))
    require(len(rows) == record["id_rows"], f"Full score row count mismatch: {name}")
    lists = {}
    for key, group in grouped_records(rows, "target", "method").items():
        require(len(group) == len(universe) and {r["stable_id"] for r in group} == universe,
                f"Incomplete candidate universe: {key}")
        require(all(r["config_hash"] == lock["original_config_sha256"] for r in group), f"Config mismatch: {key}")
        direction, = {r["score_direction"] for r in group}
        require(direction in {"descending", "ascending"}, f"Unexpected direction: {key}")
        eligible = [r for r in group if truth(r["eligible"])]
        require(all(Decimal(r["score"]).is_finite() for r in eligible), f"Nonfinite eligible score: {key}")
        sign = -1 if direction == "descending" else 1
        ordered = sorted(eligible, key=lambda r: (sign * Decimal(r["score"]), r["stable_id"]))
        offset = 0
        for _, ties in itertools.groupby(ordered, key=lambda r: Decimal(r["score"])):
            ties = list(ties)
            first, last = offset + 1, offset + len(ties)
            for rank, row in enumerate(ties, first):
                require((int(row["eval_rank"]), int(row["rank_min"]), int(row["rank_max"]),
                         int(row["tie_group_size"]), truth(row["tie_crosses_10"]), truth(row["top10"]))
                        == (rank, first, last, len(ties), first <= 10 < last, rank <= 10),
                        f"Frozen position or tie mismatch: {key}/{row['stable_id']}")
            offset = last
        require(all(r["eval_rank"] == "NA" and not truth(r["top10"])
                    for r in group if not truth(r["eligible"])), f"Ineligible rank: {key}")
        lists[key] = (ordered, {r["stable_id"]: r for r in group})
    return lists


def component_lists(data, lock, universe):
    lists = {}
    for target in TARGETS:
        for component, metric in COMPONENTS.items():
            prefix = f"compforce/{target}/{component}_NATIVE_"
            for kind in ("LIST", "RAW"):
                name = f"{prefix}{kind}.tsv"
                require(sha((data / name).read_bytes()) == lock["component_outputs"][name]["sha256"],
                        f"Native component file checksum mismatch: {name}")
            scored = {r["TF"]: r for r in tsv(data / f"{prefix}RAW.tsv")}
            native = sorted(tsv(data / f"{prefix}LIST.tsv"), key=lambda r: int(r["rank"]))
            require([int(r["rank"]) for r in native] == list(range(1, len(native) + 1)), f"Native ranks: {prefix}")
            require(len({r["TF"] for r in native}) == len(native), f"Duplicate native ID: {prefix}")
            require(all(r["TF"] in universe and Decimal(scored[r["TF"]][metric]).is_finite() for r in native),
                    f"Component ID or score: {prefix}")
            lists[target, component] = [{**r, "score": scored[r["TF"]][metric]} for r in native]
    return lists


def match_records(observed, original, key):
    """Test equivalent manuscript rows without depending on their physical order."""
    names = set(original[0])
    projected = [{k: str(r[k]) for k in names} for r in observed]
    order = lambda r: tuple(r[k] for k in key)
    require(sorted(projected, key=order) == sorted(original, key=order),
            f"Manuscript cross-check mismatch: {'/'.join(key)}")


def evaluate_primary(data, lists, references, universe):
    ordered_keys = [(t, m) for t in TARGETS for m in ("DE", *PRIMARY)]
    frozen = {(r["target"], r["method"]): r for r in tsv(data / "HITS_AT_10_AND_STATUS_21.tsv")}
    families = defaultdict(set)
    for row in tsv(data / "MOTIF_MAPPING_FROZEN.tsv"):
        families[row["gene_id"]].add(row["cluster"])
    summaries, positions, top_rows = [], [], []
    for target, method in ordered_keys:
        ranked, candidates = lists[target, method]
        status = frozen[target, method]
        refs = references[target]
        reference_ids = {r["gene_id"] for r in refs}
        eligible = {r["stable_id"] for r in ranked}
        top = ranked[:10]
        hit_genes = [r["symbol"] for r in top if r["stable_id"] in reference_ids]
        scoreable = {g for g, r in candidates.items() if truth(r["scoreable"])}
        tf_families = set().union(*(families[g] for g in reference_ids))
        top_families = set().union(*(families[r["stable_id"]] for r in top))
        matches = top_families & tf_families
        summary = {**status, "candidate_universe_n": len(universe), "scoreable_n": len(scoreable),
                   "eligible_n": len(eligible), "returned_n": sum(truth(r["returned"]) for r in candidates.values()),
                   "top10_n": len(top), "actual_recommendations_n": len(top),
                   "boundary_tie_n": sum(truth(r["tie_crosses_10"]) for r in ranked),
                   "reference_total_n": len(reference_ids), "reference_in_U_n": len(reference_ids & universe),
                   "reference_scoreable_n": len(reference_ids & scoreable), "reference_eligible_n": len(reference_ids & eligible),
                   "Hits_at_10": len(hit_genes), "exact_hits": ";".join(sorted(hit_genes)),
                   "exact_hits_in_rank_order": ";".join(hit_genes), "motif_family_Hits_at_10": len(matches),
                   "motif_family_reference_n": len(tf_families), "motif_family_hits": ";".join(sorted(matches)),
                   "motif_family_mapping_top10_n": sum(bool(families[r["stable_id"]]) for r in top)}
        require(all(str(v) == status[k] for k, v in summary.items()), f"Frozen summary: {target}/{method}")
        summaries.append(summary)
        for entry in refs:
            g = entry["gene_id"]
            r = candidates.get(g, {})
            hit = truth(r.get("top10"))
            reason = r.get("exclusion_reason", "")
            if g not in universe:
                label, reason = "outside_U", "Frozen reference gene is outside official candidate universe"
            elif hit:
                label, reason = "hit_top10", ""
            elif g in eligible:
                label, reason = "ranked_beyond_10", ""
            else:
                label = "scored_native_ineligible" if g in scoreable else "candidate_unscoreable"
            positions.append({"target": target, "method": method, "implementation": status["implementation"],
                              "stable_id": g, "reference_symbol": entry["gene_symbol"],
                              "reference_order": entry["reference_order"], "in_frozen_U": g in universe,
                              "candidate_status": label, "scoreable": g in scoreable, "eligible": g in eligible,
                              **{k: r.get(k, "NA") for k in ("score", "native_rank", "eval_rank", "rank_min", "rank_max")},
                              "touches_split_tie_at_10": r.get("tie_crosses_10", "False"), "hits_exact_top10": hit,
                              "execution_status": status["execution_status"],
                              "interpretation_status": status["interpretation_status"], "reason": reason})
        for row in top:
            top_rows.append({**{k: row[k] for k in ("target","method","stable_id","symbol","score","eval_rank","rank_min","rank_max")},
                             "exact_reference_hit": row["stable_id"] in reference_ids})
    return summaries, positions, top_rows


def random_expectation(n, r, k, hits):
    if n == 0:
        return 0.0, 1.0
    p = sum(math.comb(r, h) * math.comb(n - r, k - h)
            for h in range(max(hits, k - (n - r)), min(k, r) + 1)) / math.comb(n, k)
    return k * r / n, p


def table6(lists, components, references):
    rows = []
    counts = []
    for target in TARGETS:
        entries = references[target]
        reference_ids = {r["gene_id"] for r in entries}
        all_lists = {m: {r["stable_id"]: int(r["eval_rank"]) for r in lists[target, m][0]} for m in PRIMARY}
        all_lists.update({m: {r["TF"]: int(r["rank"]) for r in components[target, m]} for m in COMPONENTS})
        for entry in entries:
            rows.append({"target": target, "reference": entry["gene_symbol"], "stable_id": entry["gene_id"],
                         **{method: ranks.get(entry["gene_id"], "NA") for method, ranks in all_lists.items()}})
        for method, ranks in all_lists.items():
            r = len(reference_ids & ranks.keys())
            k = min(10, len(ranks))
            hits = sum(ranks[g] <= k for g in reference_ids & ranks.keys())
            expected, p = random_expectation(len(ranks), r, k, hits)
            counts.append({"target": target, "method": method, "eligible_TFs": len(ranks), "historical_references": len(entries),
                           "eligible_references": r, "actual_recommendations": k, "hits": hits,
                           "expected_hits": expected, "hypergeometric_upper_tail": p})
    return rows, counts


def supplementary(lists, summaries, references, symbols):
    status = {(r["target"], r["method"]): r["input_QC_status"] for r in summaries}
    union_rows, overlaps, common, qc = [], [], [], []
    for target in TARGETS:
        refs = {r["gene_id"] for r in references[target]}
        restricted = [m for m in PRIMARY if status[target, m] == "conditional"]
        eligible_sets = [set(r["stable_id"] for r in lists[target, m][0]) for m in PRIMARY]
        joint = set.intersection(*eligible_sets)
        diag = [r["gene_id"] for r in references[target] if r["gene_id"] in joint]
        common.append({"target": target, "shared_eligible_TFs": len(joint), "shared_eligible_references": len(diag),
                       "shared_reference_symbols": ",".join(symbols[g] for g in diag)})
        for budget in (5, 10, 20, 50):
            top = {m: {r["stable_id"] for r in lists[target, m][0][:budget]} for m in PRIMARY}
            union = set.union(*top.values())
            union_rows.append({"target": target, "budget_per_method": budget, "methods": len(PRIMARY),
                               "distinct_nominations": len(union), "distinct_reference_genes": len(union & refs),
                               "reference_total_n": len(refs), "exact_hit_symbols": ",".join(sorted(symbols[g] for g in union & refs))})
            def names(ids):
                return ";".join(sorted(symbols[g] for g in ids))
            for a, b in itertools.combinations(PRIMARY, 2):
                shared, combined = top[a] & top[b], top[a] | top[b]
                overlaps.append({"target": target, "budget": budget, "method_1": a, "method_2": b,
                                 "intersection_n": len(shared), "union_n": len(combined), "jaccard": len(shared) / len(combined),
                                 "shared_exact_hits": names(shared & refs), "only_1_hits": names((top[a] - top[b]) & refs),
                                 "only_2_hits": names((top[b] - top[a]) & refs)})
            if budget == 10:
                qualified = set().union(*(top[m] for m in PRIMARY if m not in restricted))
                qc.append({"target": target, "limited_configurations": len(restricted),
                           "limited_zero_hit_configurations": sum(not (top[m] & refs) for m in restricted),
                           "retained_workflows": len(PRIMARY) - len(restricted), "hits_all": len(union & refs),
                           "hits_after_QC_exclusion": len(qualified & refs), "reference_total_n": len(refs)})
    return union_rows, overlaps, common, qc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = args.data
    paper = data / "manuscript"
    lock = json.loads((data / "SCORE_LOCK.json").read_text(encoding="utf-8"))
    universe = {r["gene_id"] for r in tsv(data / "CANDIDATE_UNIVERSE.tsv")}
    require(len(universe) == 1374, "Unexpected official TF universe size")
    lists = blind_lists("FULL_SCORES_AND_RANKS.tsv.gz", data, lock, universe)
    sensitivity = blind_lists("COMPFORCE_SENSITIVITY_SCORES_AND_RANKS.tsv.gz", data, lock, universe)
    components = component_lists(data, lock, universe)
    # Historical outcomes enter only after the score hashes, full ranks, tie rules and native lists pass.
    references = {k[0]: v for k, v in grouped_records(tsv(data / "REFERENCE_TFS_FROZEN_ROUND2.tsv"), "task_id").items()}
    symbols = {g: r["symbol"] for g, r in lists["liver", "TRIAGE"][1].items()}
    symbols.update({r["gene_id"]: r["gene_symbol"] for refs in references.values() for r in refs})
    summaries, positions, tops = evaluate_primary(data, lists, references, universe)
    primary = lambda rows: [r for r in rows if r["method"] in PRIMARY]
    match_records(primary(summaries), tsv(paper / "primary_status_and_hits.tsv"), ("target", "method"))
    match_records(primary(tops), tsv(paper / "primary_top10.tsv"), ("target", "method", "eval_rank"))
    match_records(primary(positions), tsv(paper / "primary_reference_ranks.tsv"), ("target", "method", "stable_id"))
    match_records(positions, tsv(data / "REFERENCE_RANKS_LONG.tsv"), ("target", "method", "stable_id"))
    rank_table, counts = table6(lists, components, references)
    for row in tsv(paper / "compforce_reference_ranks.tsv"):
        ranked = next(r for r in rank_table if (r["target"], r["stable_id"]) == (row["target"], row["stable_id"]))
        require(str(ranked[row["component"]]) == row["native_rank"], f"Manuscript component rank mismatch: {row}")
    for row in tsv(paper / "compforce_component_counts.tsv"):
        current = next(r for r in counts if (r["target"], r["method"]) == (row["target"], row["component"]))
        require((current["eligible_TFs"], current["hits"], current["actual_recommendations"])
                == tuple(int(row[k]) for k in ("native_list_n", "Hits_at_10", "actual_recommendations_n")),
                f"Manuscript component counts mismatch: {row}")
    for row in tsv(paper / "single_conditional_random.tsv"):
        current = next(r for r in counts if (r["target"], r["method"]) == (row["target"], row["method"]))
        require(math.isclose(float(row["E"]), current["expected_hits"], rel_tol=1e-8)
                and math.isclose(float(row["p"]), current["hypergeometric_upper_tail"], rel_tol=1e-8),
                "Manuscript random baseline mismatch")
    unions, overlaps, common, qc = supplementary(lists, summaries, references, symbols)
    match_records([r for r in unions if r["budget_per_method"] != 50],
                  tsv(paper / "union_by_budget.tsv"), ("target", "budget_per_method"))
    require([r["limited_configurations"] for r in qc] == [2, 6, 2]
            and sum(r["limited_zero_hit_configurations"] for r in qc) == 2,
            "Unexpected Figure 6 input-quality status")
    require([(r["hits_all"], r["hits_after_QC_exclusion"]) for r in qc] == [(8, 6), (4, 0), (3, 3)],
            "Unexpected paper QC exclusion counts")
    require([r["shared_eligible_TFs"] for r in common] == [33, 19, 22], "Unexpected common candidate pools")
    overlap_summary = [
        {"target": t, "budget": b, "median_jaccard": median(r["jaccard"] for r in overlaps if (r["target"], r["budget"]) == (t, b))}
        for t in TARGETS for b in (10, 50)
    ]
    paper_overlap = {(r["target"], r["method_1"], r["method_2"]): r
                     for r in tsv(paper / "single_top10_overlap.tsv")}
    for row in (r for r in overlaps if r["budget"] == 10):
        expected = paper_overlap[row["target"], row["method_1"], row["method_2"]]
        require(math.isclose(row["jaccard"], float(expected["jaccard"]), rel_tol=1e-8)
                and all(str(row[k]) == expected[k] for k in
                        ("intersection_n", "union_n", "shared_exact_hits", "only_1_hits", "only_2_hits")),
                f"Manuscript shortlist overlap mismatch: {row}")
    require([round(r["median_jaccard"], 3) for r in overlap_summary] == [.053, .111, .053, .053, .000, .064],
            "Unexpected paper shortlist overlap medians")
    scale = []
    for target in TARGETS:
        key = target, "CompForce RRA"
        default, changed = [{r["stable_id"] for r in group[0][:10]} for group in (lists[key], sensitivity[key])]
        refs = {r["gene_id"] for r in references[target]}
        scale.append({"target": target, "top10_shared": len(default & changed),
                      "archived_RRA_hits": len(default & refs), "logCP100k_RRA_hits": len(changed & refs)})
    expected_scale = {r["target"]: r for r in tsv(paper / "compforce_input_scale_sensitivity.tsv")}
    for row in scale:
        expected = expected_scale[row["target"]]
        require(tuple(row[k] for k in ("top10_shared", "archived_RRA_hits", "logCP100k_RRA_hits"))
                == tuple(int(expected[k]) for k in ("top10_identity_intersection_n", "archived_Hits_at_10", "logCP100k_Hits_at_10")),
                "Manuscript CompForce input scale sensitivity mismatch")
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    for name, rows in (
        ("HITS_AT_10_AND_STATUS_21.tsv", summaries), ("REFERENCE_RANKS_LONG.tsv", positions), ("PRIMARY_TOP10.tsv", primary(tops)),
        ("TABLE6_RANKS.tsv", rank_table), ("TABLE6_COUNTS_AND_RANDOM_EXPECTATIONS.tsv", counts),
        ("TOP10_QC_EXCLUSION.tsv", qc), ("SHORTLIST_UNIONS.tsv", unions), ("SHORTLIST_OVERLAPS.tsv", overlaps),
        ("OVERLAP_SUMMARY.tsv", overlap_summary), ("COMMON_ELIGIBILITY.tsv", common), ("COMPFORCE_INPUT_SCALE.tsv", scale),
    ):
        write_tsv(out / name, rows)
    print("Single TF: both full score locks verified; Table 6, 18 primary cells, 10 QC limitations, 2 limited zero hits.")


if __name__ == "__main__":
    main()
