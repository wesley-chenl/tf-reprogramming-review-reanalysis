"""Rebuild all candidate ranks from native outputs, with no evaluation labels in scope."""
import argparse
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

from r2lib import BASE, METHODS, ROOT, TASKS, input_record, read_json, read_tsv, sg, sha, write_json, write_tsv


RANK_FIELDS = ("target", "method", "implementation", "stable_id", "symbol", "score", "score_direction", "native_rank",
               "eval_rank", "rank_min", "rank_max", "tie_group_size", "tie_crosses_10", "scoreable", "returned", "eligible",
               "exclusion_reason", "top10", "input_QC_status", "QC_limitations", "prior_version", "config_hash")
AVAILABLE = {"completed", "completed_empty", "completed_uninformative", "reused_verified", "exported_verified"}
PRIMARY = "CompForce_RRA5_archived"
SENSITIVITY = "CompForce_RRA5_ttest_logCP100k_sensitivity"
IMPLEMENTATIONS = {"DE": "frozen_external_DESeq2_signed_Wald", "TRIAGE": "original_discordance_score",
                   "VIPER": "VIPER_CollecTRI_mouse_adapter", "Taiji": "Taiji_1.3.0_bilateral_PageRank_ratio",
                   "ANANSE": "ANANSE_v0.5.1_influence", "CompForce RRA": PRIMARY,
                   "IMCell": "IMCell_v0.1.0_signed_singleton_3x1000"}
PRIOR = {"DE": "no prior; GENCODE M24 and frozen DESeq2 contrast", "TRIAGE": "frozen published RTS mapped to M24",
         "VIPER": "frozen mouse CollecTRI (Zenodo 8192729, HCOP)", "Taiji": "frozen CIS-BP mouse motif resources, Taiji 1.3.0",
         "ANANSE": "Gimme vertebrate v5.0 mouse motifs; archived ATAC+motif models",
         "CompForce RRA": "v0.1.0 official mouse_chipX; mapped to M24",
         "IMCell": "frozen mouse CollecTRI signed full network, unit edge weights"}


def is_true(value):
    return str(value).lower() in ("true", "1", "t")


def number(value):
    try:
        x = float(value)
    except (ValueError, TypeError):
        return None
    return x if math.isfinite(x) else None


def id_rows(path, field, *, optional=False):
    if optional and not path.is_file():
        return {}
    rows = read_tsv(path)
    result = {r[field]: r for r in rows}
    assert len(result) == len(rows), (path, field)
    return result


def dscore(value):
    return Decimal(str(value)) if number(value) is not None else None


def policy_quality(task, method):
    flags = []
    if method in ("Taiji", "ANANSE"):
        flags.append("common source pooled ATAC: 21041 peaks; 3.418% alignment records / 3.873% TLEN fragments")
        if task == "cardiomyocyte":
            flags.append("cardiomyocyte ATAC: single-end; 113 peaks; 0.030% alignment records")
    if task == "cardiomyocyte":
        flags.append("cardiomyocyte RNA: Salmon 28.615% and 27.912% mapped")
    return "conditional" if flags else "valid", " | ".join(flags)


def score_order(rows, direction):
    descending = direction == "descending"
    values = [(dscore(r["score"]), r["stable_id"], r) for r in rows if is_true(r["eligible"]) and number(r["score"]) is not None]
    values.sort(key=lambda e: (-e[0] if descending else e[0], e[1]))
    return values


def ranks(rows, direction, usability=True):
    order = score_order(rows, direction)
    position = 0
    while position < len(order):
        start = position
        position += 1
        while position < len(order) and order[position][0] == order[start][0]:
            position += 1
        for index in range(start, position):
            r = order[index][2]
            r.update(eval_rank=index + 1, rank_min=start + 1, rank_max=position, tie_group_size=position - start,
                     tie_crosses_10=start + 1 <= 10 < position, top10=usability and index < 10)
    return [r for _, _, r in order]


def attempt(data, key):
    if not data.is_file():
        return {}, "NA"
    return read_json(data), str(data)


def observation(score=None, reason="", *, native="NA", returned=None, eligible=None, scoreable=None, **extras):
    finite = number(score) is not None
    return {"score": score if finite else "NA", "native_rank": native, "returned": finite if returned is None else returned,
            "scoreable": finite if scoreable is None else scoreable, "eligible": finite if eligible is None else eligible,
            "exclusion_reason": reason, **extras}


def old_results(task, method, universe, lineage, comparisons):
    row = lineage[task, method]
    path = Path(row["original_score_file"])
    assert row["status"] == "reuse_verified" and sha(path) == row["original_score_sha256"]
    original = id_rows(path, "gene_id")
    assert set(original) == set(universe)
    items = {}
    for stable, r in original.items():
        assert r["task_id"] == task and r["config_hash"] == row["previous_config_sha256"]
        assert r["gene_symbol"] == universe[stable]["gene_symbol"], (task, method, stable)
        score = r["ranking_score"]
        denominator = number(r.get("PR_source"))
        numerator = number(r.get("PR_target"))
        if method == "Taiji" and is_true(r["eligible_for_recommendation"]):
            assert denominator and numerator and math.isclose(float(score), max(numerator / denominator, denominator / numerator), rel_tol=1e-12)
        finite = number(score) is not None
        if is_true(r["eligible_for_recommendation"]):
            assert finite
        items[stable] = observation(score, r["missing_reason"] or ("" if finite else r["candidate_status"]),
                   native=r["rank"], scoreable=is_true(r["scorable"]), returned=finite,
                   eligible=is_true(r["eligible_for_recommendation"]), original_candidate_status=r["candidate_status"],
                   original_raw_score=r["raw_score"], original_implementation=r["implementation_id"],
                   original_pre_prediction_hash=r["input_hash"], original_score_file_sha256=row["original_score_sha256"],
                   VIPER_prior_targets=r.get("n_prior_targets", "NA"), VIPER_input_targets=r.get("n_input_targets", "NA"),
                   VIPER_minsize_pass=r.get("minsize_pass", "NA"), Taiji_PR_source=r.get("PR_source", "NA"),
                   Taiji_PR_target=r.get("PR_target", "NA"), Taiji_signed_log_ratio=r.get("signed_log_ratio", "NA"))
    decisions = sorted(((int(r["rank"]), g) for g, r in original.items() if is_true(r["in_top10"])))
    comparisons[task, method] = {"old_top10": [g for _, g in decisions], "old_scores_verified": len(original),
                                 "old_ranks": {g: r["rank"] for g, r in original.items()}, "artifact_sha256": sha(path)}
    return items, "reused_verified", "", str(path)


def de_results(task, universe):
    path = BASE / "inputs/rna" / task / "DESEQ2_TARGET_V_SOURCE.tsv"
    de = id_rows(path, "gene_id")
    items = {}
    for stable, u in universe.items():
        current = de.get(stable, {})
        stat = current.get("stat")
        reason = ""
        if number(stat) is None:
            reason = "candidate_ID_unresolved" if u["mapping_status"] != "mapped_M24" else "no_finite_Wald_statistic"
        else:
            fc = number(current["log2FoldChange"])
            assert fc is not None and float(stat) * fc >= 0
        items[stable] = observation(stat, reason, returned=bool(current), DE_padj=current.get("padj", "NA"),
                                    DE_log2FoldChange=current.get("log2FoldChange", "NA"), DE_baseMean=current.get("baseMean", "NA"),
                                    DE_pvalue=current.get("pvalue", "NA"))
    return items, "exported_verified", "", str(path)


def ananse_results(task, universe):
    path = ROOT / "results" / task / "ANANSE"
    influence, at = attempt(path / "influence_RUN.json", "ANANSE")
    status = influence.get("execution_status", "blocked_dependency")
    native = id_rows(path / "influence.tsv", "factor") if status == "completed" else {}
    native_order = {tf: i+1 for i, tf in enumerate(native)}
    motifs = id_rows(ROOT / "inputs_manifest/ananse/MOUSE_MOTIF_TFS.tsv", "symbol")
    gene_map = id_rows(ROOT / "inputs_manifest/ananse/GENE_SYMBOL_MAPPING.tsv", "gene_id")
    binding = {s: (read_json(path / f"binding_{s}_RUN.json") if (path / f"binding_{s}_RUN.json").is_file() else {})
               for s in ("source", "target")}
    tf_binding = {s: set(r.get("finite_binding_factors", [])) for s, r in binding.items()}
    state_edges = {s: id_rows(path / f"{s}_NETWORK_OUTDEGREE.tsv", "symbol", optional=True) for s in ("source", "target")}
    gate = id_rows(path / "INFLUENCE_NATIVE_GATE.tsv", "symbol", optional=True)
    de = id_rows(ROOT / "inputs_manifest/ananse" / task / "DESEQ2_TARGET_V_SOURCE_gene_symbols.tsv", "gene_id")
    items = {}
    for stable, u in universe.items():
        sym = u["gene_symbol"]
        uniquely_annotated = is_true(gene_map.get(stable, {}).get("ANANSE_unique_GENCODE_M24_symbol"))
        covered_motif = sym in motifs
        own_de = de.get(sym, {})
        g = gate.get(sym, {})
        source_out = state_edges["source"].get(sym, {}).get("out_degree", 0 if state_edges["source"] else "NA")
        target_out = state_edges["target"].get(sym, {}).get("out_degree", 0 if state_edges["target"] else "NA")
        diff_out = int(g.get("differential_network_out_degree", 0)) if gate else "NA"
        result = native.get(sym) if uniquely_annotated else None
        if result is not None:
            assert int(diff_out) > 0 and is_true(g["native_DE_present"]) and float(g["native_real_fold_change"]) > 0
            reason = ""
        elif status != "completed":
            reason = f"task_{status}"
        elif not uniquely_annotated:
            reason = "unresolved_or_nonunique_M24_symbol"
        elif not covered_motif:
            reason = "no_native_mouse_motif_factor"
        elif sym not in tf_binding["target"]:
            reason = "native_target_binding_not_returned"
        elif int(target_out) == 0:
            reason = "no_target_network_out_edges"
        elif diff_out == 0:
            reason = "no_out_edges_after_native_differential_graph_gate"
        elif not is_true(g.get("native_DE_present")):
            reason = "native_DE_filtered_or_absent_requires_finite_padj_and_fold_change"
        elif float(g["native_real_fold_change"]) <= 0:
            reason = "native_own_log2FoldChange_not_positive"
        else:
            reason = "native_influence_not_returned"
        items[stable] = observation((result or {}).get("influence_score"), reason,
                                    native=native_order.get(sym, "NA"), returned=result is not None,
                                    ANANSE_unique_M24_annotation=uniquely_annotated, ANANSE_mouse_motif=covered_motif,
                                    ANANSE_binding_source=(sym in tf_binding["source"] if tf_binding["source"] else "NA"),
                                    ANANSE_binding_target=(sym in tf_binding["target"] if tf_binding["target"] else "NA"),
                                    ANANSE_source_out_degree=source_out, ANANSE_target_out_degree=target_out,
                                    ANANSE_differential_out_degree=diff_out, ANANSE_native_DE_present=g.get("native_DE_present", "NA"),
                                    ANANSE_native_own_log2FoldChange=g.get("native_real_fold_change", "NA"),
                                    DE_padj=own_de.get("padj", "NA"), DE_log2FoldChange=own_de.get("log2FoldChange", "NA"))
    reason = influence.get("reason", "")
    if not reason and status != "completed":
        markers = [(p,read_json(p)) for p in sorted(path.glob("*_RUN.json")) if p.name != "influence_RUN.json"]
        failed = [f"{p.name}: {r.get('reason',r['execution_status'])}" for p,r in markers if r.get("execution_status") not in ("completed",)]
        reason = " | ".join(failed) or "Influence stage did not produce a valid run record"
    return items, status, reason, at


def compforce_results(task, universe, variant=False):
    path = ROOT / "results" / task / "CompForce"
    run, at = attempt(path / "METHOD_RUN.json", "CompForce")
    result = run.get("sensitivity" if variant else "primary", {})
    status = result.get("execution_status", run.get("execution_status", "failed_execution"))
    native = {}
    if "output" in result:
        agg = Path(result["output"]["path"])
        assert sha(agg) == result["output"]["sha256"]
        native = id_rows(agg, "TF")
    cov = id_rows(path / "COMPFORCE_CANDIDATE_INPUT.tsv", "TF", optional=True)
    nit = id_rows(path / "CellNet_TARGET_SUBNETWORK_TFS.tsv", "TF", optional=True)
    before = id_rows(path / "Mogrify_BEFORE_NATIVE_FILTERS.tsv", "gene", optional=True)
    redundant = set(id_rows(path / "Mogrify_REDUNDANT_TFS.tsv", "TF", optional=True))
    lists = {}
    for label, meta in result.get("inputs", {}).items():
        component = next((r for r in meta.get("file_sha256", []) if r["path"].endswith("_NATIVE_RAW.tsv")), None)
        if component:
            assert sha(component["path"]) == component["sha256"]
            lists[label] = id_rows(Path(component["path"]), "TF")
    items = {}
    for stable, u in universe.items():
        row = native.get(stable, {})
        inp = cov.get(stable, {})
        mog_native = before.get(stable, {})
        source_mean = number(mog_native.get("mean_source_expr"))
        if row:
            reason = ""
        elif status not in AVAILABLE:
            reason = "task_" + status
        elif u["mapping_status"] != "mapped_M24":
            reason = "candidate_ID_unresolved"
        elif not is_true(inp.get("count_present")):
            reason = "no_count_row"
        elif not is_true(inp.get("expressed_row")):
            reason = "all_five_sample_counts_zero"
        else:
            reason = "absent_from_native_five_component_union"
        items[stable] = observation(row.get("RRA_Score"), reason, native=row.get("rank", "NA"), returned=bool(row),
                                    CF_count_present=inp.get("count_present", "NA"), CF_source_mean_raw_counts=inp.get("source_mean_raw_counts", "NA"),
                                    CF_prior_out_edges=inp.get("chipX_prior_out_edges", "NA"), CF_prior_in_edges=inp.get("chipX_prior_in_edges", "NA"),
                                    CF_CellNet_subgraph_targets=nit.get(stable, {}).get("target_n", "NA"),
                                    CF_Mogrify_prefilter=stable in before, CF_Mogrify_source_mean_pass=(source_mean <= 20 if source_mean is not None else "NA"),
                                    CF_Mogrify_redundancy_removed=stable in redundant,
                                    **{f"CF_{c}_native_rank": lists[c].get(stable, {}).get("rank", "NA") if c in lists else "NA"
                                       for c in ("DiffExp_ttest", "DiffExp_deseq", "CellNet", "Mogrify", "dAlessio")})
    return items, status, result.get("reason", run.get("reason", "CompForce component-run record unavailable" if not run else "")), at


def imcell_results(task, universe, rep_tables, stability):
    path = ROOT / "results" / task / "IMCell"
    run, at = attempt(path / "METHOD_RUN.json", "IMCell")
    status = run.get("execution_status", "failed_execution")
    source = read_tsv(path / "NATIVE_SINGLETON_SPREAD.tsv") if (path / "NATIVE_SINGLETON_SPREAD.tsv").is_file() else []
    by_gene = defaultdict(dict)
    rep_rows = defaultdict(list)
    cov = id_rows(ROOT / "inputs_manifest/IMCELL_SEED_COVERAGE.tsv", "stable_id")
    for r in source:
        if r["job_status"] != "completed":
            continue
        stable = r["stable_id"]
        replicate = int(r["replicate"])
        assert replicate not in by_gene[stable] and r["target"] == task
        v = float(r["spread"])
        milli = round(v * 1000)
        assert math.isclose(v, milli / 1000, abs_tol=1e-12, rel_tol=0)
        # The native score is a sum of two integer cascade counts divided by niter=1000.
        # Recover that exact unit to avoid splitting equal three-stream sums by floating-point rounding.
        by_gene[stable][replicate] = milli
        rep_rows[replicate].append({"target": task, "replicate": replicate, "stable_id": stable, "symbol": universe[stable]["gene_symbol"],
                                    "score": Decimal(milli) / 1000, "native_spread_text": r["spread"], "eligible": True,
                                    "job_id": r["job_id"], "seconds": r["seconds"], "candidate_is_activation_program_member": r["seed_in_activation_program"]})
    if status == "completed":
        assert run["expected_jobs"] == run["completed_jobs"] == len(source) == 2718 and len(by_gene) == 906
        assert all(len(v) == 3 for v in by_gene.values())
    rep_maps = {}
    for replicate, rows in sorted(rep_rows.items()):
        full_standalone = status == "completed" and len(rows) == 906
        order = ranks(rows, "descending", full_standalone)
        rep_tables.extend(order)
        rep_maps[replicate] = {r["stable_id"]: r for r in order}
    items = {}
    for stable, u in universe.items():
        in_graph = is_true(cov[stable]["in_network"])
        deg = int(cov[stable]["out_degree"])
        v = by_gene.get(stable, {})
        complete = len(v) == 3 and status == "completed"
        vals = [v[k] / 1000 for k in sorted(v)]
        if complete:
            reason = ""
        elif u["mapping_status"] != "mapped_M24":
            reason = "candidate_ID_unresolved"
        elif not in_graph:
            reason = "candidate_absent_from_frozen_signed_graph"
        elif not deg:
            reason = "graph_node_has_no_outgoing_regulatory_edge"
        else:
            reason = f"task_{status}_or_missing_native_stream"
        replica = [rep_maps[i].get(stable, {}) for i in sorted(rep_maps)] if v else []
        spread = Decimal(sum(v.values())) / 3000 if complete else None
        sdev = statistics.stdev(vals) if complete else None
        items[stable] = observation(spread, reason, returned=bool(v), IMCell_graph_node=in_graph, IMCell_out_degree=deg,
                                    IMCell_observed_replicates=len(v), IMCell_sample_sd=sdev if sdev is not None else "NA")
        if v:
            stability.append({"target": task, "stable_id": stable, "symbol": u["gene_symbol"],
                              "spread_mean": spread if complete else "NA", "spread_sd": sdev if sdev is not None else "NA",
                              "spread_range": max(vals) - min(vals), "top10_replicate_n": sum(is_true(r.get("top10")) for r in replica),
                              "reps_present": len(v), **{f"rep{i}_eval_rank": rep_maps.get(i, {}).get(stable, {}).get("eval_rank", "NA") for i in (1,2,3)}})
    return items, status, run.get("reason", "Native IMCell run record unavailable" if not run else ""), at


def create_rows(task, method, universe, implementation, raw_result, frozen_hash):
    native, status, reason, source = raw_result
    quality, limitations = policy_quality(task, method)
    direction = "ascending" if method == "CompForce RRA" else "descending"
    assert set(native) == set(universe)
    rows = [{"target": task, "method": method, "implementation": implementation, "stable_id": stable,
             "symbol": u["gene_symbol"], "score_direction": direction, "eval_rank": "NA", "rank_min": "NA", "rank_max": "NA",
             "tie_group_size": "NA", "tie_crosses_10": False, "top10": False, "input_QC_status": quality,
             "QC_limitations": limitations, "prior_version": PRIOR[method], "config_hash": frozen_hash, **native[stable]}
            for stable, u in universe.items()]
    eligible_scores = {dscore(r["score"]) for r in rows if r["eligible"] and number(r["score"]) is not None}
    if len(eligible_scores) == 1 and status in AVAILABLE:
        status = "completed_uninformative"
        reason = "All eligible scores equal; deterministic tie ranks retained without recommendations"
    elif not eligible_scores and status in AVAILABLE and status != "completed_empty":
        status = "completed_empty"
        reason = reason or "Native and fixed eligibility rules returned zero rankable candidates"
    usable = status in AVAILABLE and status != "completed_uninformative"
    order = ranks(rows, direction, usability=usable)
    summary = {"target": task, "method": method, "implementation": implementation, "execution_status": status,
               "interpretation_status": quality if usable else "NA", "input_QC_status": quality, "QC_limitations": limitations,
               "scope": "reused" if method in ("Taiji", "VIPER", "TRIAGE") else "exported" if method == "DE" else "new",
               "workflow_outcome_source": source, "reason": reason, "candidate_universe_n": len(rows),
               "scoreable_n": sum(r["scoreable"] for r in rows), "returned_n": sum(r["returned"] for r in rows),
               "eligible_n": len(order), "top10_n": sum(r["top10"] for r in rows),
               "boundary_tie_n": sum(r["tie_crosses_10"] for r in rows), "config_hash": frozen_hash}
    return rows, order, summary


def block_reference_open(event, args):
    if event != "open" or not isinstance(args[0], (bytes, str, os.PathLike)):
        return
    p = Path(os.path.abspath(os.fsdecode(args[0])))
    for parent in (BASE / "references", BASE / "evaluation", ROOT / "evaluation"):
        if p == parent or parent in p.parents:
            raise PermissionError(f"Rank-stage reference access barred: {p}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify-reuse", action="store_true")
    args = parser.parse_args()
    sys.addaudithook(block_reference_open)
    u = id_rows(ROOT / "inputs_manifest/CANDIDATE_UNIVERSE.tsv", "gene_id")
    assert len(u) == 1374
    cfg = ROOT / "FROZEN_CONFIG_ROUND2.yaml"
    frozen_hash = sha(cfg)
    lineage = {(r["target"], r["method"]): r for r in read_tsv(ROOT / "inputs_manifest/REUSE_LINEAGE.tsv")}
    comparisons = {}
    out = ROOT / "results/locked"
    if not args.verify_reuse:
        dispatch = read_json(ROOT / "STATUS.json")
        assert dispatch["dispatcher_complete"] and not dispatch["running"] and dispatch["config_sha256"] == frozen_hash
        lock = ROOT / "PREDICTIONS_LOCK_ROUND2.json"
        if lock.exists():
            for r in read_json(lock)["outputs"]:
                assert sha(r["path"]) == r["sha256"]
            return
    all_rows, variant_rows, tops, variant_top, all_summary, variant_summary = [], [], [], [], [], []
    individual, stability = [], []
    verified_reuse = []
    for task in TASKS:
        for method in METHODS:
            if method == "DE":
                raw = de_results(task, u)
            elif method in ("TRIAGE", "VIPER", "Taiji"):
                raw = old_results(task, method, u, lineage, comparisons)
            elif args.verify_reuse:
                continue
            elif method == "ANANSE":
                raw = ananse_results(task, u)
            elif method == "CompForce RRA":
                raw = compforce_results(task, u)
                alt, alt_order, alt_sum = create_rows(task, method, u, SENSITIVITY, compforce_results(task, u, True), frozen_hash)
                variant_rows.extend(alt)
                variant_top.extend(r for r in alt_order if r["top10"])
                variant_summary.append(alt_sum)
            else:
                raw = imcell_results(task, u, individual, stability)
            rows, order, summary = create_rows(task, method, u, IMPLEMENTATIONS[method], raw, frozen_hash)
            top = [r for r in order if r["top10"]]
            if (task, method) in comparisons:
                prev = comparisons[task, method]
                verified_reuse.append({"target": task, "method": method, "n_original_scores_checked": prev["old_scores_verified"],
                    "n_original_top10_positions_changed": sum(x != y for x, y in zip([r["stable_id"] for r in top], prev["old_top10"])) +
                                                              abs(len(top)-len(prev["old_top10"])),
                    "n_original_finite_ranks_changed": sum(str(r["eval_rank"]) != prev["old_ranks"][r["stable_id"]] for r in rows if r["eligible"]),
                    "old_score_sha256": prev["artifact_sha256"], "eligibility_preserved": True})
            all_rows.extend(rows)
            tops.extend(top)
            all_summary.append(summary)
    if args.verify_reuse:
        write_tsv(ROOT / "validation/REUSE_AND_DE_PRECHECK.tsv", all_summary)
        write_tsv(ROOT / "validation/REUSE_REBUILT_RANK_CHECK.tsv", verified_reuse)
        print("Verified", len(all_summary), "frozen exports and reuse cells without reference labels")
        return

    assert len(all_summary) == 21 and len(all_rows) == 21 * len(u) and len(variant_rows) == 3 * len(u)
    paths = {
        "FULL_SCORES_AND_RANKS.tsv": all_rows, "TOP10_BLIND.tsv": tops, "SCORE_STATUS_21.tsv": all_summary,
        "COMPFORCE_SENSITIVITY_SCORES_AND_RANKS.tsv": variant_rows, "COMPFORCE_SENSITIVITY_TOP10_BLIND.tsv": variant_top,
        "COMPFORCE_SENSITIVITY_STATUS.tsv": variant_summary, "IMCELL_ALL_REPLICATE_RANKS.tsv": individual,
        "IMCELL_SPREAD_STABILITY.tsv": stability, "REUSE_RANK_VERIFICATION.tsv": verified_reuse,
        "FILTER_REASON_COUNTS.tsv": [{"target": t, "method": m, "exclusion_reason": r, "n": n}
                                       for (t,m,r),n in sorted(Counter((x["target"],x["method"],x["exclusion_reason"]) for x in all_rows if x["exclusion_reason"]).items())],
    }
    for name, rows in paths.items():
        fields = None
        if name.endswith("_SCORES_AND_RANKS.tsv") or "TOP10_BLIND" in name:
            fields = [*RANK_FIELDS, *dict.fromkeys(k for r in rows for k in r if k not in RANK_FIELDS)]
        write_tsv(out / name, rows, fields)
    native_root = ROOT / "results"
    native_artifacts = sorted(p for p in native_root.glob("*/*/*") if p.is_file() and
                              (p.name.endswith("_RUN.json") or p.suffix==".tsv" and p.stat().st_size < 6 * 1024**2))
    sources = [*native_artifacts, cfg, ROOT / "SOURCE_MANIFEST.tsv", ROOT / "STATUS.json", Path(__file__)]
    write_json(lock, {"created_at_Singapore": sg(), "config_sha256": frozen_hash, "reference_files_opened": False,
                      "isolation": "rank-stage Python audit hook and explicit input paths; native R wrappers address frozen inputs only",
                      "workflow_scoring_terminal": True, "method_target_n": 21, "full_U_rows": len(all_rows),
                      "outputs": [input_record(out / name, "reference-independent score lock") for name in paths],
                      "lineage_sources": [input_record(p, "reference-independent source or run record") for p in sources]})
    print("Scores and eligibility sealed before reference access", lock, sha(lock))


if __name__ == "__main__":
    main()
