"""Verify the blinded score lock, then evaluate exact references and independent motif-family projections."""
import itertools
import statistics
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

from r2lib import BASE, METHODS, ROOT, TASKS, input_record, read_json, read_tsv, sg, sha, write_json, write_tsv


LOCKED = ROOT / "results/locked"
OUT = ROOT / "tables"
EVAL = ROOT / "evaluation"
COMPLETED = {"reused_verified", "exported_verified", "completed", "completed_empty"}


def truth(x):
    return str(x).lower() == "true"


def grouped(rows, *keys):
    result = defaultdict(list)
    for row in rows:
        result[tuple(row[key] for key in keys)].append(row)
    return result


def names(ids, symbol):
    return ";".join(symbol.get(g, g) or g for g in sorted(ids, key=lambda g: (symbol.get(g, g), g)))


def shortlist(rows, status):
    direction, = {r["score_direction"] for r in rows}
    assert direction in ("descending", "ascending")
    ranks = sorted((r for r in rows if truth(r["eligible"]) and r["score"] not in ("NA", "")),
                   key=lambda r: (Decimal(r["score"]) * (-1 if direction == "descending" else 1), r["stable_id"]))
    counts = Counter(Decimal(r["score"]) for r in ranks)
    values = len(ranks)
    rank = 1
    for value, iterator in itertools.groupby(ranks, key=lambda r: Decimal(r["score"])):
        group = list(iterator)
        max_rank = rank + counts[value] - 1
        for i, r in enumerate(group, start=rank):
            assert int(r["eval_rank"]) == i and int(r["rank_min"]) == rank and int(r["rank_max"]) == max_rank
            assert int(r["tie_group_size"]) == counts[value] and truth(r["tie_crosses_10"]) == (rank <= 10 < max_rank)
        rank = max_rank + 1
    evaluated = status in COMPLETED
    top = ranks[:10] if evaluated else []
    assert {r["stable_id"] for r in top} == {r["stable_id"] for r in rows if truth(r["top10"])}
    assert values == sum(r["eval_rank"] not in ("NA", "") for r in rows)
    return top


def verified_refs():
    artifacts = ("REFERENCE_TFS.tsv", "REFERENCE_EVIDENCE.tsv", "MOTIF_MAPPING.tsv")
    previous = {Path(r["path"]).name: r["sha256"] for r in read_tsv(BASE / "FREEZE_MANIFEST.tsv") if r["role"] == "evaluation_only"}
    result = {}
    for name in artifacts:
        path = BASE / "references" / name
        assert sha(path) == previous[name], path
        result[name] = input_record(path, "opened by evaluator after prediction lock")
    reference = [r for r in read_tsv(result[artifacts[0]]["path"]) if r["task_id"] in TASKS]
    reference_by_target = {task: sorted((r for r in reference if r["task_id"] == task), key=lambda r: int(r["reference_order"])) for task in TASKS}
    assert [len(reference_by_target[t]) for t in TASKS] == [11, 6, 5]
    assert all(len({r["gene_id"] for r in rs}) == len(rs) for rs in reference_by_target.values())
    motif = read_tsv(result[artifacts[2]]["path"])
    families = defaultdict(set)
    for r in motif:
        if r["gene_id"].startswith("ENSMUSG"):
            families[r["gene_id"]].add(r["cluster"])
    assert families
    evidence = [r for r in read_tsv(result[artifacts[1]]["path"]) if r["task_id"] in TASKS]
    reporting = BASE / "references/REFERENCE_EVIDENCE_REPORTING_VIEW.tsv"
    report_record = read_tsv(BASE / "OUTPUT_MANIFEST.tsv")
    published_sha, = [r["sha256"] for r in report_record if Path(r["path"]).name == reporting.name]
    assert sha(reporting) == published_sha
    context = [r for r in read_tsv(reporting) if r["task_id"] in TASKS]
    result[reporting.name] = input_record(reporting, "reporting context only, no membership effects")
    EVAL.mkdir(exist_ok=True)
    write_tsv(EVAL / "REFERENCE_TFS_FROZEN_ROUND2.tsv", reference)
    write_tsv(EVAL / "REFERENCE_EVIDENCE_FROZEN_ROUND2.tsv", evidence)
    write_tsv(EVAL / "REFERENCE_EVIDENCE_REPORTING_ROUND2.tsv", context)
    write_tsv(EVAL / "MOTIF_MAPPING_FROZEN.tsv", motif)
    write_tsv(EVAL / "REFERENCE_MANIFEST.tsv", result.values())
    return reference_by_target, families


def evaluate_group(task, method, group, state, ref_rows, families, universe):
    status = state["execution_status"]
    evaluated = status in COMPLETED
    top = shortlist(group, status)
    symbols = {r["stable_id"]: r["symbol"] for r in group}
    symbols.update({r["gene_id"]: r["gene_symbol"] for r in ref_rows})
    reference = {r["gene_id"] for r in ref_rows}
    index = {r["stable_id"]: r for r in group}
    eligible = {g for g, r in index.items() if truth(r["eligible"]) and r["score"] != "NA"}
    scoreable = {g for g, r in index.items() if truth(r["scoreable"])}
    assert len(eligible) == int(state["eligible_n"])
    top_ids = {r["stable_id"] for r in top}
    actual = reference & top_ids
    family_ref = set().union(*(families[g] for g in reference))
    motif_top = set().union(*(families[g] for g in top_ids))
    motif_hits = motif_top & family_ref
    summary = {**state, "reference_total_n": len(reference), "reference_in_U_n": len(reference & universe),
               "reference_scoreable_n": len(reference & scoreable) if evaluated else "NA",
               "reference_eligible_n": len(reference & eligible) if evaluated else "NA", "actual_recommendations_n": len(top),
               "Hits_at_10": len(actual) if evaluated else "NA", "exact_hits": names(actual, symbols) if evaluated else "NA",
               "exact_hits_in_rank_order": ";".join(symbols[r["stable_id"]] for r in top if r["stable_id"] in actual) if evaluated else "NA",
               "motif_family_Hits_at_10": len(motif_hits) if evaluated else "NA",
               "motif_family_reference_n": len(family_ref), "motif_family_hits": ";".join(sorted(motif_hits)) if evaluated else "NA",
               "motif_family_mapping_top10_n": sum(bool(families[g]) for g in top_ids) if evaluated else "NA"}
    matrix = []
    for ref in ref_rows:
        g = ref["gene_id"]
        r = index.get(g, {})
        exact = g in actual
        if g not in universe:
            state_code = "outside_U"
            reason = "Frozen reference gene is outside official candidate universe"
        elif not evaluated:
            state_code = "task_NA"
            reason = state["reason"] or status
        elif exact:
            state_code, reason = "hit_top10", ""
        elif g in eligible:
            state_code, reason = "ranked_beyond_10", ""
        elif g in scoreable:
            state_code, reason = "scored_native_ineligible", r.get("exclusion_reason", "")
        else:
            state_code, reason = "candidate_unscoreable", r.get("exclusion_reason", "")
        matrix.append({"target": task, "method": method, "implementation": state["implementation"], "stable_id": g,
                       "reference_symbol": ref["gene_symbol"], "reference_order": ref["reference_order"], "in_frozen_U": g in universe,
                       "candidate_status": state_code, "scoreable": g in scoreable, "eligible": g in eligible,
                       "score": r.get("score", "NA"), "native_rank": r.get("native_rank", "NA"), "eval_rank": r.get("eval_rank", "NA"),
                       "rank_min": r.get("rank_min", "NA"), "rank_max": r.get("rank_max", "NA"),
                       "touches_split_tie_at_10": r.get("tie_crosses_10", "False"), "hits_exact_top10": exact if evaluated else "NA",
                       "execution_status": status, "interpretation_status": state["interpretation_status"], "reason": reason})
    top_rows = [{**{k: r[k] for k in ("target","method","implementation","stable_id","symbol","score","eval_rank","native_rank","rank_min","rank_max","tie_crosses_10")},
                "exact_reference_hit": r["stable_id"] in actual,
                "assigned_frozen_families": ";".join(sorted(families[r["stable_id"]])),
                "matching_reference_families": ";".join(sorted(families[r["stable_id"]] & family_ref))} for r in top]
    per_family = [{"target": task,"method": method,"implementation": state["implementation"],"family": f,
                  "reference_members": names({g for g in reference if f in families[g]}, symbols),
                  "recommended_family_members": names({g for g in top_ids if f in families[g]}, symbols) if evaluated else "NA",
                  "family_hit": f in motif_hits if evaluated else "NA","membership_semantics": "unique cluster, separate from exact gene"}
                  for f in sorted(family_ref)]
    return summary, matrix, top_rows, per_family, actual if evaluated else None, top_ids if evaluated else None


def supplementary(sums, group_map, sensitivity_map, identity_summary, refs):
    sens_summary, detail = [], []
    cf_component, ananse_qc, edge_cutoffs = [], [], []
    for t in TASKS:
        main_state = sums[t, "CompForce RRA"]
        sensitivity_state, _, _, _, alt_hits, alt_ids = sensitivity_map[t]
        primary = ROOT / "results" / t / "CompForce/METHOD_RUN.json"
        comp = read_json(primary) if primary.is_file() else {}
        same = comp.get("identical_other_four_components", {})
        main_hits, main_ids = identity_summary[t, "CompForce RRA"]
        identity = {r["gene_id"] for r in refs[t]}
        main = {r["stable_id"]: r for r in group_map[t,"CompForce RRA"]}
        alt = {r["stable_id"]: r for r in read_tsv(LOCKED / "COMPFORCE_SENSITIVITY_SCORES_AND_RANKS.tsv") if r["target"] == t}
        if main_ids is not None and alt_ids is not None:
            main_only = main_ids - alt_ids
            alt_only = alt_ids - main_ids
            for g in sorted(main_ids | alt_ids):
                detail.append({"target": t,"stable_id": g,"symbol": main[g]["symbol"],"archived_rank": main[g]["eval_rank"],
                               "normalized_ttest_rank": alt[g]["eval_rank"], "archived_top10": g in main_ids,"normalized_ttest_top10": g in alt_ids,
                               "frozen_reference": g in identity,"archived_RRA_Score": main[g]["score"],"normalized_ttest_RRA_Score": alt[g]["score"]})
            selection_changed = main_ids != alt_ids
            n_overlap = len(main_ids & alt_ids)
            original_only = names(main_only, {g:r["symbol"] for g,r in main.items()})
            variant_only = names(alt_only, {g:r["symbol"] for g,r in alt.items()})
        else:
            selection_changed = n_overlap = original_only = variant_only = "NA"
        sens_summary.append({"target": t,"archived_status": main_state["execution_status"],"logCP100k_status": sensitivity_state["execution_status"],
                             "archived_nonempty_component_n": comp.get("primary",{}).get("native_nonempty_list_n", "NA"),
                             "sensitivity_nonempty_component_n": comp.get("sensitivity",{}).get("native_nonempty_list_n", "NA"),
                             "four_other_components_bytes_verified": len(same)==4 and all(c["component_files_available"] and c["export_sha256_verified_if_available"] and c["preserved_same_export_paths"] for c in same.values()),
                             "top10_identity_changed": selection_changed,"top10_identity_intersection_n": n_overlap,
                             "archived_only_recommendations": original_only,"logCP100k_only_recommendations": variant_only,
                             "archived_Hits_at_10": main_state["Hits_at_10"],"logCP100k_Hits_at_10": sensitivity_state["Hits_at_10"],
                             "archived_reference_hits": main_state["exact_hits"],"logCP100k_reference_hits": sensitivity_state["exact_hits"]})
        for variant, key in (('primary','primary'),('logCP100k','sensitivity')):
            for name, inp in comp.get(key, {}).get("inputs", {}).items():
                cf_component.append({"target": t,"design": variant,"component": name,"state": inp["state"],
                                    "returned_n": inp.get("number","NA"),"finite_native_metric_n": inp.get("number_finite","NA"),
                                    "reason": inp.get("reason",""),"export_hashes": ";".join(f["sha256"] for f in inp.get("file_sha256",[]))})
        for name in ('binding_source','binding_target','network_source','network_target','influence'):
            path = ROOT / "results" / t / "ANANSE" / f"{name}_RUN.json"
            r = read_json(path) if path.is_file() else {}
            argv = r.get("argv",[])
            n_slots = r.get("native_slot_option_n",int(argv[argv.index("-n")+1]) if "-n" in argv else "NA")
            ananse_qc.append({"target": t,"stage": name,"execution_status": r.get("execution_status","blocked_dependency"),
                              "shared_regions": r.get("region_n","NA"),"native_binding_TFs": len(r.get("finite_binding_factors",[])) if "finite_binding_factors" in r else "NA",
                              "network_edges": r.get("network_edges","NA"),"network_nodes": r.get("graph_nodes","NA"),
                              "network_regulators": r.get("regulators","NA"),"native_influence_TFs": r.get("influence_factors","NA"),
                              "actual_native_n_option":n_slots,"operator_threads_for_network": 2*n_slots if name.startswith("network") and isinstance(n_slots,int) else "NA",
                              "Dask_worker_memory_limit":r.get("dask_worker_memory","NA") if name.startswith("network") else "NA",
                              "run_record": str(path) if r else "NA","reason": r.get("reason","")})
        cutoff = ROOT/"results"/t/"ANANSE/INFLUENCE_NETWORK_BOUNDARY.tsv"
        if cutoff.is_file():
            edge_cutoffs.extend({"target":t,**r} for r in read_tsv(cutoff))
    write_tsv(OUT / "COMPFORCE_SENSITIVITY_COMPARISON.tsv", sens_summary)
    write_tsv(OUT / "COMPFORCE_SENSITIVITY_IDENTITY_DIFF.tsv", detail)
    write_tsv(OUT / "COMPFORCE_COMPONENTS.tsv", cf_component)
    write_tsv(OUT / "ANANSE_NATIVE_QC.tsv", ananse_qc)
    write_tsv(OUT / "ANANSE_INFLUENCE_EDGE_CUTOFFS.tsv", edge_cutoffs,
              ["target","state","total_edges","kept_edges","requested_cutoff","boundary_prob","edges_tied_at_boundary","tied_edges_selected","split_boundary_tie"])
    repeat = grouped(read_tsv(LOCKED / "IMCELL_ALL_REPLICATE_RANKS.tsv"),"target","replicate")
    stable = grouped(read_tsv(LOCKED / "IMCELL_SPREAD_STABILITY.tsv"),"target")
    stabs, windows = [], []
    for t in TASKS:
        _, average = identity_summary[t,"IMCell"]
        if average is None:
            stabs.append({"target": t,"status": sums[t,"IMCell"]["execution_status"]})
            continue
        ranked = {r["stable_id"]: r for r in group_map[t,"IMCell"] if r["eval_rank"] != "NA"}
        by_rank = sorted(ranked.values(),key=lambda r:int(r["eval_rank"]))
        n_boundary = float(by_rank[9]["score"])-float(by_rank[10]["score"]) if len(by_rank)>10 else "NA"
        sets = {i:{r["stable_id"] for r in repeat[t,str(i)] if truth(r["top10"])} for i in (1,2,3)}
        spreads = [float(r["spread_sd"]) for r in stable[t,] if r["spread_sd"]!="NA"]
        rank_order_symbols = lambda genes: ";".join(r["symbol"] for r in by_rank if r["stable_id"] in genes)
        row = {"target": t,"status":"completed","scoreable_n":len(ranked),"median_MC_sample_sd":statistics.median(spreads),
               "maximum_MC_sample_sd":max(spreads),"rank10_minus_rank11_mean_spread":n_boundary,
               "repeat12_top10_overlap_n":len(sets[1]&sets[2]),"repeat13_top10_overlap_n":len(sets[1]&sets[3]),
               "repeat23_top10_overlap_n":len(sets[2]&sets[3]),"three_repeat_top10_intersection":rank_order_symbols(sets[1]&sets[2]&sets[3]),
               "three_repeat_top10_union":rank_order_symbols(sets[1]|sets[2]|sets[3]),
               "top10_union_n":len(sets[1]|sets[2]|sets[3]),"top10_present_in_all_three_n":len(sets[1]&sets[2]&sets[3]),
               **{f"final_vs_rep{i}_top10_overlap_n":len(average&sets[i]) for i in (1,2,3)}}
        stabs.append(row)
        for r in stable[t,]:
            master = ranked[r["stable_id"]]
            if 7 <= int(master["eval_rank"]) <= 14 or r["stable_id"] in set().union(*sets.values()):
                windows.append({**r,"mean_rank":master["eval_rank"],"in_final_top10": truth(master["top10"]),
                                "is_boundary_10_or_11":int(master["eval_rank"]) in (10,11)})
    write_tsv(OUT / "IMCELL_MC_STABILITY_SUMMARY.tsv", stabs)
    write_tsv(OUT / "IMCELL_MC_BOUNDARY.tsv", windows)


def main():
    lock_path = ROOT / "PREDICTIONS_LOCK_ROUND2.json"
    lock = read_json(lock_path)
    assert lock["workflow_scoring_terminal"] and not lock["reference_files_opened"] and lock["full_U_rows"]==28854
    assert sha(ROOT / "FROZEN_CONFIG_ROUND2.yaml") == lock["config_sha256"]
    for r in lock["outputs"] + lock["lineage_sources"]:
        assert sha(r["path"]) == r["sha256"], r["path"]
    full = grouped(read_tsv(LOCKED / "FULL_SCORES_AND_RANKS.tsv"),"target","method")
    states = {(r["target"],r["method"]): r for r in read_tsv(LOCKED / "SCORE_STATUS_21.tsv")}
    universe = {r["gene_id"] for r in read_tsv(ROOT / "inputs_manifest/CANDIDATE_UNIVERSE.tsv")}
    assert len(universe) == 1374 and len(full) == len(states) == 21
    for key, group in full.items():
        assert len(group)==len(universe) and {r["stable_id"] for r in group}==universe and all(r["config_hash"]==lock["config_sha256"] for r in group)
    ref, families = verified_refs()
    assert all(truth(r["in_frozen_candidate_U"]) == (r["gene_id"] in universe) for rs in ref.values() for r in rs)
    s = grouped(read_tsv(LOCKED / "COMPFORCE_SENSITIVITY_SCORES_AND_RANKS.tsv"),"target")
    ss = {r["target"]: r for r in read_tsv(LOCKED / "COMPFORCE_SENSITIVITY_STATUS.tsv")}
    summaries, matrix, tops, projection = [], [], [], []
    all_sets, o_summaries, sensitivity = {}, {}, {}
    for t in TASKS:
        sensitivity[t] = evaluate_group(t,"CompForce RRA",s[t,],ss[t],ref[t],families,universe)
        wide, wide_status = [], []
        for m in METHODS:
            key = t,m
            summ, ranks, top, proj, genes, shortlist_ids = evaluate_group(t,m,full[key],states[key],ref[t],families,universe)
            summaries.append(summ); matrix.extend(ranks); tops.extend(top); projection.extend(proj)
            all_sets[key] = (genes,shortlist_ids); o_summaries[key] = summ
            wide.append({"method":m, **{r["reference_symbol"]:r["eval_rank"] if r["candidate_status"] in ("hit_top10","ranked_beyond_10") else r["candidate_status"] for r in ranks},
                         "Hits_at_10":summ["Hits_at_10"],"actual_recommendations_n":summ["actual_recommendations_n"]})
            wide_status.append({"method":m,**{r["reference_symbol"]: r["candidate_status"] for r in ranks}})
        write_tsv(OUT / f"REFERENCE_RANK_MATRIX_{t}.tsv", wide)
        write_tsv(OUT / f"REFERENCE_HIT_IDENTITY_MATRIX_{t}.tsv", wide_status)
    raw_top = read_tsv(LOCKED / "TOP10_BLIND.tsv")
    assert [(r["target"],r["method"],r["stable_id"],str(r["eval_rank"])) for r in tops] == [(r["target"],r["method"],r["stable_id"],r["eval_rank"]) for r in raw_top]
    pairwise, unique, totals, boundaries, multiplicity = [], [], [], [], []
    for t in TASKS:
        reference = {r["gene_id"] for r in ref[t]}
        symbol = {r["stable_id"]:r["symbol"] for r in full[t,"DE"]}
        symbol.update({r["gene_id"]:r["gene_symbol"] for r in ref[t]})
        ok = [m for m in METHODS if all_sets[t,m][0] is not None]
        union = set().union(*(all_sets[t,m][0] for m in ok))
        common = set.intersection(*(all_sets[t,m][0] for m in ok)) if ok else set()
        multis = Counter(g for m in ok for g in all_sets[t,m][0])
        for g in sorted(reference):
            obs = [m for m in ok if g in all_sets[t,m][0]]
            multiplicity.append({"target":t,"stable_id":g,"symbol":symbol[g],"methods_recovering": ";".join(obs),
                                 "method_n":len(obs),"evaluable_methods_n":len(ok),"outside_U":g not in universe,
                                 "status":"common" if obs and len(obs)==len(ok) else "unique" if len(obs)==1 else "shared" if obs else "unrecovered"})
        totals.append({"target":t,"reference_total_n":len(reference),"evaluable_methods_n":len(ok),"evaluable_methods":";".join(ok),
                       "union_reference_n":len(union) if ok else "NA","union_reference":names(union,symbol) if ok else "NA",
                       "common_to_all_evaluable_n":len(common) if ok else "NA","common_to_all_evaluable":names(common,symbol) if ok else "NA",
                       "shared_by_at_least_two":names({g for g,n in multis.items() if n>=2},symbol) if ok else "NA",
                       "unrecovered_by_evaluable_n":len(reference-union) if ok else "NA",
                       "unrecovered_by_evaluable":names(reference-union,symbol) if ok else "NA","frozen_reference_outside_U":names(reference-universe,symbol)})
        for m in METHODS:
            gs = all_sets[t,m][0]
            only = {g for g in gs if multis[g]==1} if gs is not None else None
            unique.append({"target":t,"method":m,"unique_exact_hits_n":len(only) if only is not None else "NA",
                           "unique_exact_hits":names(only,symbol) if only is not None else "NA","comparison_methods_n":len(ok)})
            boundaries.extend({"target":t,"method":m,"stable_id":r["stable_id"],"symbol":r["symbol"],"score":r["score"],"eval_rank":r["eval_rank"],
                            "rank_min":r["rank_min"],"rank_max":r["rank_max"],"tie_crosses_10":r["tie_crosses_10"],"top10":r["top10"]}
                            for r in full[t,m] if truth(r["tie_crosses_10"]) or r["eval_rank"] in ("9","10","11","12"))
        for a,b in itertools.combinations(METHODS,2):
            ma, ta = all_sets[t,a]
            mb, tb = all_sets[t,b]
            passed = ma is not None and mb is not None
            pairwise.append({"target":t,"method_1":a,"method_2":b,"both_evaluable":passed,
                "hits_intersection_n":len(ma&mb) if passed else "NA","hits_intersection":names(ma&mb,symbol) if passed else "NA",
                "hits_union":names(ma|mb,symbol) if passed else "NA","hits_only_method_1":names(ma-mb,symbol) if passed else "NA",
                "hits_only_method_2":names(mb-ma,symbol) if passed else "NA","candidate_top10_intersection_n":len(ta&tb) if passed else "NA",
                "candidate_top10_intersection":names(ta&tb,symbol) if passed else "NA",
                "candidate_top10_jaccard":len(ta&tb)/len(ta|tb) if passed and ta|tb else "NA"})
    assert len(matrix)==7*22 and len(pairwise)==3*21
    output = {"HITS_AT_10_AND_STATUS_21.tsv": summaries, "REFERENCE_RANKS_LONG.tsv": matrix, "TOP10_IDENTITIES.tsv": tops,
              "FAMILY_RECOVERY.tsv": projection,"PAIRWISE_IDENTITY_OVERLAP.tsv":pairwise,"REFERENCE_IDENTITY_UNIQUE.tsv":unique,
              "REFERENCE_RECOVERY_UNION_AND_MISSING.tsv":totals,"REFERENCE_RECOVERY_BY_IDENTITY.tsv":multiplicity,"RANK10_BOUNDARIES.tsv":boundaries,
              "COMPFORCE_SENSITIVITY_REFERENCE_RANKS.tsv":[r for t in TASKS for r in sensitivity[t][1]]}
    for p, rows in output.items():
        write_tsv(OUT/p, rows)
    supplementary(o_summaries,full,sensitivity,all_sets,ref)
    write_json(EVAL / "EVALUATION_VALIDATION.json", {"evaluation_time_Singapore": sg(),
         "score_lock_sha256": sha(lock_path),"prediction_scores_checked_before_reference_read":True,
         "full_U_gated_and_independently_sorted":True,"top10_order_matches_blind_artifact":True,"exact_gene_ranks_written":len(matrix),
         "method_target_states":len(summaries),"references_by_target":{t:len(ref[t]) for t in TASKS},"family_deduplicated": True,
         "evaluation_outputs": [input_record(p,"evaluated result table") for p in sorted(OUT.glob("*.tsv"))]})
    print("Completed independent evaluation for", len(summaries), "target–method tasks across", sum(map(len,ref.values())), "target-associated exact references")


if __name__ == "__main__":
    main()
