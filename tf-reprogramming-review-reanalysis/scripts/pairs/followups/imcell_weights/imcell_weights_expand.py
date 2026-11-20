"""Expand independent seed simulations and inspect blinded rank stability."""
from collections import Counter, defaultdict
import csv
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import sqrt
from pathlib import Path
from statistics import median, stdev
import sys

inp, out = map(Path, sys.argv[1:3])
mode = sys.argv[3]
ENDPOINTS = ("GSE53951", "GSE114809")
METHODS = ("p010", "in_degree", "mean")
SCORE_FIELDS = (
    "endpoint", "variant", "pair_id", "method", "score", "eligible", "exclusion_reason",
    "tie_key", "score_num", "score_den", "precision_trials", "rank", "rank_min", "rank_max",
    "midrank", "tie_n", "work_id", "bhlh_id", "partner_id", "activation_sum", "repression_sum",
    "mc_se", "mc_rep_sd", "rep1_score_num", "rep2_score_num", "rep3_score_num",
    "n_activate_readouts", "n_repress_readouts",
)

def read(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name, rows, fields=None):
    path = out / name
    with path.with_suffix(".partial").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    path.with_suffix(".partial").replace(path)


def make_scores(acts, reps, squares, n_trials, rep_totals, active, repress):
    n = 3*n_trials
    total = acts + reps
    variance = float(Fraction(squares * n - total**2, n * (n-1)))
    assert variance >= 0
    return {"activation_sum":acts, "repression_sum":reps, "n_trials":n_trials, "num":total,
            "den":n, "rep_totals":rep_totals, "mcse":sqrt(variance/n),
            "n_active":active, "n_repress":repress}


def sums(endpoint, method, precision):
    if precision == 1000 and method == "p010":
        legacy = read(inp / f"mc_block_{ENDPOINTS.index(endpoint)+1}.tsv")
        quantities = {}
        for row in legacy:
            assert row["endpoint"] == endpoint and row["method"] == "imcell_p010" and row["score_den"] == "3000"
            rep_totals = tuple(int(Fraction(row[f"rep{i}"])*1000) for i in range(1,4))
            assert Fraction(row["rep1"])*1000==rep_totals[0] and sum(rep_totals)==int(row["score_num"])
            counts = make_scores(int(row["activation_sum"]),int(row["repression_sum"]),int(row["sum_square_total"]),
                                 precision,rep_totals,int(row["n_active_goals"]),int(row["n_repress_goals"]))
            assert abs(counts["mcse"]-float(row["mc_se"]))<1e-14
            quantities[row["job"]] = counts
        return quantities
    block = read(out / f"mc_{endpoint}_{method}_{precision}.tsv")
    groups = defaultdict(list)
    for row in block:
        assert row["endpoint"] == endpoint and row["method"] == method and int(row["trial_count"]) == precision
        groups[row["job"]].append(row)
    quantities = {}
    for job, rows in groups.items():
        assert [row["replicate"] for row in rows] == ["1","2","3"]
        cols = ("activation_sum","repression_sum","sum_square_total","n_active_goals","n_repress_goals")
        act, repress, squares, n_active, n_repress = [list(map(int,(r[col] for r in rows))) for col in cols]
        assert len(set(n_active))==1 and len(set(n_repress))==1
        quantities[job] = make_scores(sum(act),sum(repress),sum(squares),precision,
                                     tuple(a+b for a,b in zip(act,repress)),n_active[0],n_repress[0])
    assert len(quantities)==296
    return quantities


def ranks(rows, key):
    order = sorted(rows, key=lambda x:(-key(x),x["tie_key"],x["pair_id"]))
    index = 0
    while index < len(order):
        end = index+1
        while end < len(order) and key(order[end]) == key(order[index]):
            end += 1
        for k in range(index,end):
            order[k]["rank"] = k+1
            order[k]["rank_min"] = index+1
            order[k]["rank_max"] = end
            order[k]["midrank"] = (index+1+end)/2
            order[k]["tie_n"] = end-index
        index = end
    return order


def annotate(endpoint, method, n_trials):
    runs = sums(endpoint,method,n_trials)
    slots = read(out/"imcell_expand.tsv")
    records = []
    by_variant = defaultdict(lambda: {"joint":[], "additive":[]})
    for slot in slots:
        meta = {col:slot[col] for col in ("variant","pair_id","tie_key","work_id","bhlh_id","partner_id")}
        meta.update(endpoint=endpoint, eligible=int(not slot["reason"]), exclusion_reason=slot["reason"],
                    precision_trials=n_trials)
        joint = dict(meta, method=f"imcell_weight_{method}_joint")
        add = dict(meta, method=f"imcell_weight_{method}_additive")
        records.extend((joint,add))
        if slot["reason"]:
            continue
        a,b = slot["work_id"].split("|")
        readouts = runs["pair:"+slot["work_id"]]
        for record, keys in ((joint,("pair:"+slot["work_id"],)), (add,("node:"+a,"node:"+b))):
            terms = [runs[x] for x in keys]
            n_total = terms[0]["den"]
            assert all(t["den"]==n_total for t in terms)
            total = sum(t["num"] for t in terms)
            rep_totals = tuple(map(sum,zip(*(t["rep_totals"] for t in terms))))
            assert sum(rep_totals)==total
            record.update(score=float(Fraction(total,n_total)), score_num=total, score_den=n_total,
                          activation_sum=sum(t["activation_sum"] for t in terms),
                          repression_sum=sum(t["repression_sum"] for t in terms),
                          mc_se=sqrt(sum(t["mcse"]**2 for t in terms)),
                          mc_rep_sd=stdev(rep_totals)/n_trials,
                          rep1_score_num=rep_totals[0], rep2_score_num=rep_totals[1], rep3_score_num=rep_totals[2],
                          n_activate_readouts=readouts["n_active"],n_repress_readouts=readouts["n_repress"],
                          simulation_keys=keys)
        by_variant[slot["variant"]]["joint"].append(joint)
        by_variant[slot["variant"]]["additive"].append(add)

    diagnostics = []
    for variant, two_models in by_variant.items():
        for model, block in two_models.items():
            ordered = ranks(block, lambda r:r["score_num"])
            for budget in (10,20):
                left,right = ordered[budget-1:budget+1]
                gap = Fraction(left["score_num"]-right["score_num"],left["score_den"])
                assert gap>=0 and left["score_den"]==right["score_den"]
                subtotal = Counter(left["simulation_keys"])
                subtotal.subtract(right["simulation_keys"])
                combined_se = sqrt(sum((weight*runs[k]["mcse"])**2 for k,weight in subtotal.items()))
                finalist = sorted(rep for rep in ("rep1_score_num", "rep2_score_num","rep3_score_num"))
                rep_sets = [
                    {r["pair_id"] for r in sorted(block, key=lambda v:(-v[rep],v["tie_key"],v["pair_id"]))[:budget]}
                    for rep in finalist
                ]
                jac = min(len(a & b)/len(a | b) for a,b in combinations(rep_sets,2))
                diag = dict(endpoint=endpoint, variant=variant, weight=method, score_model=model,
                            precision_trials=n_trials, budget=budget, n_scoreable=len(block),
                            boundary_pair_id=left["pair_id"], next_pair_id=right["pair_id"],
                            boundary_gap=float(gap), strict_tie=int(gap==0),
                            shared_terms_cancelled=len(set(left["simulation_keys"]) & set(right["simulation_keys"])),
                            combined_mc_se=combined_se,
                            within_1p96_se=int(float(gap)<=1.96*combined_se),
                            min_rep_top_jaccard=jac, replicate_top_unstable=int(jac<1),
                            median_mc_se=median(r["mc_se"] for r in block))
                for rep, ids in enumerate(rep_sets,1):
                    diag[f"rep{rep}_top_pair_ids"] = "|".join(sorted(ids))
                diagnostics.append(diag)
    for record in records:
        record.pop("simulation_keys",None)
    return records, diagnostics


def run():
    complete_initial = []
    stability_initial = []
    receipts = []
    for endpoint in ENDPOINTS:
        comparisons = []
        for method in METHODS:
            rows, diag = annotate(endpoint,method,1000)
            complete_initial.extend(rows)
            comparisons.extend(diag)
        stability_initial.extend(comparisons)
        trigger_n = sum(r["within_1p96_se"] or r["replicate_top_unstable"] for r in comparisons)
        receipts.append({"endpoint":endpoint,"trigger":int(trigger_n>0),"n_boundary_checks":len(comparisons),
                         "n_triggered":trigger_n,"trigger_rule":"any unlabeled variant x weight x J/A x B=10/20 flags",
                         "n_trials_each_replicate_initial":1000,
                         "n_trials_each_replicate_final":5000 if trigger_n else 1000})
    write("IMCELL_PRECISION_TRIGGER.tsv",receipts)
    write("IMCELL_WEIGHT_SCORES_INITIAL.tsv",complete_initial,SCORE_FIELDS)
    if mode == "trigger":
        write("IMCELL_MC_STABILITY.tsv",stability_initial)
        return
    final = []
    full_stability = stability_initial.copy()
    for r in receipts:
        endpoint,precision = r["endpoint"],r["n_trials_each_replicate_final"]
        if precision == 1000:
            final.extend(x for x in complete_initial if x["endpoint"]==endpoint)
            continue
        for method in METHODS:
            records, diagnostics = annotate(endpoint,method,precision)
            final.extend(records)
            full_stability.extend(diagnostics)
    assert len(final)==598*9*2*3*2
    assert len({tuple(x[k] for k in ("endpoint","variant","pair_id","method")) for x in final})==len(final)
    write("IMCELL_WEIGHT_SCORES.tsv",final,SCORE_FIELDS)
    write("IMCELL_MC_STABILITY.tsv",full_stability)


def verify():
    scores = read(out/"IMCELL_WEIGHT_SCORES.tsv")
    by_model = defaultdict(dict)
    for r in scores:
        by_model[(r["endpoint"],r["variant"],r["method"])][r["pair_id"]] = r
    design = [r for r in read(out/"imcell_expand.tsv") if r["variant"]=="main" and not r["reason"]]
    sample = sorted(design,key=lambda r:sha256(("section4-weight-v1|"+r["pair_id"]).encode()).digest())[:12]
    checks = []
    for endpoint in ENDPOINTS:
        for weight in METHODS:
            j_group = by_model[(endpoint,"main",f"imcell_weight_{weight}_joint")]
            a_group = by_model[(endpoint,"main",f"imcell_weight_{weight}_additive")]
            precision = int(next(iter(j_group.values()))["precision_trials"])
            if precision == 1000 and weight == "p010":
                path = inp / f"mc_block_{ENDPOINTS.index(endpoint)+1}.tsv"
                raw_ratios = {r["job"]:Fraction(int(r["score_num"]),int(r["score_den"])) for r in read(path)}
            else:
                path = out / f"mc_{endpoint}_{weight}_{precision}.tsv"
                raw_parts = defaultdict(Fraction)
                for r in read(path):
                    raw_parts[r["job"]] += Fraction(int(r["activation_sum"])+int(r["repression_sum"]),3*int(r["trial_count"]))
                raw_ratios = dict(raw_parts)

            def exact_score(r, joint):
                if joint:
                    return raw_ratios["pair:"+r["work_id"]]
                left,right = r["work_id"].split("|")
                return raw_ratios["node:"+left] + raw_ratios["node:"+right]

            def validate_rank(group, joint, pair_id):
                rows = [r for r in group.values() if r["eligible"]=="1"]
                current = group[pair_id]
                value = exact_score(current,joint)
                better = [r for r in rows if exact_score(r,joint) > value]
                eq = [r for r in rows if exact_score(r,joint) == value]
                tie_key = (current["tie_key"],pair_id)
                rank = len(better)+sum((r["tie_key"],r["pair_id"])<tie_key for r in eq)+1
                assert rank==int(current["rank"])
                assert len(better)+1==int(current["rank_min"]) and len(better)+len(eq)==int(current["rank_max"])
                assert value==Fraction(int(current["score_num"]),int(current["score_den"]))
                assert abs(float(value)-float(current["score"]))<=max(abs(float(value))*1e-15,1e-15)
                return value,rank

            for item in sample:
                pair_id = item["pair_id"]
                j, j_rank = validate_rank(j_group,True,pair_id)
                a, a_rank = validate_rank(a_group,False,pair_id)
                checks.append({"endpoint":endpoint,"weight":weight,"precision_trials":precision,
                               "pair_id":pair_id,"canonical_seed_key":item["work_id"], "joint_exact":str(j),
                               "additive_exact":str(a),"joint_minus_additive_exact":str(j-a),
                               "recomputed_joint_rank":j_rank,"recomputed_additive_rank":a_rank,
                               "source_checkpoint":path.name,"pass":True})
    write("IMCELL_INDEPENDENT_REBUILD.tsv",checks)


if mode=="verify":
    verify()
else:
    assert mode in {"trigger","publish"}
    run()

