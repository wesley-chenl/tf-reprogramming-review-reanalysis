"""Expand joint-seed scores to all unlabeled physical-slot hypotheses."""
from collections import defaultdict
import csv
from fractions import Fraction
from pathlib import Path
from statistics import median
from math import sqrt
import sys

directory = Path(sys.argv[1])
def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

models = [r for block in range(1,7) for r in read(directory / f"mc_block_{block}.tsv")]
slots = read(directory/"imcell_expand.tsv")
by_model = defaultdict(dict)
for row in models:
    by_model[(row["endpoint"],row["method"])][row["job"]] = row
assert len(by_model) == 6
scored = []
counts = []
boundary = []
for (endpoint,method), runs in by_model.items():
    groups = defaultdict(list)
    for slot in slots:
        values = {"method":method,"endpoint":endpoint,"variant":slot["variant"],"pair_id":slot["pair_id"],
                  "reason":slot["reason"],"details":"fixed_probability_signed_TRRUST_native_joint_kernel",
                  "tie_key":slot["tie_key"]}
        if not slot["reason"]:
            left,right = slot["work_id"].split("|")
            joint = runs["pair:"+slot["work_id"]]
            a,b = [runs["node:"+symbol] for symbol in (left,right)]
            single_a = Fraction(int(a["score_num"]),int(a["score_den"]))
            single_b = Fraction(int(b["score_num"]),int(b["score_den"]))
            if left != slot["bhlh_id"].removeprefix("GeneID:"):
                single_a,single_b = single_b,single_a
            values |= {k:joint[k] for k in ("score_num","score_den","activation_sum","repression_sum","sum_square_total","mc_se",
                                           "rep1","rep2","rep3","n_active_goals","n_repress_goals")}
            values["score"] = str(float(Fraction(int(values["score_num"]), int(values["score_den"]))))
            values["bhlh_single"] = str(float(single_a))
            values["partner_single"] = str(float(single_b))
            values["single_sum"] = str(float(single_a+single_b))
            values["joint_minus_single_sum"] = str(float(Fraction(int(values["score_num"]),int(values["score_den"])) - single_a - single_b))
            groups[slot["variant"]].append(values)
        scored.append(values)
    for variant, block in groups.items():
        counts.append({"endpoint":endpoint,"method":method,"variant":variant,"n_scoreable":len(block),
                       "distinct_mc_means":len({(x["score_num"],x["score_den"]) for x in block}),
                       "n_nonzero":sum(float(x["score"])>0 for x in block),
                       "median_mc_se":median(float(x["mc_se"]) for x in block)})
        ordered = sorted(block,key=lambda x:(-Fraction(int(x["score_num"]),int(x["score_den"])),x["tie_key"],x["pair_id"]))
        additive = sorted(block,key=lambda x:(-float(x["single_sum"]),x["tie_key"],x["pair_id"]))
        for bval in (10,20):
            a,b = ordered[bval-1:bval+1]
            gap = float(a["score"])-float(b["score"])
            noise = sqrt(float(a["mc_se"])**2+float(b["mc_se"])**2)
            sets = [set(x["pair_id"] for x in sorted(block,key=lambda p:(-float(p[f"rep{i}"]),p["tie_key"],p["pair_id"]))[:bval])
                    for i in range(1,4)]
            jac = min(len(x & y)/len(x | y) for k,x in enumerate(sets) for y in sets[k+1:])
            boundary.append({"endpoint":endpoint,"method":method,"variant":variant,"b":bval,
                             "boundary_score_gap":gap,"combined_mc_se":noise,"criterion":"gap <= 1.96 * combined_mc_se",
                             "boundary_unstable":int(gap<=1.96 * noise),"minimum_replicate_top_jaccard":jac,
                             "joint_vs_additive_top_overlap":len(set(x["pair_id"] for x in ordered[:bval]) & set(x["pair_id"] for x in additive[:bval]))})

write(directory/"imcell.tsv",scored,("method","endpoint","variant","pair_id","score","reason","details","tie_key","score_num","score_den",
  "activation_sum","repression_sum","sum_square_total","mc_se","rep1","rep2","rep3","n_active_goals","n_repress_goals",
  "bhlh_single","partner_single","single_sum","joint_minus_single_sum"))
write(directory/"IMCELL_COVERAGE.tsv",counts,list(counts[0]))
write(directory/"IMCELL_MC_BOUNDARY.tsv",boundary,list(boundary[0]))
print(f"expanded {len(scored)} IMCell rows and {len(boundary)} MC stability diagnostics")
