"""Reconstruct Taiji K=2 gate and products using independent Decimal arithmetic."""
import argparse
import csv
from decimal import Decimal, localcontext
from hashlib import sha256
import json
from pathlib import Path

def table(p):
    with p.open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def valid(txt):
    if not txt:
        return False
    d=Decimal(txt)
    return d.is_finite() and d>0

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score",type=Path)
    parser.add_argument("member",type=Path)
    parser.add_argument("blind_design",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    design={(r["variant"],r["pair_id"]):r for r in table(args.blind_design)}
    identities={(r["candidate_symbol"],r["candidate_id"]):r for r in table(args.member)}
    lookup={}
    with localcontext() as ctx:
        ctx.prec=50
        for ident, row in identities.items():
            pr_a,pr_b=row["PR_source"],row["PR_target"]
            t_a,t_b=row["RNA_source_mean_TPM"],row["RNA_target_TPM"]
            if ident[1] in ("","unresolved","NA") or not(valid(pr_a) and valid(pr_b)):
                lookup[ident]=None
            elif not t_a or not t_b or Decimal(t_b)<=Decimal(t_a):
                lookup[ident]=None
            else:
                a,b=Decimal(pr_a),Decimal(pr_b)
                lookup[ident]=max(a/b,b/a)
        recomputed={}
        emitted={}
        for row in table(args.score):
            key=(row["variant"],row["pair_id"])
            p=design[key]
            a=(p["bhlh_symbol"],p["bhlh_id"])
            b=(p["partner_symbol"],p["partner_id"])
            lhs,rhs=lookup[a],lookup[b]
            projected=lhs*rhs if lhs is not None and rhs is not None and a[1]!=b[1] else None
            observed=Decimal(row["score"]) if row["score"] else None
            assert key not in recomputed and (observed is None)==(projected is None),key
            assert row["eligible"]==str(int(projected is not None)),key
            if observed is not None:
                assert abs(observed-projected)<=Decimal("3e-15")*projected,(key,observed,projected)
            recomputed[key]=projected
            emitted[key]=observed
    assert len(recomputed)==len(design)==9*598 and set(recomputed)==set(design)
    keys=sorted((k for k,r in recomputed.items() if k[0]=="main" and r is not None),
                key=lambda k:sha256(k[1].encode("ascii")).digest())
    checks={"full_slots_verified":len(recomputed),"main_finite_count":len(keys),
            "independent_product":"50-digit Decimal ratio, RNA comparison and two-member gate",
            "blind_hash_selected_main":[{"pair_id":k[1], "observed":str(emitted[k]),
                                         "recomputed":str(recomputed[k])} for k in keys[:12]]}
    args.output.write_text(json.dumps(checks,indent=2)+"\n")

if __name__=="__main__":
    main()
