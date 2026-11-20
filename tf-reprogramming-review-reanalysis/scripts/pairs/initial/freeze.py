"""Aggregate only blinded output artifacts and seal full design rankings."""
import argparse
from collections import defaultdict
from decimal import Decimal, ROUND_HALF_EVEN
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
from math import isfinite

from shared import js,json_write,now,sha,tsv,write,verify_lock

ENDPOINTS=("GSE53951","GSE114809")
RAW_FILES=("scores/core/core.tsv","scores/imcell.tsv","scores/native_modules.tsv","scores/preflight_modules.tsv")
INTEGRATION={"integration","integration_expression","integration_response"}
FALLBACK={"imcell_p010","sfa"}
BASE="expression"
STATUS=["已跑且可解释","已跑但有条件","退化无区分","仅集合/证书","输入不足","实现失败"]


def fraction(row):
    if row.get("score","") in {"","NA","None","NaN"}:
        return None
    if row.get("score_num") and row.get("score_den"):
        return Fraction(int(row["score_num"]),int(row["score_den"]))
    decimal=Decimal(row["score"])
    if not decimal.is_finite():
        raise ValueError(f"Nonfinite score should be marked NA: {row}")
    return Fraction(decimal.quantize(Decimal("0.000000000001"),rounding=ROUND_HALF_EVEN))


def order(numbers,ties):
    keys=[pid for pid,n in numbers.items() if n is not None]
    keys.sort(key=lambda pid:(-numbers[pid],ties[pid],pid))
    ranks={}
    i=0
    while i<len(keys):
        j=i+1
        while j<len(keys) and numbers[keys[j]]==numbers[keys[i]]:
            j+=1
        for k in range(i,j):
            ranks[keys[k]] = (k+1,i+1,j)
        i=j
    return ranks


def pct(numbers,ties):
    ranks=order(numbers,ties)
    n=len(ranks)
    return {pid:Fraction(2*n+1-lo-hi,2*n) for pid,(_,lo,hi) in ranks.items()}


def main(root):
    lock=root/"SCORES_LOCK.sha256"
    if lock.exists():
        verify_lock(root)
        return
    frozen=js(root/"locks/REGISTRATION.json")
    assert all(sha(root/name)==digest for name,digest in frozen["policy"].items())
    assert sha(root/"scores/core/all_pair_scores.tsv")==frozen["old_scores_sha256"]
    by_method={r["method"]:r for r in tsv(root/"METHOD_REGISTRY.tsv")}
    design=tsv(root/"inputs/CLEAN_PAIRS.tsv")
    salt=js(root/"inputs/SCORING_POLICY.json")["tie_salt"]
    variants=sorted({r["variant"] for r in design})
    pairs={(r["variant"],r["pair_id"]):r for r in design}
    ties={(v,p):sha256(f'{salt}|{r["bhlh_id"]}|{r["partner_id"]}'.encode()).hexdigest() for (v,p),r in pairs.items()}
    scores=defaultdict(dict)
    sources=[]
    for file in RAW_FILES:
        p=root/file
        if not p.is_file():
            continue
        sources.append(file)
        for row in tsv(p):
            key=row["endpoint"],row["variant"],row["method"]
            pid=row["pair_id"]
            assert key[0] in ENDPOINTS and key[2] in by_method and (key[1],pid) in pairs,(key,pid)
            if pid in scores[key]:
                raise ValueError(f"Duplicate score row {file}/{key}/{pid}")
            details=row.get("details","")
            if key[2].startswith("imcell_"):
                diag_fields=("activation_sum","repression_sum","sum_square_total","mc_se","rep1","rep2","rep3",
                             "n_active_goals","n_repress_goals","bhlh_single","partner_single","single_sum","joint_minus_single_sum")
                details=json.dumps({"source":details,**{c:row[c] for c in diag_fields if row.get(c)}},separators=(",",":"))
            scores[key][pid]=(fraction(row),row.get("reason",""),details,row.get("score",""))
    eligible_imcell=any(row[0] is not None for (ep,v,m),cache in scores.items() if m=="imcell_p010" for row in cache.values())
    selected="imcell_p010" if eligible_imcell else "sfa"
    statuses=[]
    for ep in ENDPOINTS:
        for v in variants:
            A=scores[ep,v,BASE]
            B=scores[ep,v,selected]
            common={p for p in A.keys() & B.keys() if A[p][0] is not None and B[p][0] is not None}
            tk={p:ties[v,p] for p in common}
            P=pct({p:A[p][0] for p in common},tk)
            Q=pct({p:B[p][0] for p in common},tk)
            for pid in common:
                details=json.dumps({"response_module":selected,"expression_pct":str(P[pid]),"response_pct":str(Q[pid]),
                                    "expression_score":A[pid][3],"response_score":B[pid][3]},sort_keys=True,separators=(",",":"))
                values=(("integration",(P[pid]+Q[pid])/2),("integration_expression",P[pid]),("integration_response",Q[pid]))
                for method,value in values:
                    scores[ep,v,method][pid]=(value,"",details,str(float(value)))
    paths=["PROTOCOL.yaml","METHOD_REGISTRY.tsv","INPUTS_BY_METHOD.tsv","FIGURE_POLICY.yaml","locks/REGISTRATION.json",
           "validation/BLIND_DESIGN_AUDIT.tsv","validation/FVS_SELF_LOOPS.tsv","scores/core/all_pair_scores.tsv",
           "scores/core/program.tsv","scores/core/CORE_AUDIT.json","scores/core/READ_ACCESS.json",*sources]
    paths += [f"inputs/{p.name}" for p in (root/"inputs").glob("*") if p.is_file()]
    audit_files=("validation/VIPER_STATUS.tsv","validation/TFCOmb_STATUS.tsv","external/preflight/METHOD_PREFLIGHT_STATUS.tsv")
    extra_status={}
    for p in root.rglob("*STATUS.tsv"):
        if p.name in {"METHOD_STATUS.tsv"} or ".git" in p.parts:
            continue
        for row in tsv(p):
            if row.get("method") in by_method:
                extra_status[row["method"]]=row
                if p.relative_to(root).as_posix() not in paths:
                    paths.append(p.relative_to(root).as_posix())
    entries=[]
    for method,meta in by_method.items():
        n_all=sum(n is not None for (e,v,m),cache in scores.items() if m==method for n,*_ in cache.values())
        has_row=any(m==method and cache for (e,v,m),cache in scores.items())
        unique_scores={n for (e,v,m),cache in scores.items() if m==method for n,*_ in cache.values() if n is not None}
        proof=extra_status.get(method,{})
        hint=proof.get("reason","") or proof.get("details","") or proof.get("output_and_rule","")
        state="仅集合/证书" if method=="fvs_certificate" else "已跑且可解释" if n_all and method in {"expression","expr_product","network_add","network_union"} else "已跑但有条件" if n_all else "输入不足"
        if n_all and len(unique_scores)==1:
            state="退化无区分"
        if not n_all and proof.get("execution_status","") in {"execution_failed","implementation_failed"}:
            state="实现失败"
        if method=="fvs_certificate":
            hint="16 个自环必须逐一覆盖，任意两节点都不是全图的完整 FVS。属于结构证书，不是实验的负标签。"
        if not hint and not n_all and has_row:
            hint="无标签资格预检未产生可逐对排名结果，查看逐对缺失原因和原生日志。"
        if not hint and not n_all and not has_row and method!="fvs_certificate":
            hint="查看官方数据或实现预检记录。"
        statuses.append({"method":method,"family":meta["family"],"state":state,"scored_rows":n_all,
                         "implementation":proof.get("implementation","") or meta["implementation"],"source":meta["source"],
                         "notes":hint,"detail_status":proof.get("execution_status","")})
        for ep in ENDPOINTS:
            for v in variants:
                cache=scores[ep,v,method]
                numbers={p:r[0] for p,r in cache.items()}
                ranks=order(numbers,{p:ties[v,p] for p in numbers})
                for pid in sorted(p for vv,p in pairs if vv==v):
                    pair=pairs[v,pid]
                    n,reason,details,original=cache.get(pid,(None,"unavailable_registered_module","",""))
                    if method in INTEGRATION and pid not in ranks:
                        reason="component_intersection_missing"
                    if method=="fvs_certificate":
                        reason="two_members_cannot_hit_16_self_loops"
                    ordinal,lo,hi=ranks.get(pid,("","",""))
                    entries.append({**pair,"method":method,"endpoint":ep,"operation":meta["operation"],"source_direction":meta["direction"].strip(),
                                    "training_exposure":"independent_external_or_not_applicable","scoreable":int(n is not None),
                                    "reason":reason,"score":format(float(n),".17g") if n is not None else "",
                                    "raw_score":original,"score_num":n.numerator if n is not None else "","score_den":n.denominator if n is not None else "",
                                    "rank":ordinal,"rank_min":lo,"rank_max":hi,"midrank":(lo+hi)/2 if lo else "",
                                    "N":len(ranks),"tie_key":ties[v,pid],"details":details})
    write(root/"scores/BLIND_RANKS.tsv",entries)
    write(root/"METHOD_STATUS.tsv",statuses)
    json_write(root/"locks/FREEZE_META.json",{"sealed_utc":now(),"grid_rows":len(entries),"response_for_integration":selected,"registered_methods":len(by_method),"score_artifacts":sources})
    paths.extend(("scores/BLIND_RANKS.tsv","METHOD_STATUS.tsv","locks/FREEZE_META.json","code/freeze.py","code/core_score.py","code/historical_score.py",
                  "code/run_core.sh","external/sfa/sfa/utils.py","external/sfa/sfa/algorithms/sp.py"))
    lock.write_text("".join(f"{sha(root/file)}  {file}\n" for file in sorted(set(paths))),encoding="utf-8")
    print(f"Locked {len(entries)} rows. Integration uses {selected}.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path)
    main(p.parse_args().root)
