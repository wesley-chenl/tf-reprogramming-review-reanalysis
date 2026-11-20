"""Join frozen scores with verified outcomes and evaluate ranks on named panels."""
import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, groupby
from math import comb, fsum, sqrt
from pathlib import Path
import shutil

import numpy as np

from shared import js, json_write, now, sha, tsv, verify_lock, write

MAIN = ("expression","two_union","network_union","imcell_p010","boolean","integration")
PRIMARY="GSE53951"
SEED=20261119
BUDGETS=(10,20)
DRAWS=10000
UTILITY=("binary","gated")
NA={"","NA","unknown"}
OUTCOME_FILES=("design_pairs.tsv","outcomes_private.tsv","mapping_conflicts.tsv","constructs.tsv")


def key(r):
    return r["score_num"],r["score_den"]


def condition(r):
    return r["allowed"] and r["scoreable"]=="1" and r["binary"] is not None


def safe_number(s):
    return None if s in NA else float(s)


def seed(x):
    return np.random.default_rng(np.random.SeedSequence([SEED,int.from_bytes(sha256(repr(x).encode()).digest()[:8],"big")]))


def pct(v):
    return np.quantile(v,[.025,.975],method="inverted_cdf").tolist()


def ordered(rows):
    return sorted(rows,key=lambda r:int(r["rank"]))


def group_rank(rows):
    out={}
    pos=1
    for _, block in groupby(rows,key=key):
        group=list(block)
        hi=pos+len(group)-1
        mid=(pos+hi)/2
        value=(mid-1)/(len(rows)-1) if len(rows)>1 else None
        for offset,r in enumerate(group):
            out[r["pair_id"]]=(pos+offset,pos,hi,mid,value)
        pos=hi+1
    return out


def binary_metrics(rows):
    grouped=[(len(z),sum(r["binary"]==1 for r in z)) for _,g in groupby(rows,key=key) if (z:=list(g))]
    n=len(rows)
    k=sum(p for _,p in grouped)
    l=n-k
    if not k or not l:
        return n,k,None,None,None,tuple(c for c,p in grouped)
    correct=0
    ap=0.
    expected=0.
    cp,cn=0,0
    for count,p in grouped:
        q=count-p
        cp+=p
        cn+=q
        correct+=p*(l-cn+0.5*q)
        ap+=p*cp/(cp+cn)/k
        expected+=(count/n+count*(cp+cn-1)*(k-1)/(n*(n-1)))/(cp+cn)
    return n,k,correct/(k*l),ap,expected,tuple(c for c,p in grouped)


@lru_cache(maxsize=256)
def permutation_reference(sizes,k):
    n=sum(sizes)
    if k==0 or k==n:
        return None,None,None,None
    g=seed(("rank-place",sizes,k))
    p=g.multivariate_hypergeometric(sizes,k,size=DRAWS,method="marginals")
    groups=np.asarray(sizes)
    cp=p.cumsum(axis=1)
    cn=groups.cumsum() - cp
    score_auc=(p*((n-k)-cn+0.5*(groups-p))).sum(axis=1)/(k*(n-k))
    score_ap=(p*cp/groups.cumsum()).sum(axis=1)/k
    return *pct(score_auc),*pct(score_ap)


def stratify(rows):
    groups=defaultdict(list)
    for row in rows:
        groups[row["bhlh_id"]].append(row)
    info=[]
    for bhlh, subgroup in sorted(groups.items()):
        n,k,auc,*_=binary_metrics(subgroup)
        if auc is not None:
            info.append({"bhlh_id":bhlh,"symbol":subgroup[0]["bhlh_symbol"],"N":n,"K":k,"positive_negative_pairs":k*(n-k),"concordance":auc})
    contrasts=sum(r["positive_negative_pairs"] for r in info)
    return {"bhlh_informative_n":len(info),"within_pair_comparisons":contrasts,
            "within_equal_bhlh":fsum(r["concordance"] for r in info)/len(info) if info else None,
            "within_weighted_bhlh":fsum(r["concordance"]*r["positive_negative_pairs"] for r in info)/contrasts if info else None},info


def metric(rows, detail=False):
    n,k,auc,ap,expected,shape=binary_metrics(rows)
    out={"N":n,"K":k,"pos_fraction":k/n if n else None,"AUROC":auc,"AP":ap,"null_AUROC_expected":0.5 if auc is not None else None,
         "null_AP_expected":expected,"distinct_scores":len(shape),"ties_pairs_fraction":sum(x*(x-1) for x in shape)/(n*(n-1)) if n>1 else None}
    if detail and auc is not None:
        out.update(dict(zip(("rand_AUC_q025","rand_AUC_q975","rand_AP_q025","rand_AP_q975"),permutation_reference(shape,k))))
    strat,sub=stratify(rows)
    out.update(strat)
    return out,sub


def hypergeom(n,k,b):
    return np.array([(comb(k,i)*comb(n-k,b-i))/comb(n,b) if i<=k and 0<=b-i<=n-k else 0 for i in range(b+1)])


def qdist(p, scale=1):
    c=np.cumsum(p)
    return [int(np.searchsorted(c,b))/scale for b in (.025,.975)]


def take(n,b,rng):
    if b==n:
        return np.broadcast_to(np.arange(n),(DRAWS,n))
    if n<3*b:
        return np.argpartition(rng.random((DRAWS,n)),b-1,axis=1)[:,:b]
    out=np.empty((DRAWS,b),dtype=np.int32)
    for j in range(b):
        x=rng.integers(0,n,size=DRAWS,dtype=np.int32)
        conflict=np.any(x[:,None]==out[:,:j],axis=1)
        while np.any(conflict):
            x[conflict]=rng.integers(0,n,size=conflict.sum(),dtype=np.int32)
            conflict=np.any(x[:,None]==out[:,:j],axis=1)
        out[:,j]=x
    return out


def num_summary(rows,label,ranks):
    subset=[ranks[r["pair_id"]] for r in rows if r["binary"]==label]
    if not subset:
        return None,None,None,None
    rr=[r[3] for r in subset]
    low,mid,hi=np.quantile(rr,[.25,.5,.75]).tolist()
    norm=[r[4] for r in subset if r[4] is not None]
    return low,mid,hi,float(np.median(norm)) if norm else None


@lru_cache(maxsize=192)
def gated_uniform(ys,b):
    series=np.asarray(ys)
    choices=series[take(len(ys),b,seed(("gated",ys,b)))]
    return (*pct(choices.mean(axis=1)),*pct(series.max()-choices.max(axis=1)))


def budget(rows, B, utilities=True):
    n=len(rows)
    b=min(B,n)
    k=sum(r["binary"]==1 for r in rows)
    out={"B":B,"B_eff":b,"N":n,"K":k,
         "budget_status":"none" if n==0 else "truncated" if b<B else "exhaustive" if n==b else "complete"}
    if not b:
        return out,[]
    picked=rows[:b]
    h=sum(r["binary"]==1 for r in picked)
    eq=[r for r in rows if key(r)==key(picked[-1])]
    before=next(i for i,r in enumerate(rows) if key(r)==key(picked[-1]))
    slots=b-before
    blockh=sum(r["binary"]==1 for r in eq)
    fixed_h=sum(r["binary"]==1 for r in rows[:before])
    qlo,qhi=qdist(hypergeom(n,k,b))
    out.update(hits=h,negatives=b-h,uniform_expected_hits=b*k/n,uniform_hit_q025=qlo,uniform_hit_q975=qhi,
               boundary_group_n=len(eq),boundary_slots=slots,boundary_fraction=slots/len(eq),
               tie_expected_hits=fixed_h+slots*blockh/len(eq),
               tie_min_hits=fixed_h+max(0,slots-(len(eq)-blockh)),tie_max_hits=fixed_h+min(slots,blockh),
               mean_binary=h/b,best_regret_binary=int(k>0 and h==0),
               uniform_p_zero_binary_regret=1 if not k or n-k<b else 1-comb(n-k,b)/comb(n,b))
    quotas=Counter(r["bhlh_id"] for r in picked)
    by=defaultdict(list)
    for r in rows:
        by[r["bhlh_id"]].append(r)
    expectation=fsum(quotas[g]*sum(r["binary"] for r in by[g])/len(by[g]) for g in quotas)
    reference=np.array([1.])
    for g,c in sorted(quotas.items()):
        reference=np.convolve(reference,hypergeom(len(by[g]),sum(r["binary"]==1 for r in by[g]),c))
    q1,q2=qdist(reference)
    out.update(fixed_bhlh_expected_hits=expectation,fixed_bhlh_hit_q025=q1,fixed_bhlh_hit_q975=q2,
               fixed_bhlh_informative_rows=sum(0<sum(r["binary"] for r in by[g])<len(by[g]) and c<len(by[g]) for g,c in quotas.items()),
               fixed_bhlh_delta_hits=h-expectation,bhlh_quotas=",".join(f"{g}:{c}" for g,c in sorted(quotas.items())))
    if utilities and all(r["gated"] is not None for r in rows):
        ys=tuple(sorted(r["gated"] for r in rows))
        vals=[r["gated"] for r in picked]
        lo,hi,rlo,rhi=gated_uniform(ys,b)
        g_rng=seed(("partner-gated",tuple(sorted(quotas.items())),tuple(sorted(r["pair_id"] for r in rows))))
        sums=np.zeros(DRAWS)
        pbest=np.zeros(DRAWS)
        q_expect=0.
        for g,c in sorted(quotas.items()):
            counts=np.asarray([r["gated"] for r in by[g]])
            draws=counts[take(len(counts),c,g_rng)]
            sums+=draws.sum(axis=1)
            pbest=np.maximum(pbest,draws.max(axis=1))
            q_expect+=c*counts.mean()
        pq=pct(sums/b)
        pr=pct(max(ys)-pbest)
        out.update(mean_gated_pct=fsum(vals)/b,best_regret_gated_pp=max(ys)-max(vals),
                   uniform_mean_gated_pct=fsum(ys)/n,uniform_gated_q025=lo,uniform_gated_q975=hi,
                   uniform_gated_regret_q025=rlo,uniform_gated_regret_q975=rhi,
                   fixed_bhlh_mean_gated_pct=q_expect/b,fixed_bhlh_gated_q025=pq[0],fixed_bhlh_gated_q975=pq[1],
                   fixed_bhlh_gated_regret_q025=pr[0],fixed_bhlh_gated_regret_q975=pr[1])
    return out,picked


def copy_labels(root,old):
    dest=root/"evaluation/label_source"
    dest.mkdir(parents=True,exist_ok=True)
    original=js(old/"locks/LABEL_LOCK.json")
    for name in OUTCOME_FILES:
        src=old/"sources_private"/name
        check=original["files"][f"sources_private/{name}"]
        assert sha(src)==check
        p=dest/name
        if p.is_file() and sha(p)!=check:
            raise ValueError("Prepared evaluation label copy changed")
        if not p.exists():
            shutil.copyfile(src,p)
    json_write(dest/"LABEL_LINEAGE.json",{"checked_utc":now(),"old_label_lock":sha(old/"locks/LABEL_LOCK.json"),
                                          "four_label_files":{name:sha(dest/name) for name in OUTCOME_FILES}})
    return dest


def main(root,old):
    verify_lock(root)
    stage=root/"evaluation"
    start=stage/"JOIN_STARTED.json"
    lock_digest=sha(root/"SCORES_LOCK.sha256")
    if not start.exists():
        json_write(start,{"at":now(),"score_lock_sha256":lock_digest,"sealed_utc":js(root/"locks/FREEZE_META.json")["sealed_utc"]})
    else:
        assert js(start)["score_lock_sha256"]==lock_digest
    dest=copy_labels(root,old)
    d={(r["variant"],r["pair_id"]):r for r in tsv(dest/"design_pairs.tsv")}
    obs={r["pair_slot_id"]:r for r in tsv(dest/"outcomes_private.tsv") if r["source_version"]=="nature2018"}
    assert len(obs)==598 and Counter(r["binary_call"] for r in obs.values())=={"1":76,"0":522}
    grouped=defaultdict(list)
    data=tsv(root/"scores/BLIND_RANKS.tsv")
    for r in data:
        local=d[r["variant"],r["pair_id"]]
        for col in ("bhlh_id","partner_id","bhlh_symbol","partner_symbol"):
            if local[col]!=r[col]:
                raise ValueError(f"Score/label design mismatch {r['pair_id']} {col}")
        y=obs[local["pair_slot_id"]]
        r.update(slot=local["pair_slot_id"],bhlh_construct=local["bhlh_axis"],partner_construct=local["partner_axis"],
                 allowed=local["eligible"]=="1",design_exclusion=local["exclusion_reason"],
                 binary=safe_number(y["binary_call"]),gated=safe_number(y["q_gated"]),raw_S1=safe_number(y["raw_s1_pct"]),
                 fig_printed=safe_number(y["fig1b_pct_printed"]),source_locator=y["source_cell_or_bbox"],
                 raw_locator=y["raw_s1_cell"],call_basis=y["call_basis"])
        grouped[r["endpoint"],r["variant"],r["method"]].append(r)
    native={k:ordered([r for r in rows if condition(r)]) for k, rows in grouped.items()}
    response=js(root/"locks/FREEZE_META.json")["response_for_integration"]
    reps=tuple(response if m=="imcell_p010" else m for m in MAIN)
    scopes={}
    common_members={}
    for ep,v in sorted({(e,v) for e,v,m in grouped}):
        active=[m for m in reps if native[ep,v,m]]
        common_members[ep,v]=set.intersection(*[{r["pair_id"] for r in native[ep,v,m]} for m in active]) if active else set()
        for m in [method for e,version,method in grouped if (e,version)==(ep,v)]:
            rows=native[ep,v,m]
            scopes[ep,v,m,"native"]=rows
            scopes[ep,v,m,"common"]=[r for r in rows if r["pair_id"] in common_members[ep,v]]
    rankings={}
    summaries,budgets,top,strata,flows=[],[],[],[],[]
    for (ep,v,m,scope),rows in scopes.items():
        full=grouped[ep,v,m]
        distinct_relation = "native" if scope=="native" else "same_common_panel" if len(rows)==len(common_members[ep,v]) and rows else "restricted_common_intersection"
        identity={"endpoint":ep,"variant":v,"method":m,"panel":scope,"panel_relation":distinct_relation}
        ranks=group_rank(rows)
        rankings[ep,v,m,scope]=ranks
        met,sub=metric(rows,v=="main")
        sizes={}
        for label,tag in ((1,"positive"),(0,"negative")):
            q1,med,q3,norm=num_summary(rows,label,ranks)
            count=sum(r["binary"]==label for r in rows)
            denom=76 if label else 522
            allowed_n=sum(r["allowed"] and r["binary"]==label for r in full)
            sizes.update({f"{tag}_eligible":count,f"{tag}_all":denom,f"{tag}_coverage":count/denom,
                          f"{tag}_resolved_design_n":allowed_n,f"{tag}_resolved_coverage":count/allowed_n if allowed_n else None,
                          f"{tag}_rank_q1":q1,f"{tag}_rank_median":med,f"{tag}_rank_q3":q3,f"{tag}_normalized_median":norm})
        summaries.append({**identity,**met,**sizes})
        flows.append({**identity,"physical_design_n":len(full),"blind_ranked_n":sum(r["scoreable"]=="1" for r in full),
                      "identity_allowed_n":sum(r["allowed"] for r in full),"known_binary_n":sum(r["binary"] is not None for r in full),
                      "native_eval_N":len(native[ep,v,m]),"panel_eval_N":len(rows),"panel_positive":met["K"]})
        strata.extend({**identity,**r} for r in sub)
        for b in BUDGETS:
            stats,chosen=budget(rows,b,utilities=True)
            budgets.append({**identity,**stats})
            top.extend({**identity,"B":b,"rank":i,"slot":r["slot"],"pair_id":r["pair_id"],"bhlh":r["bhlh_symbol"],
                            "partner":r["partner_symbol"],"binary":r["binary"],"gated_pct":r["gated"],"score":r["score"],
                            "source_locator":r["source_locator"],"details":r["details"]} for i,r in enumerate(chosen,1))

    metrics_by_key={(r["endpoint"],r["variant"],r["method"],r["panel"]):r for r in summaries}
    pairs_test,pair_variation=[],[]
    for ep in sorted({e for e,v,m in grouped}):
        tested=[m for e,v,m in grouped if e==ep and v=="main" and native[e,v,m]]
        for a,b in combinations(tested,2):
            x=native[ep,"main",a]
            y=native[ep,"main",b]
            xy={r["pair_id"] for r in x} & {r["pair_id"] for r in y}
            if not xy:
                pairs_test.append({"endpoint":ep,"method_a":a,"method_b":b,"N":0})
                continue
            pools=([r for r in x if r["pair_id"] in xy],[r for r in y if r["pair_id"] in xy])
            ar,br=[group_rank(p) for p in pools]
            pids=sorted(xy)
            rv0=[ar[p][3] for p in pids]
            rv1=[br[p][3] for p in pids]
            rho=float(np.corrcoef(rv0,rv1)[0,1]) if len(set(rv0))>1 and len(set(rv1))>1 else None
            ax,bx=[metric(p,False)[0] for p in pools]
            descr={"endpoint":ep,"method_a":a,"method_b":b,"N":len(xy),"K":ax["K"],"spearman_midrank":rho,
                   "AUC_a":ax["AUROC"],"AUC_b":bx["AUROC"],"AP_a":ax["AP"],"AP_b":bx["AP"],
                   "within_a":ax["within_weighted_bhlh"],"within_b":bx["within_weighted_bhlh"],
                   "tie_fraction_a":ax["ties_pairs_fraction"],"tie_fraction_b":bx["ties_pairs_fraction"]}
            for k in (10,20):
                left=budget(pools[0],k,utilities=False)[0]
                right=budget(pools[1],k,utilities=False)[0]
                descr.update({f"top{k}_overlap":len({r["pair_id"] for r in pools[0][:k]} & {r["pair_id"] for r in pools[1][:k]}),
                              f"top{k}_a_hits":left["hits"],f"top{k}_b_hits":right["hits"],
                              f"top{k}_random_expected_hits":left["uniform_expected_hits"]})
            pairs_test.append(descr)
            if len(xy)>1:
                lookup=[{r["pair_id"]:r for r in p} for p in pools]
                for p in sorted(pids,key=lambda p:(-abs(ar[p][4]-br[p][4]),p))[:5]:
                    ra,rb=(table[p] for table in lookup)
                    pair_variation.append({**{k:descr[k] for k in ("endpoint","method_a","method_b","N","K")},
                                           "slot":ra["slot"],"pair_id":p,"binary":ra["binary"],"gated_pct":ra["gated"],
                                           "bhlh":ra["bhlh_symbol"],"partner":ra["partner_symbol"],"rank_a":ar[p][3],"rank_b":br[p][3],
                                           "normalized_delta":ar[p][4]-br[p][4],"source_locator":ra["source_locator"],
                                           "decomposition_a":ra["details"],"decomposition_b":rb["details"]})
    long,wide_by_id,examples=[],{},[]
    for (ep,v,m),rows in grouped.items():
        if v!="main":
            continue
        rank_a=rankings[ep,v,m,"native"]
        rank_b=rankings[ep,v,m,"common"]
        for r in rows:
            p=r["pair_id"]
            a=rank_a.get(p)
            b=rank_b.get(p)
            base={field:r[field] for field in ("slot","pair_id","bhlh_construct","partner_construct","bhlh_id","partner_id","bhlh_symbol",
                   "partner_symbol","binary","gated","raw_S1","fig_printed","source_locator","raw_locator","call_basis")}
            out={"endpoint":ep,"variant":v,"method":m,**base,"score":r["score"],"scoreable_blind":r["scoreable"],
                 "native_N":len(rank_a),"common_N":len(rank_b),"native_eligible":int(a is not None),"common_eligible":int(b is not None),
                 "reason":r["design_exclusion"] or r["reason"],"details":r["details"],
                 "training_exposure":r["training_exposure"]}
            for prefix,rank in (("native",a),("common",b)):
                for name,value in zip(("rank","lo","hi","mid","normalized"),rank if rank else (None,)*5):
                    out[prefix+"_"+name]=value
            long.append(out)
            record=wide_by_id.setdefault(r["slot"],base.copy())
            record[f"{ep}__{m}"]=(f"{a[1]}-{a[2]} / {len(rank_a)}" if a[1]!=a[2] else f"{a[0]} / {len(rank_a)}") if a else "NA "+out["reason"]
        native_pool=native[ep,v,m]
        for name,label,seq in (("highest_negative",0,native_pool),("lowest_positive",1,list(reversed(native_pool)))):
            for r in [r for r in seq if r["binary"]==label][:5]:
                p=r["pair_id"]
                examples.append({"endpoint":ep,"method":m,"category":name,"slot":r["slot"],"bhlh":r["bhlh_symbol"],
                                 "partner":r["partner_symbol"],"rank":rank_a[p][0],"mid":rank_a[p][3],"N":len(native_pool),
                                 "binary":r["binary"],"gated_pct":r["gated"],"source_locator":r["source_locator"],"details":r["details"]})
    write(stage/"PAIR_RANKS_LONG.tsv",long)
    write(stage/"PAIR_RANKS_WIDE.tsv",sorted(wide_by_id.values(),key=lambda r:r["slot"]))
    for tag,label in (("POSITIVE",1),("NEGATIVE",0)):
        write(stage/f"SCREEN_{tag}_ALL.tsv",[r for r in long if r["binary"]==label])
    write(stage/"RANK_METRICS.tsv",summaries)
    write(stage/"PANEL_FLOW.tsv",flows)
    write(stage/"BUDGET_UTILITY.tsv",budgets)
    write(stage/"TOP_B_PAIRS.tsv",top)
    write(stage/"PAIRWISE_SCOPE.tsv",pairs_test)
    write(stage/"PARTNER_STRATA.tsv",strata)
    write(stage/"DISCORDANT_PAIRS.tsv",pair_variation)
    write(stage/"RANK_EXAMPLES.tsv",examples)
    json_write(stage/"ANALYSIS_META.json",{"evaluated_utc":now(),"score_lock":lock_digest,"ranks_long":len(long),
             "physical_slots":len(wide_by_id),"metrics_rows":len(summaries),"budget_rows":len(budgets),"story_representatives":list(reps),
             "seed":SEED,"draws":DRAWS,"labels_n":len(obs),
             "common_sizes":{f"{ep}__{v}":len(pairs) for (ep,v),pairs in common_members.items()}})
    print(f"{len(long)} main slot-method rows, {len(summaries)} metric rows, {len(pairs_test)} pairwise panels")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path)
    p.add_argument("old",type=Path)
    a=p.parse_args()
    main(a.root,a.old)
