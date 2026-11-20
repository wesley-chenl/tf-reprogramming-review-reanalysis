"""Rebuild the three figure panels from the two frozen slot-level rank files."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
import csv
import gzip
from fractions import Fraction as F
from hashlib import sha256, file_digest
from itertools import groupby
import json
from pathlib import Path
from statistics import median

A = ("expression", "network_union", "two_union", "imcell_p010", "sfa", "integration")
B = ("expression_GSE43916", "taiji_pair")
C = ("expression", "tfcomb")
LABELS = dict(expression="Expression mean", network_union="Direct signed union", two_union="Two-step signed union",
              imcell_p010="IMCell dual seed", sfa="SFA", integration="Expression + IMCell",
              expression_GSE43916="Expression mean", taiji_pair="Taiji-reprogram",
              tfcomb="TFcomb scoring module (bulk adaptation)")
FIELDS = ("pair_id","bhlh_id","partner_id","bhlh_symbol","partner_symbol","bhlh_construct","partner_construct")
REF_COLS = {"old": ("common_rank","common_lo","common_hi"),
            "followup": ("old_reference_rank","old_reference_lo","old_reference_hi")}

def digest(path):
    with path.open("rb") as f:
        return file_digest(f,"sha256").hexdigest()

def rows(path):
    open_input = gzip.open if path.suffix == ".gz" else open
    with open_input(path, "rt", newline="",encoding="utf-8-sig") as f:
        yield from csv.DictReader(f,delimiter="\t")

def write(path, items):
    assert items, path
    fields = list(dict.fromkeys(k for row in items for k in row))
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(items)

def numeric(v):
    return int(F(v)) if v else None

def load(path, cohort):
    groups=defaultdict(dict)
    for line,r in enumerate(rows(path),2):
        if r["variant"]!="main":
            continue
        label_name="binary" if cohort=="old" else "screen_binary"
        label=F(r[label_name])
        assert label in (0,1), (cohort,line,label)
        n_rank=numeric(r["native_rank"])
        if cohort=="old":
            assert (r["native_eligible"]=="1")==bool(n_rank),(line,r["slot"])
            assert (r["common_eligible"]=="1")==bool(r["common_rank"]),(line,r["slot"])
            assert not n_rank or (r["scoreable_blind"]=="1" and r["score"]), (line,r["slot"])
        else:
            assert not n_rank or (r["eligible_design"]=="1" and r["exact_score"]), (line,r["slot"])
        a,b,c=REF_COLS[cohort]
        rec={k:r[k] for k in FIELDS}
        rec.update(slot=r["slot"],screen_positive=int(label),source_line=line,source_locator=r["source_locator"],
                   source_native_N=int(r["native_N"]),
                   native_rank=n_rank,native_lo=numeric(r["native_lo"]),native_hi=numeric(r["native_hi"]),
                   reference_rank=numeric(r[a]),reference_lo=numeric(r[b]),reference_hi=numeric(r[c]),
                   reference_full_N=r["common_N"] if cohort=="old" else r["old_reference_N"],
                   reference_intersection_N=r.get("old_reference_intersection_N",""),
                   design_eligible=r.get("eligible_design",""),
                   score=r["score"] if cohort=="old" else r["blind_score"],
                   exact_score=r.get("exact_score",""),
                   exclusion=r["reason"] if cohort=="old" else r["not_ranked_reason"],
                   q_gated=r["gated"] if cohort=="old" else r["gated_pct"],
                   raw_yield_pct=r["raw_S1"] if cohort=="old" else r["raw_s1_pct"],
                   measurement_status=r.get("measurement_status",""), call_basis=r["call_basis"])
        group=groups[r["endpoint"],r["method"]]
        assert rec["slot"] not in group,(cohort,line,rec["slot"])
        group[rec["slot"]]=rec
    return groups

def check_slots(old,new):
    anchor=old["GSE53951","expression"]
    assert Counter(r["screen_positive"] for r in anchor.values())=={0:522,1:76}
    design=new["GSE43916","expression_GSE43916"]
    assert Counter(r["screen_positive"] for r in design.values() if r["design_eligible"]=="1")=={0:475,1:71}
    for cohort,groups in (("old",old),("followup",new)):
        for key,group in groups.items():
            diff=group.keys()^anchor.keys()
            assert not diff,(cohort,key,next(iter(diff),None))
            for p,r in group.items():
                for f in (*FIELDS,"screen_positive","q_gated","raw_yield_pct"):
                    first,second=anchor[p][f],r[f]
                    if f in {"q_gated","raw_yield_pct"} and first and second:
                        first,second=F(first),F(second)
                    assert first==second,(cohort,key,p,f,first,second)
                if cohort=="followup":
                    assert r["design_eligible"]==design[p]["design_eligible"],(key,p)
                assert not r["native_rank"] or design[p]["design_eligible"]=="1",(cohort,key,p)
            nn=sum(bool(x["native_rank"]) for x in group.values())
            assert all(x["source_native_N"]==nn for x in group.values()),(cohort,key)

def valid(group,scope):
    return {p for p,r in group.items() if r[scope+"_rank"]}

def rerank(group, slots, scope="native"):
    ordered=sorted((group[p] for p in slots),key=lambda r:r[scope+"_rank"])
    n=len(ordered)
    rebuilt=[]
    for _,tie in groupby(ordered,key=lambda r:(r[scope+"_lo"],r[scope+"_hi"])):
        tie=list(tie)
        lo,hi=len(rebuilt)+1,len(rebuilt)+len(tie)
        mid=F(lo+hi,2)
        for i,r in enumerate(tie,lo):
            x=1-2*(mid-1)/(n-1) if n>1 else None
            rebuilt.append(dict(r,rank=i,rank_lo=lo,rank_hi=hi,midrank=float(mid),normalized_position=float(x) if x is not None else "",
                                normalized_position_exact=str(x) if x is not None else "",
                                rank_tie_size=len(tie),source_rank_scope=scope))
    return rebuilt

def pair_auc(records):
    p=[r["midrank"] for r in records if r["screen_positive"]]
    m=[r["midrank"] for r in records if not r["screen_positive"]]
    if not p or not m:
        return None,0
    total=F(sum(2*(a<b)+(a==b) for a in p for b in m),2)
    u=F(0)
    cursor=1
    for value in sorted(set(p+m)):
        size=p.count(value)+m.count(value)
        u+=p.count(value)*F(2*cursor+size-1,2)
        cursor+=size
    rank_auc=1-(u-F(len(p)*(len(p)+1),2))/(len(p)*len(m))
    assert rank_auc==total/(len(p)*len(m))
    return rank_auc,len(p)*len(m)

def w_recall(records,low,high):
    """Exact membership probability, expected count, and attainable bounds."""
    results={g:{"frozen":0,"expected":F(0),"minimum":0,"maximum":0} for g in (0,1)}
    weights={}
    crossing=[]
    for (lo,hi), group in groupby(records,key=lambda r:(r["rank_lo"],r["rank_hi"])):
        block=list(group)
        n=len(block)
        assert n==hi-lo+1
        k=max(0,min(hi,high)-max(lo,low)+1)
        for r in block:
            cls=r["screen_positive"]
            results[cls]["frozen"]+=int(low<=r["rank"]<=high)
            weights[r["slot"]]=F(k,n)
        c1=sum(r["screen_positive"] for r in block)
        for cls,c in ((0,n-c1),(1,c1)):
            q=results[cls]
            q["expected"]+=F(k*c,n)
            q["minimum"]+=max(0,k-(n-c))
            q["maximum"]+=min(k,c)
        if 0<k<n:
            crossing.append(dict(tie_lo=lo,tie_hi=hi,tie_size=n,budget_take=k,tie_positive=c1,tie_negative=n-c1,
                                 tie_member_slots=" ".join(x["slot"] for x in block),
                                 frozen_selected_slots=" ".join(x["slot"] for x in block if low<=x["rank"]<=high)))
    assert sum(weights.values())==high-low+1
    assert sum(x["expected"] for x in results.values())==high-low+1
    for x in results.values():
        assert x["minimum"]<=x["frozen"]<=x["maximum"]
        assert x["minimum"]<=x["expected"]<=x["maximum"]
    return results,weights,crossing

def metric(meta,rr, detail=True):
    n=len(rr)
    class_n=Counter(r["screen_positive"] for r in rr)
    out=dict(**meta,N=n,K=class_n[1],negative_N=class_n[0],score_levels=len({r["rank_lo"] for r in rr}),
             execution_status="evaluable" if n else "NA")
    for cls,name in ((0,"negative"),(1,"positive")):
        value=[r["midrank"] for r in rr if r["screen_positive"]==cls]
        out[f"{name}_median_midrank"]=median(value) if value else ""
    if not n:
        return out,[],[],[],[],[]
    auc,_=pair_auc(rr)
    out.update(AUROC=float(auc) if auc is not None else "", AUROC_exact=str(auc) if auc is not None else "")
    by_role=defaultdict(list)
    for r in rr:
        by_role[r["bhlh_id"]].append(r)
    bhlh=[]
    concord=F(0)
    comp=0
    for ident,subset in sorted(by_role.items()):
        ca,pairs=pair_auc(subset)
        if ca is not None:
            concord+=ca*pairs
            comp+=pairs
        bhlh.append(dict(**meta,bhlh_id=ident,bhlh_symbol=subset[0]["bhlh_symbol"],N=len(subset),
                         K=sum(r["screen_positive"] for r in subset),
                         positive_negative_comparisons=pairs,concordance=float(ca) if ca is not None else "",
                         concordance_exact=str(ca) if ca is not None else ""))
    out.update(bhlh_informative_groups=sum(r["positive_negative_comparisons"]>0 for r in bhlh),
               bhlh_positive_negative_comparisons=comp,bhlh_weighted_concordance=float(concord/comp) if comp else "",
               bhlh_weighted_concordance_exact=str(concord/comp) if comp else "")
    tails=[]; tie_rows=[]; quota_rows=[]; bins=[]
    for budget in (10,20):
        b=min(budget,n)
        for win,low,high in (("top",1,b),("bottom",n-b+1,n)):
            counts,weights,cross=w_recall(rr,low,high)
            if detail:
                for r in rr:
                    v=weights[r["slot"]]
                    r[f"{win}{budget}_frozen"]=int(low<=r["rank"]<=high)
                    r[f"{win}{budget}_tie_weight"]=float(v)
                    r[f"{win}{budget}_tie_weight_exact"]=str(v)
            for t in cross:
                tie_rows.append(dict(**meta,N=n,budget=budget,window=win,window_lo=low,window_hi=high,**t))
            quotas=Counter(r["bhlh_id"] for r in rr if low<=r["rank"]<=high)
            fixed_pos=F(0)
            for ident,quota in sorted(quotas.items()):
                qset=by_role[ident]
                gk=sum(r["screen_positive"] for r in qset)
                exp=F(quota*gk,len(qset))
                fixed_pos+=exp
                quota_rows.append(dict(**meta,budget=budget,window=win,bhlh_id=ident,bhlh_symbol=qset[0]["bhlh_symbol"],
                                      quota=quota,group_N=len(qset),group_K=gk,fixed_positive_expected=float(exp),
                                      fixed_positive_expected_exact=str(exp)))
            for cls,t in counts.items():
                k=class_n[cls]
                exp=t["expected"]
                fx=fixed_pos if cls else b-fixed_pos
                stats=dict(**meta,N=n,K=class_n[1],budget=budget,effective_budget=b,window=win,
                           win_lo=low,win_hi=high,windows_overlap=int(2*b>n),screen_positive=cls,class_denominator=k,
                           frozen_count=t["frozen"],tie_expected_count=float(exp),tie_expected_count_exact=str(exp),
                           tie_min_count=t["minimum"],tie_max_count=t["maximum"],crossing_groups=len(cross),
                           tie_count_variable=int(t["minimum"]!=t["maximum"]),tie_recall_pct=float(100*exp/k) if k else "",
                           frozen_recall_pct=float(F(100*t["frozen"],k)) if k else "",
                           precision_positive_frozen=float(F(t["frozen"],b)) if cls==1 and win=="top" else "",
                           precision_positive_tie_expected=float(exp/b) if cls==1 and win=="top" else "",
                           uniform_random_recall_pct=float(F(100*b,n)) if k else "",
                           uniform_expected_count=float(F(b*k,n)),fixed_bhlh_expected_count=float(fx),
                           fixed_bhlh_expected_count_exact=str(fx))
                tails.append(stats)
                prefix=f"{win}{budget}_{'positive' if cls else 'negative'}"
                for field in ("frozen_count","tie_expected_count","tie_min_count","tie_max_count","fixed_bhlh_expected_count"):
                    out[f"{prefix}_{field}"]=stats[field]
                if budget==10 and ((win=="top" and cls) or (win=="bottom" and not cls)):
                    out[f"{prefix}_recall_pct"]=stats["tie_recall_pct"]
                    out[f"{prefix}_variable"]=stats["tie_count_variable"]
    out["random_B10_recall_pct"]=float(F(100*min(10,n),n))
    if detail and n>1:
        for cls in (0,1):
            hist=Counter()
            for r in rr:
                if r["screen_positive"]==cls:
                    x=F(r["normalized_position_exact"]); ix=min(19,int((x+1)*10))
                    assert 0<=ix<=19
                    hist[ix]+=1
            assert sum(hist.values())==class_n[cls]
            for j in range(20):
                k=class_n[cls]
                bins.append(dict(**meta,screen_positive=cls,bin=j,lo=float(F(j,10)-1),hi=float(F(j+1,10)-1),
                                 count=hist[j],class_size=k,class_mass=float(F(hist[j],k)) if k else ""))
    return out,tails,bhlh,quota_rows,tie_rows,bins

def slot_sha(slots):
    return sha256(("\n".join(sorted(slots))+"\n").encode()).hexdigest()

def scope_reason(cohort,ep,m,scope,n):
    if not n:
        return "no eligible pair ranking in this scope"
    if cohort=="old" and ep=="GSE53951" and m in A and scope=="reference":
        return "main panel A"
    if cohort=="followup" and ep=="GSE43916" and m in B and scope=="native":
        return "panel B uses shared physical slots"
    if cohort=="old" and ep=="GSE53951" and m in C and scope=="native":
        return "panel C uses shared physical slots"
    if cohort=="followup" and m in A:
        return "archived-method carry-over or cross-input comparison"
    if m.startswith(("collectri_","imcell_weight_","imcell_p0")) or m.endswith(("_joint","_additive")):
        return "prior, propagation, or precision sensitivity"
    if m in ("integration_expression","integration_response","network_add","two_add","expr_product") or "diag" in m:
        return "component, alternate scoring rule, or support diagnostic"
    return "auxiliary method, reference subset, or endpoint outside the three fixed panels"

def main(opt=None):
    if opt is None:
        ap=ArgumentParser(description=__doc__)
        for arg in ("old","new","out"):
            ap.add_argument("--"+arg,type=Path,required=True)
        opt=ap.parse_args()
    out=opt.out
    out.mkdir(parents=True,exist_ok=True)
    old=load(opt.old,"old")
    new=load(opt.new,"followup")
    check_slots(old,new)
    mappings={"A":(old,"old","GSE53951",A),"B":(new,"followup","GSE43916",B),"C":(old,"old","GSE53951",C)}
    all_slots=[]; panel_slots=[]; summary=[]; tails=[]; strata=[]; quotas=[]; ties=[]; bins=[]; checks=[]
    full={}
    for pane,(groups,cohort,ep,methods) in mappings.items():
        selected=[groups[ep,m] for m in methods]
        slots=(valid(selected[0],"reference") if pane=="A" else valid(selected[0],"native")&valid(selected[1],"native"))
        n=len(slots)
        assert n>1
        keyhash=slot_sha(slots)
        for p in sorted(slots):
            r=selected[0][p]
            panel_slots.append(dict(panel=pane,endpoint=ep,slot=p,pair_id=r["pair_id"],bhlh=r["bhlh_symbol"],
                                    partner=r["partner_symbol"],screen_positive=r["screen_positive"],source_locator=r["source_locator"]))
        for m,group in zip(methods,selected):
            assert slots<=valid(group,"native"),(pane,m)
            rank=rerank(group,slots)
            if pane=="A":
                assert slots==valid(group,"reference"),(pane,m)
                for r in rank:
                    assert (r["rank"],r["rank_lo"],r["rank_hi"])==(r["reference_rank"],r["reference_lo"],r["reference_hi"]),(pane,m,r["slot"])
            meta=dict(panel=pane,endpoint=ep,method=m,display_name=LABELS[m],source_round=cohort)
            met,t,st,q,ti,bn=metric(meta,rank)
            summary.append(met); tails+=t; strata+=st; quotas+=q; ties+=ti; bins+=bn
            all_slots.extend(dict(**meta,**r) for r in rank)
            full[pane,m]={r["slot"]:r for r in rank}
            checks.append(dict(**meta,N=n,K=met["K"],negative_N=n-met["K"],set_sha256=keyhash,
                               method_native_N=len(valid(group,"native")),reference_intersection_N=len(valid(group,"reference")),
                               check="native ranks reprojected; common ranks match" if pane=="A" else "native eligible slot intersection reranked"))
    old_c=old["GSE53951","tfcomb"]
    c65=set(full["C","tfcomb"]); c60=valid(old_c,"reference")
    assert c60<c65,(c60-c65,c65-c60)
    crows=[]
    for p,r in sorted(full["C","tfcomb"].items(),key=lambda q:q[1]["rank"]):
        crows.append(dict(slot=p,pair_id=r["pair_id"],bhlh=r["bhlh_symbol"],partner=r["partner_symbol"],
                          screen_positive=r["screen_positive"],in_old_common=int(p in c60),tfcomb_rank_on_65=r["rank"],
                          tfcomb_rank_on_old_common_60=r["reference_rank"],source_locator=r["source_locator"]))
    appendix=[]; appendix_tails=[]
    for source,groups in (("old",old),("followup",new)):
        for (ep,m),group in sorted(groups.items()):
            for scope in ("native","reference"):
                valid_slots=valid(group,scope)
                # GSE43916 follow-up reference is the older GSE53951 panel; explicitly keep its intersection size.
                if scope=="reference" and {r["reference_full_N"] for r in group.values()}=={""}:
                    continue
                meta=dict(source_round=source,endpoint=ep,method=m,scope=scope,
                          reference_full_N=next(iter(group.values()))["reference_full_N"] if scope=="reference" else "",
                          scope_in_figure=scope_reason(source,ep,m,scope,len(valid_slots)))
                met,t,*_=metric(meta,rerank(group,valid_slots,scope),detail=False)
                if not valid_slots:
                    reasons=Counter(r["exclusion"] for r in group.values() if r["exclusion"])
                    met["most_common_nonrank_reasons"]=" | ".join(f"{x} ({y})" for x,y in reasons.most_common(3))
                appendix.append(met)
                appendix_tails.extend(t)
    for name,data in (("PLOT_SLOTS",all_slots),("PANEL_SLOTS",panel_slots),("PANEL_SCOPE_CHECKS",checks),
                      ("PANEL_METHOD_SUMMARY",summary),("HEAD_TAIL_STATS",tails),("BOUNDARY_TIES",ties),
                      ("PLOT_BINS",bins),("BHLH_CONCORDANCE",strata),("BHLH_QUOTAS",quotas),
                      ("TFComb_NATIVE65_VS_COMMON60",crows),("METHOD_SCOPE_APPENDIX",appendix),
                      ("METHOD_SCOPE_HEAD_TAIL",appendix_tails)):
        write(out/(name+".tsv"),data)
    inputs=[dict(role=k,file=p.name,bytes=p.stat().st_size,sha256=digest(p)) for k,p in (("old",opt.old),("followup",opt.new))]
    (out/"CORE_REBUILD.json").write_text(json.dumps(dict(inputs=inputs,slots=len(all_slots),main_budgets=len(tails),
                                                         score_levels_definition="unique original tie groups in each projected scope",
                                                         tie_validation="rational membership and attainable bounds; AUROC pair comparisons and exact rank-sum"),
                                                   indent=2)+"\n",encoding="utf-8")
    print([(r["panel"],r["method"],r["N"],r["K"]) for r in summary])
    return old, new, full

if __name__=="__main__":
    main()
