"""Score the registered deterministic extensions inside the offline allowlist."""
from collections import Counter, defaultdict
import csv
from fractions import Fraction
import json
from math import fsum, sqrt
from pathlib import Path
import sys

READS=set()
NAMES={"CLEAN_PAIRS.tsv","SCORING_POLICY.json","GSE53951_EFFECT.tsv","GSE114809_EFFECT.tsv","TRRUST_MOUSE_SIGNED.tsv"}
GENES=("bhlh_id","partner_id")

def restrict(event,args):
    if event.startswith(("socket.","subprocess.")):
        raise PermissionError(event)
    if event=="open" and isinstance(args[0],str) and isinstance(args[1],str):
        p=Path(args[0]).absolute()
        if args[1]=="r" and ((p.parent==Path("/input") and p.name in NAMES)
                            or p in {Path("/out/all_pair_scores.tsv"),Path("/out/program.tsv")}):
            READS.add(str(p))
        elif args[1]=="r" and str(p).startswith("/usr/"):
            pass
        elif p.parent==Path("/out") and "w" in args[1]:
            pass
        else:
            raise PermissionError(p)

def read(p):
    with Path(p).open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def sid(s):
    return s.removeprefix("GeneID:")

def fformat(x):
    return "" if x is None else format(float(x),".17g")

def result(r,method,value,reason="",details=None):
    is_exact=isinstance(value,Fraction)
    return {"endpoint":r["endpoint"],"variant":r["variant"],"method":method,"pair_id":r["pair_id"],
            "score":fformat(value),"score_num":value.numerator if is_exact else "","score_den":value.denominator if is_exact else "",
            "reason":reason,"details":json.dumps(details or {},separators=(",",":"),sort_keys=True)}


def path_members(edges, panel, measured, programme):
    adj=defaultdict(list)
    for a,b,s in edges:
        adj[a].append((b,s))
    support={}
    for src in sorted(panel):
        directions=defaultdict(set)
        for middle,s in adj[src]:
            if middle==src:
                continue
            directions[middle].add(s)
            for target,t in adj[middle]:
                if target!=src and target!=middle:
                    directions[target].add(s*t)
        valid={g:next(iter(s)) for g,s in directions.items() if g in measured and g not in panel and len(s)==1}
        support[src]={"degree":len(valid),"hits":{g for g,sgn in valid.items() if programme.get(g) == sgn},
                      "conflicting":sum(len(v)==2 and k in measured and k not in panel for k,v in directions.items())}
    return support


class SignalFlow:
    def __init__(self,edges):
        self.nodes=sorted({a for a,b,s in edges}|{b for a,b,s in edges})
        self.index={g:i for i,g in enumerate(self.nodes)}
        i_degree=Counter(b for a,b,s in edges)
        o_degree=Counter(a for a,b,s in edges)
        self.links=[(self.index[a],self.index[b],0.5*s/sqrt(o_degree[a]*i_degree[b])) for a,b,s in edges]

    def multiply(self,v):
        nxt=[0.0]*len(v)
        for a,b,w in self.links:
            nxt[b]+=v[a]*w
        return nxt

    def solve(self,b):
        base=[0.5*x for x in b]
        answer=base.copy()
        for _ in range(48):
            answer=[x+y for x,y in zip(base,self.multiply(answer))]
        residual=max(abs(x-y-z) for x,y,z in zip(answer,base,self.multiply(answer)))
        return answer,residual

    def prepare(self, panel, measured, programme):
        sham_b=[measured.get(g, (0, 0))[1] for g in self.nodes]
        sham,residual=self.solve(sham_b)
        keep=[(self.index[g],measured[g][0]) for g in programme if g in self.index and g not in panel]
        effects={}
        for g in sorted(panel & self.index.keys()):
            data=[0.0]*len(sham)
            data[self.index[g]]=1
            delta,error=self.solve(data)
            residual=max(residual,error)
            effects[g]=[delta[i] for i,_ in keep]
        if residual>1e-11:
            raise ArithmeticError(f"SFA residual {residual}")
        return effects, [d for _,d in keep], residual

    @staticmethod
    def pair(effects,target,a,b):
        if a not in effects or b not in effects:
            return None,"member_missing_from_graph",{}
        if not target or not any(target):
            return None,"no_comparable_programme_nodes",{}
        z=[x+y for x,y in zip(effects[a],effects[b])]
        length=sqrt(fsum(x*x for x in z))
        if length==0:
            return 0.0,"valid_zero_response",{"target_n":len(target),"l2":0.0}
        value=fsum(x*y for x,y in zip(z,target))/(length*sqrt(fsum(x*x for x in target)))
        return value,"",{"target_n":len(target),"l2":length}


def main():
    sys.addaudithook(restrict)
    design=read("CLEAN_PAIRS.tsv")
    panel={sid(r[c]) for r in design for c in GENES if r[c]}
    edges=[(r["tf_gene_id"],r["target_gene_id"],int(r["sign"])) for r in read("TRRUST_MOUSE_SIGNED.tsv")]
    old=read("/out/all_pair_scores.tsv")
    old_programme=read("/out/program.tsv")
    scores,diagnostics=[],[]
    sfa=SignalFlow(edges)
    for endpoint in ("GSE53951","GSE114809"):
        members=read(endpoint+"_EFFECT.tsv")
        measured={sid(r["stable_id"]):(float(r["d_g"]),float(r["mean_source"])) for r in members}
        program={sid(r["stable_id"]):int(r["direction"]) for r in old_programme if r["endpoint"]==endpoint}
        by_gene=path_members(edges,panel,measured,program)
        effect,target,residual=sfa.prepare(panel,measured,program)
        diagnostics.append({"endpoint":endpoint,"sfa_residual":residual,"SFA_programme_n":len(target),
                            "program_n":len(program),"nodes":len(sfa.nodes),"edges":len(edges),
                            "node_with_conflicting_measured_paths_n":sum(v["conflicting"]>0 for v in by_gene.values())})
        for row in (r for r in old if r["endpoint"]==endpoint):
            flat=Fraction(int(row["num"]),int(row["den"])) if row["scoreable"]=="1" else None
            fields={"bhlh_p":row["bhlh_p"],"partner_p":row["partner_p"],"bhlh_hits_n":int(row["bhlh_hits_n"]),
                    "partner_hits_n":int(row["partner_hits_n"]),"intersection_n":int(row["intersection_n"]),
                    "bhlh_measured_degree":int(row["bhlh_measured_degree"]),
                    "partner_measured_degree":int(row["partner_measured_degree"]),"program_n":int(row["program_n"])}
            scores.append(result(row,row["method"],flat,row["reason"],fields))
            if row["method"]!="expression":
                continue
            product=Fraction(row["bhlh_p"])*Fraction(row["partner_p"]) if flat is not None else None
            scores.append(result(row,"expr_product",product,row["reason"],fields))
            a,b=(sid(row[c]) for c in GENES)
            both=a in by_gene and b in by_gene and a!=b
            if both:
                x,y=by_gene[a],by_gene[b]
                p,q=x["hits"],y["hits"]
                info={"member_degree":[x["degree"],y["degree"]],"member_hits":[sorted(p),sorted(q)],
                      "intersection_ids":sorted(p&q),"conflicting_measured":[x["conflicting"],y["conflicting"]],"program_n":len(program)}
                reason="" if x["degree"] and y["degree"] else "member_has_no_unambiguous_measured_reachable_target"
                for method,n in (("two_add",len(p)+len(q)),("two_union",len(p|q))):
                    scores.append(result(row,method,Fraction(n,len(program)) if not reason else None,reason,info))
                score,note,extra=sfa.pair(effect,target,a,b)
                rev,_,_=sfa.pair(effect,target,b,a)
                if score!=rev:
                    raise ValueError("SFA joint response order dependence")
                scores.append(result(row,"sfa",score,note,extra))
            else:
                for method in ("two_add","two_union","sfa"):
                    scores.append(result(row,method,None,"unresolved_or_equal_gene_id"))
    with Path("/out/core.tsv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(scores[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(scores)
    Path("/out/CORE_AUDIT.json").write_text(json.dumps({"opened":sorted(READS),"SFA":diagnostics},indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":
    main()
