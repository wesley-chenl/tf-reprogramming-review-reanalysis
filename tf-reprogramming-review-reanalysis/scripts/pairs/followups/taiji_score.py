"""Apply the frozen Taiji-reprogram K=2 adaptation to native Taiji PageRank."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re

from taiji_rna import ROOT, DESIGN, OUT, SCORE, SCORE_COLUMNS, BENCH, emit, load, sha

def finite_pr(x):
    try:
        p = float(x)
    except (TypeError, ValueError):
        return None
    return p if math.isfinite(p) and p > 0 else None

def source_table(path, motif_path):
    native = load(path)
    key, source, target = native[0]
    if {source,target} != {"source","target"}:
        raise ValueError(f"Expected two states named source, target: {list(native[0])}")
    motifs = {line.split()[1].split("+",1)[0].casefold() for line in motif_path.read_text().splitlines()
              if line.startswith("MOTIF ")}
    dup = Counter(x[key].casefold() for x in native)
    if max(dup.values()) != 1:
        raise ValueError("Native Taiji PageRank symbol keys are not unique case-insensitively")
    q = {}
    for raw in native:
        symbol = raw[key]
        src, dst = finite_pr(raw["source"]), finite_pr(raw["target"])
        ratio = dst/src if src is not None and dst is not None else None
        m = max(ratio,1/ratio) if ratio and math.isfinite(ratio) else None
        q[symbol.casefold()] = {
            "native_symbol": symbol, "PR_source_raw": raw["source"], "PR_target_raw": raw["target"],
            "PR_source": src, "PR_target": dst, "PR_ratio_target_source": ratio, "PR_ratio_abs": m,
            "log_ratio_abs": math.log(m) if m is not None else None,
            "native_cisBP_motif_name": int(symbol.casefold() in motifs),
        }
    valid = [c for c in q.values() if c["PR_ratio_abs"] is not None]
    assert valid and any(r["PR_source"] != r["PR_target"] for r in valid), "Invalid or identical state ranks"
    emit(ROOT / "taiji" / "TAIJI_TF_PAGERANK.tsv", list(q.values()))
    return q

def annotate_members(q):
    expr = {r["gene_symbol"]: r for r in load(OUT / "FULL_GENE_TPM.tsv")}
    configured = load(OUT / "TAIJI_RNA_CANDIDATE_MAP.tsv")
    annotation = load(BENCH / "inputs" / "annotation" / "GENES_M24.tsv")
    symbols = Counter(p["gene_symbol"].casefold() for p in annotation)
    out = {}
    for c in configured:
        rname = c["rna_gene_symbol"]
        en = expr.get(rname, {})
        kd = c["candidate_symbol"].casefold()
        # For an old RefSeq symbol, only fall back to that already Entrez-verified
        # alias, provided it identifies an unambiguous native annotation symbol.
        alt = rname.casefold() if rname else ""
        pr = q.get(kd) if symbols[kd] == 1 else None
        if not pr and alt and symbols[alt] == 1:
            pr = q.get(alt)
        srna, trna = (float(en[f]) if f in en else None for f in ("source_mean_tpm","target_tpm"))
        oe = int(trna>srna) if srna is not None and trna is not None else None
        if c["candidate_id"] in {"","NA","unresolved"}:
            cause = "identity_unresolved"
        elif oe is None:
            cause = "RNA_missing_or_ambiguous"
        elif not pr or pr["PR_ratio_abs"] is None:
            cause = "native_PR_missing_or_nonpositive"
        elif not oe:
            cause = "fails_RNA_double_OE_gate"
        else:
            cause = ""
        valid = bool(not cause)
        out[(c["candidate_symbol"],c["candidate_id"])] = {
            **c, "native_symbol": pr["native_symbol"] if pr else "",
            "PR_source": pr["PR_source"] if pr else None,
            "PR_target": pr["PR_target"] if pr else None,
            "PR_ratio_target_source": pr["PR_ratio_target_source"] if pr else None,
            "PR_ratio_abs": pr["PR_ratio_abs"] if pr else None,
            "log_ratio_abs": pr["log_ratio_abs"] if pr else None,
            "RNA_source_mean_TPM": srna, "RNA_target_TPM": trna, "RNA_OE": oe, 
            "eligible": int(valid), "exclusion_reason": cause,
        }
    emit(ROOT / "taiji" / "TAIJI_CANDIDATE_MEMBERS.tsv", list(out.values()))
    return out

def pairs(members):
    output = []
    for p in load(DESIGN):
        a,b = (members[(p[f"{side}_symbol"], p[f"{side}_id"])] for side in ("bhlh","partner"))
        eligible = bool(a["eligible"] and b["eligible"] and a["candidate_id"] != b["candidate_id"])
        product = a["PR_ratio_abs"]*b["PR_ratio_abs"] if eligible else None
        logsum = a["log_ratio_abs"] + b["log_ratio_abs"] if eligible else None
        if eligible:
            if not math.isfinite(product) or not math.isclose(math.log(product),logsum,rel_tol=2e-12,abs_tol=2e-15):
                raise ArithmeticError("Non-finite or non-equivalent pair product")
        reasons = ";".join(f"{side}:{r['exclusion_reason']}" for side,r in (("bhlh",a),("partner",b))
                           if r["exclusion_reason"])
        if a["candidate_id"] == b["candidate_id"]:
            reasons = (reasons + ";" if reasons else "") + "identical_resolved_member_IDs"
        output.append({
            "endpoint": "GSE43916", "variant": p["variant"], "pair_id": p["pair_id"],
            "method": "taiji_double_oe", "score": "" if product is None else format(product,".17g"),
            "score_num": "", "score_den": "", "eligible": int(eligible), "exclusion_reason": reasons,
            "details": json.dumps({"bhlh":a,"partner":b,"log_score":logsum},separators=(",",":")),
        })
    assert len({(p["variant"],p["pair_id"]) for p in output}) == 598*9 == len(output)
    return output

def main():
    parse = argparse.ArgumentParser(description=__doc__)
    parse.add_argument("gene_rank",type=Path)
    parse.add_argument("motif",type=Path)
    a = parse.parse_args()
    assert (ROOT / "ADDENDUM_PROTOCOL.yaml").is_file()
    q = source_table(a.gene_rank,a.motif)
    members = annotate_members(q)
    result = pairs(members)
    emit(SCORE / "TAIJI_SCORES.tsv",result,SCORE_COLUMNS)
    validation = {"native_pagerank_sha256":sha(a.gene_rank), "motif_sha256":sha(a.motif),
                  "native_rows":len(q), "candidate_rows":len(members), "returned_rows":len(result),
                  "main_eligible":sum(r["eligible"] for r in result if r["variant"]=="main"),
                  "pair_score_test":"native ratio multiplication versus log-ratio addition"}
    (ROOT / "taiji/TAIJI_RANK_VALIDATION.json").write_text(json.dumps(validation,indent=2)+"\n")

if __name__ == "__main__":
    main()
