"""TFcomb exact native combination helper with external graph/bulk feature adapter."""
import argparse
import ast
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy
from scipy.sparse import csr_matrix
import sklearn
from sklearn.linear_model import Ridge


def source_fn(path, name, context):
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    methods = [v for v in tree.body if isinstance(v, ast.FunctionDef) and v.name == name]
    if len(methods) != 1:
        raise ValueError(f"Expected exactly one official {name} function")
    fn = methods[0]
    if fn.decorator_list:
        raise ValueError("Unexpected decorator on selected official function")
    exec(compile(ast.Module(body=[fn], type_ignores=[]), str(path), "exec"), context)
    return context[name]


def records(path):
    with path.open(encoding="utf-8-sig", newline="") as s:
        return list(csv.DictReader(s, delimiter="\t"))


def write(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as s:
        w = csv.DictWriter(s, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def only_whitelist(event, args):
    if event.startswith(("socket.", "subprocess.")):
        raise PermissionError("Network and subprocess are not needed by this score")
    if event == "open":
        p, mode, _ = args
        if not isinstance(p, str) or not isinstance(mode, str):
            return
        q = Path(p).absolute()
        if str(q).startswith(("/inputs/", "/authors/", "/out/", str(Path(sys.prefix).absolute()) + "/", "/usr/")):
            return
        raise PermissionError(f"Unexpected file access: {p}")


def score(inp, out, source, endpoint):
    sys.addaudithook(only_whitelist)
    pairs = records(inp / "CLEAN_PAIRS.tsv")
    effect = records(inp / f"{endpoint}_EFFECT.tsv")
    edge = records(inp / "TRRUST_MOUSE_SIGNED.tsv")
    reads = ["CLEAN_PAIRS.tsv", f"{endpoint}_EFFECT.tsv", "TRRUST_MOUSE_SIGNED.tsv"]
    assert len(edge) == 3964
    assert len({r["stable_id"] for r in effect}) == len(effect)
    effect = {r["stable_id"]: r for r in effect if r["d_g"] not in ("", "NA", "NaN", "Inf", "-Inf")}
    genes = sorted(effect)
    idx = {gene: i for i, gene in enumerate(genes)}
    source_avg = np.asarray([float(effect[g]["mean_source"]) for g in genes])
    target_avg = np.asarray([float(effect[g]["mean_target"]) for g in genes])
    y = target_avg - source_avg
    frozen_d = np.asarray([float(effect[g]["d_g"]) for g in genes])
    np.testing.assert_allclose(y, frozen_d, rtol=2e-14, atol=2e-14)
    coll = [(idx[e["target_gene_id"]], idx[e["tf_gene_id"]], int(e["sign"])) for e in edge
            if e["target_gene_id"] in idx and e["tf_gene_id"] in idx]
    row, col, val = zip(*coll)
    assert len(set(zip(row, col))) == len(row) and set(val) <= {1, -1}
    n = len(genes)
    adjacency = csr_matrix((val, (row, col)), shape=(n, n), dtype=float)
    regulator_indices = sorted(set(col))
    regulators = [genes[j] for j in regulator_indices]
    rix = {g: j for j, g in enumerate(regulators)}
    direct = adjacency[:, regulator_indices].tocsr()
    hop2 = adjacency @ direct
    hop3 = adjacency @ hop2
    design = 0.8 * direct + 0.15 * hop2 + 0.05 * hop3
    fitted = Ridge(alpha=1).fit(design.toarray(), y)
    del design, adjacency, hop2, hop3
    candidate_ids = sorted({r[c].removeprefix("GeneID:") for r in pairs for c in ("bhlh_id", "partner_id") if r[c]})
    candidate_with_columns = [g for g in candidate_ids if g in rix]
    positions = [rix[g] for g in candidate_with_columns]
    native_cols = pd.DataFrame(direct[:, positions].toarray().T, index=candidate_with_columns, columns=genes)
    native_fit = SimpleNamespace(coef_=fitted.coef_[positions], intercept_=float(fitted.intercept_))
    # Upstream tqdm only reports progress. The module AST keeps both functions byte-identical.
    ns = {"np": np, "pd": pd, "tqdm": lambda iterable: iterable}
    direct_native = source_fn(source / "tf_inference.py", "get_directing_score", ns)
    helper = source_fn(source / "utils.py", "get_multi_TF", ns)
    native = helper(candidate_with_columns, native_fit, native_cols, {}, source_avg, target_avg, "pos", 2)
    scored = {tuple(sorted(pair)): value for pair, value in zip(native["TF_combination"], native["directing_score"])}
    eligible_unordered = set()
    rows = []
    audits = {}
    for p in pairs:
        a, b = (p[c].removeprefix("GeneID:") for c in ("bhlh_id", "partner_id"))
        ga, gb = tuple(sorted((a, b)))
        aa = float(fitted.coef_[rix[a]]) if a in rix else None
        bb = float(fitted.coef_[rix[b]]) if b in rix else None
        value = None
        if not a or not b:
            reason = "unresolved_identity"
        elif aa is None or bb is None:
            reason = "member_no_measured_outgoing_in_model"
        elif a == b:
            reason = "identical_members_native_distinct_filter"
        elif aa <= 0 or bb <= 0:
            reason = "native_direction_positive_filter"
        elif (ga, gb) not in scored:
            reason = "native_helper_missing_pair"
        elif not np.isfinite(scored[ga, gb]):
            reason = "native_PCC_nonfinite_constant_or_invalid_readout"
        else:
            value = float(scored[ga, gb])
            reason = ""
            eligible_unordered.add((ga, gb))
        details = {
            "bhlh_coef": aa, "partner_coef": bb, "measured_gene_n": n,
            "fitted_tf_n": len(regulators), "raw_direct_nonzero_n": len(coll),
            "intercept": float(fitted.intercept_), "native_X_none": True,
        }
        item = {"method": "tfcomb", "endpoint": endpoint, "variant": p["variant"], "pair_id": p["pair_id"],
                "score": "" if value is None else format(value, '.17g'), "reason": reason,
                "details": json.dumps(details, sort_keys=True, separators=(",",":"))}
        rows.append(item)
        audit_key = hashlib.sha256(("native_module_audit_v1|" + p["pair_id"] + "|" + a + "|" + b).encode()).hexdigest()
        audits[audit_key] = item, a, b
    checks = []
    chosen = []
    for good, count in ((True,12),(False,8)):
        chosen.extend([(key, record) for key, record in sorted(audits.items()) if bool(record[0]["score"]) == good][:count])
    for key, (item, a, b) in chosen:
        value = float(item["score"]) if item["score"] else None
        reverse = reconstructed = None
        if value is not None:
            reverse = float(direct_native([b,a], native_fit, native_cols, {}, y, mode="multi", X=None))
            pred = direct[:, rix[a]].toarray().ravel() * fitted.coef_[rix[a]] + direct[:, rix[b]].toarray().ravel() * fitted.coef_[rix[b]]
            pred -= pred.mean()
            yc = y - y.mean()
            reconstructed = float(pred.dot(yc) / (np.linalg.norm(pred)*np.linalg.norm(yc)))
            np.testing.assert_allclose((reverse,reconstructed),value,rtol=1e-12,atol=1e-13)
        elif a and b and a in rix and b in rix and a != b:
            assert fitted.coef_[rix[a]] <= 0 or fitted.coef_[rix[b]] <= 0 or not np.isfinite(scored.get(tuple(sorted((a,b))),np.nan))
        checks.append({"audit_hash":key,"a":a,"b":b,"score":value,"reverse_native":reverse,
                       "independent_pearson_sparse":reconstructed,"reason":item["reason"]})
    write(out/f"tfcomb_{endpoint}.tsv",rows,["method","endpoint","variant","pair_id","score","reason","details"])
    predicted = []
    for a,b in sorted(eligible_unordered):
        pred = direct[:, rix[a]].toarray().ravel() * fitted.coef_[rix[a]] + direct[:, rix[b]].toarray().ravel() * fitted.coef_[rix[b]]
        for k in np.flatnonzero(pred):
            predicted.append({"a":a,"b":b,"gene_id":genes[k],"predicted_delta_without_intercept":format(pred[k],".17g"),
                              "observed_delta":format(y[k],".17g"),"intercept":format(fitted.intercept_,".17g")})
    write(out/f"tfcomb_predicted_{endpoint}.tsv",predicted,list(predicted[0]) if predicted else
          ["a","b","gene_id","predicted_delta_without_intercept","observed_delta","intercept"])
    (out/f"tfcomb_audit_{endpoint}.json").write_text(json.dumps({
        "score_rows":len(rows),"positive_direction_native_pairs":len(scored),
        "variant_score_rows":sum(bool(r["score"]) for r in rows),"scored_study_gene_pairs":len(eligible_unordered),
        "reads":reads, "source_files":[f.name for f in source.iterdir()],
        "source_sha256":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in source.glob("*.py")},
        "python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,
        "sklearn":sklearn.__version__,"hash_selected_checks":checks,
    },indent=2,ensure_ascii=False,allow_nan=False),encoding="utf-8")
    print(f"{endpoint}: {len(rows)} variant rows; {len(eligible_unordered)} distinct studied pairs with valid positive TFcomb PCC; predicted vector stored sparsely.")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("inp",type=Path)
    parser.add_argument("out",type=Path)
    parser.add_argument("source",type=Path)
    parser.add_argument("endpoint")
    ns=parser.parse_args()
    score(ns.inp,ns.out,ns.source,ns.endpoint)
