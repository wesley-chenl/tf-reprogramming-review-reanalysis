"""Freeze the exact archive prior and explicit IMCell target programs."""
from collections import Counter, defaultdict
import math
from pathlib import Path
import subprocess

from r2lib import BASE, BASE_ENV, ROOT, TASKS, input_record, read_json, read_tsv, run_env, sha, write_json, write_tsv

output = ROOT / "inputs_manifest"
release = read_json(ROOT / "config/SOURCE_RELEASES.json")
genes = read_tsv(BASE / "inputs/annotation/GENES_M24.tsv")
ids = {r["gene_id"] for r in genes}
exact = defaultdict(set)
for row in genes:
    exact[row["gene_symbol"]].add(row["gene_id"])
alias = read_json(BASE / "inputs/annotation/SYMBOL_TO_STABLE_ID.json")


def resolve(symbol):
    match = exact.get(symbol)
    if match:
        return (next(iter(match)), "exact_M24_symbol") if len(match) == 1 else (None, "ambiguous_exact_M24")
    match = alias.get(symbol) or []
    if isinstance(match, str):
        match = [match]
    match = set(match) & ids
    return (next(iter(match)), "frozen_unique_alias") if len(match) == 1 else (None, "unresolved_or_ambiguous_alias")


archive = ROOT / "sources/CompForce/data/mouse_chipX.rda"
archive_sha = sha(archive)
assert archive_sha == release["CompForce"]["file_hashes"]["data/mouse_chipX.rda"]
raw = ROOT / "checkpoints/mouse_chipX_native_extraction.tsv"
subprocess.run([str(BASE_ENV / "bin/Rscript"), str(ROOT / "scripts/export_chipx.R"), str(archive), str(raw)],
               check=True, env=run_env())
edges = read_tsv(raw)
symbols = sorted({r[c] for r in edges for c in ("TF", "TG")})
mappings = {s: resolve(s) for s in symbols}
write_tsv(output / "CHIPX_SYMBOL_MAPPING.tsv", [
    {"native_mouse_symbol": s, "stable_id": mappings[s][0], "resolution": mappings[s][1]} for s in symbols])
mapped = []
lost = Counter()
for row in edges:
    tf, tf_state = mappings[row["TF"]]
    tg, tg_state = mappings[row["TG"]]
    if not tf or not tg:
        lost["TF=" + tf_state + ";TG=" + tg_state] += 1
        continue
    mapped.append({"TF": tf, "TG": tg})
assert mapped
write_tsv(output / "COMPFORCE_MOUSE_CHIPX_M24.tsv", mapped, ["TF", "TG"])
distinct = len({(r["TF"], r["TG"]) for r in mapped})
signed_path = BASE / "inputs/priors/COLLECTRI_MOUSE_M24.tsv"
pairs = defaultdict(set)
disc = Counter()
for row in read_tsv(signed_path):
    tf, tg, weight = row["source"], row["target"], float(row["weight"])
    if tf not in ids or tg not in ids:
        disc["unresolved_ID"] += 1
    elif not math.isfinite(weight) or weight == 0:
        disc["nonfinite_or_zero_mor"] += 1
    else:
        sign = 1 if weight > 0 else -1
        disc["duplicate_identical_signed_edge"] += sign in pairs[(tf, tg)]
        pairs[(tf, tg)].add(sign)
signed = []
conflicts = []
for (tf, tg), signs in sorted(pairs.items()):
    if len(signs) > 1:
        conflicts.append({"TF": tf, "TG": tg, "reason": "conflicting activation and repression"})
        continue
    signed.append({"TF": tf, "TG": tg, "weight": 1, "type": signs.pop()})
assert {r["type"] for r in signed} == {-1, 1}
write_tsv(output / "IMCELL_COLLECTRI_MOUSE_SIGNED.tsv", signed, ["TF", "TG", "weight", "type"])
write_tsv(output / "IMCELL_SIGN_CONFLICTS.tsv", conflicts, ["TF", "TG", "reason"])
nodes = {r[k] for r in signed for k in ("TF", "TG")}
out_degree = Counter(r["TF"] for r in signed)
u = read_tsv(BASE / "config/CANDIDATE_UNIVERSE.tsv")
coverage = [{"stable_id": row["gene_id"], "symbol": row["gene_symbol"], "in_network": row["gene_id"] in nodes,
             "out_degree": out_degree[row["gene_id"]]} for row in u]
write_tsv(output / "IMCELL_SEED_COVERAGE.tsv", coverage)
programs = []
for task in TASKS:
    for row in read_tsv(BASE / f"inputs/rna/{task}/DESEQ2_TARGET_V_SOURCE.tsv"):
        lfc, padj = float(row["log2FoldChange"].replace("NA", "nan")), float(row["padj"].replace("NA", "nan"))
        if math.isfinite(lfc) and math.isfinite(padj) and padj < .05 and abs(lfc) > 1:
            stable_id = row["gene_id"]
            programs.append({"target": task, "stable_id": stable_id, "direction": "activate" if lfc > 1 else "repress",
                             "log2FoldChange": row["log2FoldChange"], "padj": row["padj"],
                             "mapped_to_graph": stable_id in nodes})
write_tsv(output / "IMCELL_TARGET_PROGRAM.tsv", programs)
decisions = [
    {"method": "CompForce RRA", "order": 1, "candidate_prior": "official mouse_chipX.rda / chipX_gs",
     "selected": True, "source": str(archive), "source_sha256": archive_sha,
     "reason": "Zenodo 20563029 archive and v0.1.0 commit verified; official example uses TF and TG only"},
    {"method": "IMCell", "order": 1, "candidate_prior": "official mouse_chipX.rda / chipX_gs",
     "selected": False, "source": str(archive), "source_sha256": archive_sha,
     "reason": "Only TF and TG supplied; signed edge types and positive weights cannot be verified"},
    {"method": "IMCell", "order": 2, "candidate_prior": "frozen VIPER mouse CollecTRI, Zenodo 8192729 / HCOP",
     "selected": True, "source": str(signed_path), "source_sha256": sha(signed_path),
     "reason": "Frozen mouse mapped regulator/target, nonzero signed MOR converted to type and positive unit weight"},
]
write_tsv(ROOT / "config/PRIOR_DECISION.tsv", decisions)
audit = {"chipx_rows": len(edges), "chipx_mapped_rows": len(mapped), "chipx_mapped_unique_pairs": distinct,
         "duplicate_mapped_pairs": len(mapped) - distinct, "duplicate_policy": "archive order and multiplicity retained",
         "chipx_lost_rows": dict(lost), "raw_chipx_fields": ["TF", "TG"], "IMCell_edges": len(signed),
         "IMCell_nodes": len(nodes), "IMCell_positive": sum(r["type"] == 1 for r in signed),
         "IMCell_negative": sum(r["type"] == -1 for r in signed), "IMCell_exclusions": dict(disc),
         "IMCell_conflicting_pairs": len(conflicts), "IMCell_candidate_nodes": sum(r["in_network"] for r in coverage),
         "IMCell_eligible_seed_n": sum(r["out_degree"] > 0 for r in coverage),
         "program_counts": {task: {direction: {"before_graph_mapping": sum(r["target"] == task and r["direction"] == direction for r in programs),
                           "after_graph_mapping": sum(r["target"] == task and r["direction"] == direction and r["mapped_to_graph"] for r in programs)}
                      for direction in ("activate", "repress")} for task in TASKS}}
write_json(output / "NETWORK_AND_PROGRAM_AUDIT.json", audit)
print(audit, flush=True)
