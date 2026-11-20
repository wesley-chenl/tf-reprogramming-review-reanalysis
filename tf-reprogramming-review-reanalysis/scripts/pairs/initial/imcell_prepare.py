"""Create exact endpoint programs and canonical seed keys from isolated unlabeled inputs."""
import csv
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import sys

inp, out = map(Path, sys.argv[1:3])
def get(name):
    with (inp / name).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def put(name, rows, names):
    with (out / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, names, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def identifier(value):
    return value.removeprefix("GeneID:") if value.startswith("GeneID:") else ""

policy = json.loads((inp / "SCORING_POLICY.json").read_text(encoding="utf-8-sig"))
design = get("CLEAN_PAIRS.tsv")
assert len(design) == 598 * 9
assert len({(x["variant"], x["pair_id"]) for x in design}) == len(design)
edges = get("TRRUST_MOUSE_SIGNED.tsv")
assert len(edges) == 3964
assert {x["sign"] for x in edges} == {"1", "-1"}
assert len({(x["tf_gene_id"], x["target_gene_id"]) for x in edges}) == len(edges)
node_ids = {x[k] for x in edges for k in ("tf_gene_id", "target_gene_id")}
origins = {x["tf_gene_id"] for x in edges}
assert sum(x["tf_gene_id"] == x["target_gene_id"] for x in edges) == 16
panel = {identifier(x[k]) for x in design for k in ("bhlh_id", "partner_id")} - {""}
panel_symbols = {x[k] for x in design for k in ("bhlh_symbol", "partner_symbol")} - {""}

program = []
for endpoint in policy["endpoint_order"]:
    source = get(endpoint + "_EFFECT.tsv")
    candidate = []
    for x in source:
        if x["stable_id"] not in panel and x["symbol"] not in panel_symbols and x["d_g"] not in {"", "NA", "NaN", "Inf", "-Inf"}:
            candidate.append((x["stable_id"], Fraction(x["d_g"])))
    positive = sorted((x for x in candidate if x[1] > 0), key=lambda x:(-x[1], x[0]))
    negative = sorted((x for x in candidate if x[1] < 0), key=lambda x:(x[1], x[0]))
    n = policy["program_each_direction"]
    for direction, rows in ((1, positive[:n]), (-1, negative[:n])):
        for stable_id, _ in rows:
            program.append({"endpoint":endpoint, "stable_id":stable_id, "direction":direction,
                            "graph_available":int(stable_id in node_ids)})
put("imcell_target_program.tsv", program, ("endpoint", "stable_id", "direction", "graph_available"))

pairs = []
valid = set()
for x in design:
    a,b = (identifier(x[z]) for z in ("bhlh_id","partner_id"))
    if not a or not b:
        why = "unresolved_construct_identity"
    elif a == b:
        why = "distinct_two_seeds_required"
    elif a not in node_ids or b not in node_ids:
        why = "seed_absent_from_signed_graph"
    elif a not in origins or b not in origins:
        why = "seed_has_no_outgoing_edge_in_signed_graph"
    else:
        why = ""
    left,right = sorted((a,b))
    work_id = "|".join((left,right))
    if not why:
        valid.add((work_id,left,right))
    row = {k:x[k] for k in ("variant","pair_id","bhlh_id","partner_id")}
    row["tie_key"] = sha256("|".join((policy["tie_salt"], x["bhlh_id"],x["partner_id"])).encode()).hexdigest()
    pairs.append({**row,"work_id":work_id,"reason":why})
put("imcell_expand.tsv", pairs, ("variant","pair_id","bhlh_id","partner_id","work_id","reason","tie_key"))
jobs = [{"work_id":wid, "left":a, "right":b,
         "audit_key":sha256(("section4-imcell-audit-v1|" + wid).encode()).hexdigest()}
        for wid,a,b in sorted(valid)]
put("imcell_jobs.tsv", jobs, ("work_id","left","right","audit_key"))
print(f"{len(edges)} signed edges; {len(node_ids)} graph nodes; {len(jobs)} distinct eligible unordered pairs")
