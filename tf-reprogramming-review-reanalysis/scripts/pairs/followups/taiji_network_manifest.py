"""Check and inventory the real source/target graphs before producing pair scores."""
import csv
import hashlib
import math
from pathlib import Path
import sys

network = Path(sys.argv[1]) / "Network"
dest = Path(sys.argv[2])
rows = []
node_sets = {}
for state in ("source","target"):
    folder = network / state
    nodes_path = folder / "nodes.csv"
    with nodes_path.open(newline="") as nodesfile:
        nodes = list(csv.DictReader(nodesfile))
    identifiers = {row["geneName:ID"] for row in nodes}
    assert identifiers and len(identifiers) == len(nodes), f"{state} incomplete/duplicate native nodes"
    node_sets[state] = identifiers
    for name in ("nodes.csv","edges_binding.csv","edges_combined.csv"):
        path = folder / name
        if name=="nodes.csv":
            records=len(nodes)
        else:
            edge_column="weight" if name=="edges_combined.csv" else "affinity"
            with path.open(newline="") as stream:
                records=0
                for edge in csv.DictReader(stream):
                    assert math.isfinite(float(edge[edge_column])), f"{state} nonfinite {edge_column}"
                    if name=="edges_combined.csv":
                        assert edge[":START_ID"] and edge[":END_ID"], (state,edge)
                    records+=1
        assert records>0, f"{state} empty {name}"
        with path.open("rb") as source:
            digest=hashlib.file_digest(source,"sha256").hexdigest()
        rows.append({"state":state,"type":name,"path":str(path),"records":records,
                     "bytes":path.stat().st_size,"sha256":digest})
with dest.open("w",newline="") as file:
    w=csv.DictWriter(file,list(rows[0]),delimiter="\t",lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print({state:len(ids) for state,ids in node_sets.items()}, "shared", len(node_sets["source"] & node_sets["target"]))
