"""Seal independently generated scores and positional labels before joins."""
import argparse
from pathlib import Path
import json
from common import now, read, record, sha

SCORE = [
    "protocol.yaml", "code/score.py", "code/run_scoring.sh", "design_public/SCORING_POLICY.json",
    "design_public/CLEAN_PAIRS.tsv", "scores/all_pair_scores.tsv", "scores/top_blind.tsv", "scores/program.tsv",
    "scores/READ_ACCESS.json", "scores/SCORING_INPUT_SHA256.txt", "scores/SCORER_SHA256.txt", "scores/ISOLATION.txt",
]
LABEL = [
    "protocol.yaml", "sources_private/design_pairs.tsv", "sources_private/outcomes_private.tsv",
    "sources_private/mapping_conflicts.tsv", "sources_private/constructs.tsv",
]


def seal(root, side):
    output = root / "locks" / (side + "_LOCK.json")
    if output.exists():
        from common import check_lock
        check_lock(root, side)
        return
    paths = list(SCORE if side == "SCORE" else LABEL)
    if side == "SCORE":
        for line in (root / "scores/SCORING_INPUT_SHA256.txt").read_text().splitlines():
            digest, name = line.split(maxsplit=1)
            name = name.lstrip("*")
            folder = "design_public" if name in {"CLEAN_PAIRS.tsv","SCORING_POLICY.json"} else "inputs_scoring"
            relative = f"{folder}/{name}"
            if sha(root / relative) != digest:
                raise ValueError(f"Source bytes changed since scoring: {relative}")
            if relative not in paths:
                paths.append(relative)
        access = json.loads((root / "scores/READ_ACCESS.json").read_text())
        if not any(branch["state"] == "scored" for branch in access["endpoints"]):
            raise ValueError("No expression endpoint actually scored")
    else:
        supplemental = tuple("sources_private/"+name for name in (
            "SOURCE_OUTCOMES.tsv","OUTCOME_AUDIT.json","figure_grid_audit.tsv","figure_ocr_audit.tsv",
            "s1_control_cells.tsv","s1_raw_outcomes.tsv","ORTHOLOGY.tsv"))
        paths.extend(p for p in supplemental if (root / p).is_file())
        lines = read(root / "sources_private/outcomes_private.tsv")
        if any(r["binary_call"] not in {"", "0", "1", "NA"} for r in lines if r["source_version"] == "nature2018"):
            raise ValueError("Invalid binary label")
    record(output, {"sealed_utc": now(), "side": side, "files": {p:sha(root / p) for p in paths}})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("side", choices=["SCORE","LABEL"])
    args = p.parse_args()
    seal(args.root, args.side)
