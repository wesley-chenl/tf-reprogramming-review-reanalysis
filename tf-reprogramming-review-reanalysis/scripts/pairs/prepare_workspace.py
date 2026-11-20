"""Stage copied inputs for an isolated score, then independently stage labels."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data/pairs"
CONFIG = REPO / "config/pairs"
CODE = REPO / "scripts/pairs"
LABELS = ("design_pairs.tsv", "outcomes_private.tsv", "constructs.tsv", "mapping_conflicts.tsv")
PRIMARY_INPUTS = {
    "design_public": ("CLEAN_PAIRS.tsv", "SCORING_POLICY.json"),
    "inputs_scoring": ("GSE53951_EFFECT.tsv", "GSE114809_EFFECT.tsv", "TRRUST_MOUSE_SIGNED.tsv"),
}


def sha(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def put(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() and sha(src) != sha(dst):
        raise ValueError(f"Existing byte mismatch in replay workspace {dst}")
    if not dst.exists():
        shutil.copyfile(src, dst)


def score(root, handoff):
    for folder, names in PRIMARY_INPUTS.items():
        for name in names:
            put(handoff / folder / name, root / folder / name)
    for name in ("score.py", "common.py", "run_scoring.sh", "seal.py", "evaluate.py"):
        put(CODE / "utility" / name, root / "code" / name)
    put(CONFIG / "utility/protocol.yaml", root / "protocol.yaml")
    (root / "locks").mkdir(exist_ok=True)
    (root / "evaluation").mkdir(exist_ok=True)


def labels(root, labels_dir):
    receipt = root / "locks/SCORE_LOCK.json"
    if not receipt.is_file():
        raise ValueError(f"Seal scoring before placing outcomes in {root}")
    for name, expected in json.loads(receipt.read_text(encoding="utf-8"))["files"].items():
        if sha(root / name) != expected:
            raise ValueError(f"Changed since score freeze: {name}")
    for name in LABELS:
        put(labels_dir / name, root / "sources_private" / name)


def initial(root, old):
    receipt = json.loads((old / "locks/SCORE_LOCK.json").read_text(encoding="utf-8"))
    for folder, names in PRIMARY_INPUTS.items():
        for name in names:
            artifact = old / folder / name
            assert sha(artifact) == receipt["files"][f"{folder}/{name}"]
            put(artifact, root / "inputs" / name)
    for name in ("PROTOCOL.yaml", "FIGURE_POLICY.yaml", "INPUTS_BY_METHOD.tsv", "METHOD_REGISTRY.tsv"):
        put(CONFIG / "initial" / name, root / name)
    for name in ("BLIND_DESIGN_AUDIT.tsv", "FVS_SELF_LOOPS.tsv"):
        put(CONFIG / "initial" / name, root / "validation" / name)
    for path in (CODE / "initial").iterdir():
        if path.is_file():
            put(path, root / "code" / path.name)
    document = {
        "policy": {n: sha(root / n) for n in ("PROTOCOL.yaml", "FIGURE_POLICY.yaml", "INPUTS_BY_METHOD.tsv", "METHOD_REGISTRY.tsv")},
        "old_scores_sha256": receipt["files"]["scores/all_pair_scores.tsv"],
    }
    (root / "locks").mkdir(exist_ok=True)
    registration = root / "locks/REGISTRATION.json"
    if registration.exists() and json.loads(registration.read_text(encoding="utf-8")) != document:
        raise ValueError(f"Registration differs from {registration}")
    registration.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("phase", choices=("score", "labels", "initial"))
    ap.add_argument("workspace", type=Path)
    ap.add_argument("--source", type=Path)
    args = ap.parse_args()
    if args.phase == "score":
        score(args.workspace, args.source or DATA / "score_inputs")
    elif args.phase == "labels":
        labels(args.workspace, args.source or DATA / "labels")
    elif args.phase == "initial":
        if args.source is None:
            ap.error("--source must identify the sealed utility staging workspace")
        initial(args.workspace, args.source)
