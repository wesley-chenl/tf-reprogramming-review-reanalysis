"""Rebuild manuscript tables and Figures 4–6 from the bundled frozen inputs."""

import argparse
import csv
import hashlib
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def check_inputs():
    with (ROOT / "INPUT_SHA256.tsv").open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            path = (ROOT / row["path"]).resolve()
            if not path.is_relative_to(ROOT):
                raise ValueError(f"Input path leaves the repository: {row['path']}")
            with path.open("rb") as source:
                digest = hashlib.file_digest(source, "sha256").hexdigest()
            if digest != row["sha256"]:
                raise ValueError(f"Input checksum differs: {row['path']}")


def run(title, script, *args):
    print(f"\n{title}", flush=True)
    subprocess.run(
        [sys.executable, "-B", str(ROOT / script), *map(str, args)],
        cwd=ROOT,
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build")
    args = parser.parse_args()
    output = args.output.resolve()
    pairs = output / "tables" / "pairs"
    single = output / "tables" / "single"
    figures = output / "figures"
    check_inputs()
    run("Recompute paired-TF tables", "scripts/pairs/reproduce.py", "--out", pairs)
    run("Recompute single-TF tables", "scripts/single/recompute.py", "--output", single)
    run("Render Figures 4–6", "scripts/figures/render.py", "--output", figures, "--rank-tables", pairs)
    print(f"\nReproduction complete in {output}")


if __name__ == "__main__":
    main()
