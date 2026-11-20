"""Extract the manuscript diagnostic's missing inputs from frozen upstream artifacts."""
import argparse
import csv
import gzip
import io
import json
from fractions import Fraction
from pathlib import Path

from rebuild_rank_panels import rows


def gzip_table(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(records)


def extract(data, out):
    key = "GSE53951", "main", "expression"
    expression = {r["pair_id"]: r for r in rows(data / "utility/scores/all_pair_scores.tsv")
                  if (r["endpoint"], r["variant"], r["method"]) == key and r["scoreable"] == "1"}
    design = {r["pair_id"]: r for r in rows(data / "labels/design_pairs.tsv")
              if r["variant"] == "main" and r["eligible"] in {"1", "True", "true"}}
    coefficient_table = data / "initial/tfcomb_primary_adapter.tsv"
    score_inputs = data / "score_inputs/inputs_scoring"
    effects = sorted(
        (dict(stable_id=r["stable_id"], delta=r["d_g"]) for r in rows(score_inputs / "GSE53951_EFFECT.tsv")
         if r["d_g"] not in {"", "NA", "NaN", "Inf", "-Inf"}),
        key=lambda r: r["stable_id"],
    )
    genes = {r["stable_id"] for r in effects}
    eligible = []
    for r in rows(coefficient_table):
        pair = r["pair_id"]
        if not r["score"] or pair not in expression or pair not in design:
            continue
        e = expression[pair]
        c = json.loads(r["details"])
        eligible.append(dict(
            slot=design[pair]["pair_slot_id"], bhlh_symbol=e["bhlh_symbol"], partner_symbol=e["partner_symbol"],
            bhlh_id=e["bhlh_id"].removeprefix("GeneID:"), partner_id=e["partner_id"].removeprefix("GeneID:"),
            bhlh_coef=c["bhlh_coef"], partner_coef=c["partner_coef"],
            tfcomb_pair_score=r["score"], bhlh_fraction=str(Fraction(e["bhlh_p"])),
            expression_pair_fraction=str(Fraction(int(e["num"]), int(e["den"]))), tie_key=e["tie_key"],
        ))
    assert len(eligible) == len({r["slot"] for r in eligible}) == 65
    assert len(genes) == len(effects) == 20478
    members = {r[key] for r in eligible for key in ("bhlh_id", "partner_id")}
    edges = [dict(tf_gene_id=r["tf_gene_id"], target_gene_id=r["target_gene_id"], sign=r["sign"])
             for r in rows(score_inputs / "TRRUST_MOUSE_SIGNED.tsv")
             if r["tf_gene_id"] in members and r["target_gene_id"] in genes]
    inputs = out / "score_inputs"
    gzip_table(inputs / "panel_inputs.tsv.gz", sorted(eligible, key=lambda r: r["slot"]))
    gzip_table(inputs / "gene_effect.tsv.gz", effects)
    gzip_table(inputs / "direct_edges.tsv.gz", edges)
    return {r["slot"] for r in eligible}


def extract_labels(data, out, slots):
    labels = [dict(slot=r["pair_slot_id"], positive=r["binary_call"], gated_yield_pct=r["q_gated"])
              for r in rows(data / "labels/outcomes_private.tsv")
              if r["source_version"] == "nature2018" and r["pair_slot_id"] in slots]
    assert len(labels) == 65 and sum(int(r["positive"]) for r in labels) == 18
    gzip_table(out / "labels/screen_outcomes.tsv.gz", sorted(labels, key=lambda r: r["slot"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    slots = extract(args.data, args.out)
    extract_labels(args.data, args.out, slots)
