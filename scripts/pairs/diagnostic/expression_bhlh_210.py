"""Reconstruct the post hoc expression ablation on the fixed 210-pair GSE53951 panel."""
import csv
import gzip
from collections import Counter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent
ENDPOINT = "GSE53951"


def records(path):
    open_input = gzip.open if path.suffix == ".gz" else open
    with open_input(path, "rt", encoding="utf-8", newline="") as source:
        yield from csv.DictReader(source, delimiter="\t")


def primary_expression(row):
    return (row["endpoint"], row["variant"], row["method"]) == (ENDPOINT, "main", "expression")


def write_table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def score_pairs(eligible_slots, expression_scores):
    membership = {
        row["pair_id"]: row["slot"] for row in records(eligible_slots)
    }
    ranked = []
    for row in records(expression_scores):
        if not primary_expression(row) or row["pair_id"] not in membership:
            continue
        score = Fraction(row["bhlh_p"])
        ranked.append(dict(slot=membership[row["pair_id"]], pair_id=row["pair_id"],
                           bhlh=row["bhlh_symbol"], partner=row["partner_symbol"],
                           score=score, tie_key=row["tie_key"]))
    assert len(ranked) == len(membership) == len({row["slot"] for row in ranked}) == 210
    assert list(Counter(row["bhlh"] for row in ranked).values()) == [10] * 21
    ranked.sort(key=lambda row: (-row["score"], row["tie_key"], row["pair_id"]))
    for rank, row in enumerate(ranked, 1):
        row["rank"] = rank
    return ranked


def summarize(ranked, screened_slots):
    labels = {
        row["slot"]: int(row["screen_positive"]) for row in records(screened_slots)
        if row.get("panel") == "A" and row.get("method") == "expression"
    }
    positives = [row for row in ranked if labels[row["slot"]] == 1]
    negatives = [row for row in ranked if labels[row["slot"]] == 0]
    assert (len(positives), len(negatives)) == (37, 173)
    global_twice_wins = 0
    within_twice_wins = 0
    group_sizes = Counter()
    for pos in positives:
        for neg in negatives:
            twice_credit = 2 * (pos["score"] > neg["score"]) + (pos["score"] == neg["score"])
            global_twice_wins += twice_credit
            if pos["bhlh"] == neg["bhlh"]:
                group_sizes[pos["bhlh"]] += 1
                within_twice_wins += twice_credit
    global_pairs = len(positives) * len(negatives)
    within_pairs = sum(group_sizes.values())
    return dict(N=len(ranked), positive=len(positives), informative_bhlh=len(group_sizes),
                informative_within_comparisons=within_pairs,
                AUROC=global_twice_wins / (2 * global_pairs),
                within_concordance=within_twice_wins / (2 * within_pairs),
                hits10=sum(labels[row["slot"]] for row in ranked[:10]),
                hits20=sum(labels[row["slot"]] for row in ranked[:20]))


if __name__ == "__main__":
    from argparse import ArgumentParser
    from hashlib import sha256
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--eligibility", required=True, type=Path)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    ranking = score_pairs(args.eligibility, args.scores)
    rank_file = args.out / "expression_bhlh_210_ranks.tsv"
    write_table(rank_file, ranking)
    (args.out / "EXPRESSION_DIAGNOSTIC_SCORE.sha256").write_text(
        sha256(rank_file.read_bytes()).hexdigest() + "  " + rank_file.name + "\n", encoding="ascii"
    )
    write_table(args.out / "expression_bhlh_210_summary.tsv", [summarize(ranking, args.labels)])
