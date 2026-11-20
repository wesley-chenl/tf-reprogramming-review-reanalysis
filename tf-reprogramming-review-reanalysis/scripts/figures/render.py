"""Rebuild the three manuscript figures and check them against the paper tables."""

import argparse
from collections import Counter
from decimal import Decimal
from math import isclose
from pathlib import Path
from xml.etree import ElementTree

import matplotlib
from PIL import Image

import plot_top10_count_and_rank as figure6
import render_budget10_panel_c_bolder as figure4
import render_mirrored_larger as figure5
from render_rounded import read


ROOT = Path(__file__).resolve().parents[2]
FIGURE_DATA = ROOT / "data/figures"
PAPER_TABLES = ROOT / "data/pairs/manuscript"


def check_figure4(ranks, values):
    expected = read(PAPER_TABLES / "budget10_first12.tsv")
    observed = [
        (
            row["pair_slot_id"], row["bhlh_mouse_symbol"], row["partner_mouse_symbol"],
            row["eval_rank"], row["binary_call"], Decimal(row["q_gated_pct"]),
        )
        for row in ranks
    ]
    manuscript = [
        (
            row["slot"], row["bhlh_symbol"], row["partner_symbol"], int(row["native_rank"]),
            int(Decimal(row["binary"])), Decimal(row["gated"]),
        )
        for row in expected
    ]
    if observed != manuscript:
        raise ValueError("Figure 4 nominations differ from the manuscript table.")
    if (values["N"], values["hits"], values["mean"], values["rank11"]) != (
        533, 3, Decimal("0.011"), Decimal("1.60"),
    ):
        raise ValueError("Figure 4 eligibility, hits, or yield check failed.")
    summary = read(PAPER_TABLES / "budget10_figure_summary.tsv")[0]
    keys = {
        "N": "eligible", "K": "positives", "hits": "selected_positive",
        "uniform_hits": "random_expected_hits", "mean": "mean_yield_pct",
        "maximum": "best_selected_pct", "uniform_mean": "uniform_mean_yield_pct",
        "rank11": "rank11_yield_pct", "rank11_delta_mean": "rank11_minus_mean_pp",
        "rank11_delta_max": "rank11_minus_best_pp", "oracle": "panel_max_pct",
        "regret": "regret_pp",
    }
    for computed, published in keys.items():
        if not isclose(float(values[computed]), float(summary[published]), abs_tol=1e-12):
            raise ValueError(f"Figure 4 variable {computed} differs from the manuscript.")
    print("Figure 4 PASS: 533 pairs, 3 top-ten positives, 0.011% mean, 1.60% at rank 11")


def check_figure5(summary):
    keys = {
        "N": "N", "K": "K", "negative_median_midrank": "neg_median",
        "positive_median_midrank": "pos_median", "AUROC": "auroc",
        "bhlh_weighted_concordance": "within_bhlh",
        "bhlh_positive_negative_comparisons": "within_contrasts",
        "top10_positive_frozen_count": "hits10", "top20_positive_frozen_count": "hits20",
        "bottom10_negative_frozen_count": "negative_bottom10",
        "top10_positive_fixed_bhlh_expected_count": "quota10",
        "top20_positive_fixed_bhlh_expected_count": "quota20",
    }
    original = {(row["panel"], row["method"]): row for row in read(PAPER_TABLES / "mirrored_panel_metrics.tsv")}
    observed = {(row["panel"], row["method"]): row for row in summary}
    if original.keys() != observed.keys() or len(observed) != 10:
        raise ValueError("Figure 5 panel/method roster differs from the manuscript.")
    for key, row in observed.items():
        for computed, published in keys.items():
            if not isclose(float(row[computed]), float(original[key][published]), abs_tol=1e-12):
                raise ValueError(f"Figure 5 metric {computed} differs from the manuscript for {key}.")
    panels = {(row["panel"], int(row["N"]), int(row["K"])) for row in summary}
    if panels != {("A", 210, 37), ("B", 90, 7), ("C", 65, 18)}:
        raise ValueError("Figure 5 panel sizes or positive counts differ from the manuscript.")
    print("Figure 5 PASS: 210, 90 and 65 pairs; class medians and tie-aware AUROC agree")


def check_figure6(cells, output):
    chart_data = (output / "Figure_6_DATA.tsv").read_bytes()
    expected = (FIGURE_DATA / "TOP10_HIT_COUNT_AND_RANK_DATA.tsv").read_bytes()
    if chart_data != expected:
        raise ValueError("Figure 6 hits, ranks or interpretation status differ from the final data table.")
    restricted = {key for key, cell in cells.items() if cell["input_qc"] == "conditional"}
    if len(cells) != 18 or Counter(target for _, target in restricted) != {
        "cardiomyocyte": 6, "liver": 2, "betapanc": 2,
    }:
        raise ValueError("Figure 6 QC cell counts differ from the final design.")
    zero_hits = {key for key in restricted if cells[key]["count"] == 0}
    if zero_hits != {("IMCell", "cardiomyocyte"), ("ANANSE", "betapanc")}:
        raise ValueError("The two zero-hit QC cells differ from the final design.")
    ids = [
        node.get("id", "")
        for node in ElementTree.parse(output / "Figure_6.svg").iter()
        if node.get("id", "").startswith("qc-hatch-")
    ]
    expected_ids = {f"qc-hatch-{m.lower().replace(' ', '_')}--{t}" for m, t in restricted}
    if len(ids) != 10 or set(ids) != expected_ids:
        raise ValueError("Figure 6 vector glyphs differ from the frozen QC annotations.")
    print("Figure 6 PASS: 18 cells, 10 QC hatches, 2 zero-hit QC cells")


def check_exports(output):
    for number in (4, 5, 6):
        stem = output / f"Figure_{number}"
        ElementTree.parse(stem.with_suffix(".svg"))
        if not stem.with_suffix(".pdf").read_bytes().startswith(b"%PDF-"):
            raise ValueError(f"Missing PDF export for Figure {number}.")
        with Image.open(stem.with_suffix(".png")) as png:
            if min(png.info["dpi"]) < 299:
                raise ValueError(f"PNG export for Figure {number} is below 300 dpi.")
            png.load()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results/figures")
    parser.add_argument(
        "--rank-tables", type=Path, default=FIGURE_DATA,
        help="Directory with rebuilt PLOT_SLOTS.tsv and PANEL_METHOD_SUMMARY.tsv for Figure 5.",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    with matplotlib.rc_context():
        ranks, values = figure4.main(["--output", str(output)])
    check_figure4(ranks, values)
    with matplotlib.rc_context():
        summary = figure5.main(["--tables", str(args.rank_tables), "--out", str(output)])
    check_figure5(summary)
    with matplotlib.rc_context():
        cells = figure6.main(ROOT / "data/single", output)
    check_figure6(cells, output)
    check_exports(output)


if __name__ == "__main__":
    main()
