"""Draw exact top-ten hit counts and recovered TF ranks from the locked evaluation."""
import argparse
import csv
import hashlib
import io
import json
from math import radians, tau
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

METHODS = {
    "TRIAGE": "TRIAGE",
    "VIPER": "VIPER",
    "Taiji": "Taiji†",
    "ANANSE": "ANANSE†",
    "IMCell": "IMCell\nsingleton",
    "CompForce RRA": "CompForce\nRRA",
}

# Unmodified fill swatches sampled from the user's Figure 1.
TARGETS = {
    "liver": ("Liver", "#CAE4EA", "#4E8195", 300, (2725, 2959)),
    "cardiomyocyte": ("Cardiomyocytes‡", "#F0E2B6", "#A48441", 180, (451, 2594)),
    "betapanc": ("β-TC6", "#AD90D0", "#74559A", 60, (3603, 2431)),
}
FROZEN_TABLES = {
    "HITS_AT_10_AND_STATUS_21.tsv": "65a16a1f7fa12d827b6875d5bcec4a789709b6a3861727c1c9a81d9503544a96",
    "REFERENCE_RANKS_LONG.tsv": "9f278d7eac02c797a4954e591e849f8917629b08ac35bfc4681e9b478df433d9",
}
PALETTE_SHA256 = "33eeb4ca010f660f4177f66527862c92d66fa02c2f81af8643bc375e861b8fd7"

INK = "#172932"
MUTED = "#515E6C"
GRID = "#D9E1E6"
RANK_BAND = "#F4F7F9"
AXES = (0.069, 0.098, 0.862, 0.766)
HUB_R = 33.5
FIRST_RANK_R = 40
LAST_RANK_R = 72
COUNT_BASE = 83
COUNT_STEP = 12
QC_TRACK_BOTTOM = 76
QC_TRACK_TOP = 81
QC_HATCH_STYLE = {
    "facecolor": "#F5F6F7",
    "edgecolor": "#555F67",
    "hatch": "////",
    "linewidth": 0.45,
}
RADIAL_LIMIT = 157
SECTOR = 60
SECTOR_BOUNDARIES = tuple(index * SECTOR + SECTOR // 2 for index in range(len(METHODS)))
RANK_LABELS = (1, 5, 10)
BAR_WIDTH = 10.4
BAR_OFFSETS = (-13.0, 0, 13.0)
STEM = "Figure_6"


class RankRangeLegend(HandlerBase):
    def create_artists(self, legend, handle, x, y, width, height, fontsize, transform):
        del legend, handle, fontsize
        center = width / 2 - x
        bottom = -y
        top = height - y
        cap = width * 0.23
        return [
            Line2D([center, center], [bottom, top], color=INK, lw=0.9, transform=transform),
            Line2D([center - cap, center + cap], [bottom, bottom], color=INK, lw=0.9, transform=transform),
            Line2D([center - cap, center + cap], [top, top], color=INK, lw=0.9, transform=transform),
        ]


def frozen_tsv(directory, name):
    data = (directory / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != FROZEN_TABLES[name]:
        raise ValueError(f"The input does not match the frozen Round 2 table: {name}")
    return csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t")


def load_cells(tables):
    cells = {}
    for row in frozen_tsv(tables, "HITS_AT_10_AND_STATUS_21.tsv"):
        key = row["method"], row["target"]
        if key[0] not in METHODS or key[1] not in TARGETS:
            continue
        if key in cells:
            raise ValueError(f"Duplicate evaluation row: {key}")
        cells[key] = {
            "count": int(row["Hits_at_10"]),
            "reference_n": int(row["reference_total_n"]),
            "picks": int(row["actual_recommendations_n"]),
            "qc": row["interpretation_status"],
            "input_qc": row["input_QC_status"],
            "qc_limitations": row["QC_limitations"],
            "hits": [],
        }

    expected = {(method, target) for method in METHODS for target in TARGETS}
    if cells.keys() != expected:
        raise ValueError(f"Incomplete design: {expected - cells.keys()}")

    for row in frozen_tsv(tables, "REFERENCE_RANKS_LONG.tsv"):
        key = row["method"], row["target"]
        if key in cells and row["hits_exact_top10"] == "True":
            rank = int(row["eval_rank"])
            if not 1 <= rank <= 10:
                raise ValueError(f"Invalid rank for historical hit: {key}/{rank}")
            cells[key]["hits"].append((rank, row["stable_id"], row["reference_symbol"]))

    for key, cell in cells.items():
        cell["hits"].sort()
        cell["ranks"] = [hit[0] for hit in cell["hits"]]
        if len(set(cell["ranks"])) != len(cell["ranks"]) or len(cell["ranks"]) != cell["count"]:
            raise ValueError(f"Inconsistent hit count or evaluated ranks: {key}")
        if not 0 <= cell["count"] <= min(cell["reference_n"], cell["picks"], 10):
            raise ValueError(f"Invalid count: {key}")
    return cells


def write_sources(output, tables, cells):
    with (output / f"{STEM}_DATA.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow([
            "method", "target", "hits_at_10", "hit_ranks", "hit_symbols_in_rank_order",
            "hit_stable_ids_in_rank_order", "min_hit_rank", "median_hit_rank", "max_hit_rank",
            "reference_total_n", "actual_recommendations_n", "interpretation_status",
        ])
        for method in METHODS:
            for target in TARGETS:
                cell = cells[method, target]
                ranks = cell["ranks"]
                low, middle, high = (ranks[0], median(ranks), ranks[-1]) if ranks else ("NA", "NA", "NA")
                writer.writerow([
                    method, target, cell["count"], ",".join(map(str, ranks)),
                    ",".join(hit[2] for hit in cell["hits"]), ",".join(hit[1] for hit in cell["hits"]),
                    low, middle, high, cell["reference_n"], cell["picks"], cell["qc"],
                ])
    provenance = {
        "source_tables": FROZEN_TABLES,
        "bar": "Exact stable-gene count of historical reference TFs in the first 10 recommendations.",
        "rank_marks": "Per method and target, min, median and max of locked evaluation ranks among exact hits. Dots mark each individual hit rank.",
        "rank_guides": {
            "dashed_at_ranks": list(range(1, 11)),
            "labels": RANK_LABELS,
            "label_axes_clockwise_from_top_degrees": SECTOR_BOUNDARIES,
        },
        "qc_hatch": {
            "rule": "Hatch marks input_QC_status == conditional.",
            "angular_width_degrees": BAR_WIDTH,
            "radial_interval": [QC_TRACK_BOTTOM, QC_TRACK_TOP],
            "restricted_cells": [
                {"method": method, "target": target, "limitation": cell["qc_limitations"]}
                for (method, target), cell in cells.items()
                if cell["input_qc"] == "conditional"
            ],
        },
        "palette_source": "figure 1 (1).png, 6519 × 5697 pixels, supplied in chat",
        "palette_source_sha256": PALETTE_SHA256,
        "palette": {
            key: {"fill_hex": fill, "stroke_hex": edge, "sample_xy": point}
            for key, (_, fill, edge, _, point) in TARGETS.items()
        },
        "matplotlib_version": matplotlib.__version__,
    }
    (output / f"{STEM}_SOURCE.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def arc(ax, angle, radius, half_width, **style):
    angles = [angle + half_width * (part - 12) / 12 for part in range(25)]
    return ax.plot(angles, [radius] * len(angles), **style)[0]


def rank_radius(rank):
    return FIRST_RANK_R + (rank - 1) * (LAST_RANK_R - FIRST_RANK_R) / 9


def draw_ranks(ax, key, angle, edge, ranks):
    if not ranks:
        return
    if len(ranks) > 1:
        low, high = rank_radius(ranks[0]), rank_radius(ranks[-1])
        span = ax.plot([angle, angle], [low, high], color=edge, lw=1.0, zorder=4)[0]
        span.set_gid(f"hit-range-{key}")
        for end in (low, high):
            arc(ax, angle, end, radians(2.05), color=edge, lw=1.0, zorder=4)
    center = arc(ax, angle, rank_radius(median(ranks)), radians(2.75), color=edge, lw=2.05, zorder=5)
    center.set_gid(f"hit-median-{key}")
    for rank in ranks:
        dot = ax.plot(angle + radians(3.85), rank_radius(rank), marker="o", ms=3.4,
                      mfc="white", mec=edge, mew=0.9, linestyle="none", zorder=6)[0]
        dot.set_gid(f"hit-dot-{key}-{rank}")


def draw_axes(ax, max_count):
    count_top = COUNT_BASE + COUNT_STEP * max_count
    ax.bar(0, LAST_RANK_R - FIRST_RANK_R, width=tau, bottom=FIRST_RANK_R,
           color=RANK_BAND, linewidth=0, zorder=0)
    for index, axis in enumerate(SECTOR_BOUNDARIES):
        if index % 2 == 0:
            ax.bar(radians(index * SECTOR), count_top - COUNT_BASE, width=radians(SECTOR),
                   bottom=COUNT_BASE, color="#F7F9FA", linewidth=0, zorder=0)
        boundary = radians(axis)
        ax.plot([boundary, boundary], [FIRST_RANK_R, LAST_RANK_R], color=GRID, lw=0.67, zorder=1)
        ax.plot([boundary, boundary], [COUNT_BASE, count_top], color=GRID, lw=0.78, zorder=1)

    angles = [tau * i / 720 for i in range(721)]
    rank_tick_style = {"facecolor": RANK_BAND, "edgecolor": "none", "pad": 0.40}
    for rank in range(1, 11):
        radius = rank_radius(rank)
        major = rank in RANK_LABELS
        guide = ax.plot(
            angles, [radius] * len(angles), color="#B8C8D1" if major else "#CCD9DF",
            lw=0.67 if major else 0.51, linestyle=(0, (2.2, 2.3)) if major else (0, (1.6, 2.5)),
            zorder=1,
        )[0]
        guide.set_gid(f"rank-guide-{rank:02d}")
        if major:
            for axis in SECTOR_BOUNDARIES:
                label = ax.text(radians(axis), radius, str(rank), ha="center", va="center",
                                fontsize=7.1, color=MUTED, bbox=rank_tick_style, zorder=7)
                label.set_gid(f"rank-scale-{axis:03d}-{rank:02d}")

    tick_style = {"facecolor": "white", "edgecolor": "none", "pad": 0.6}
    for hit_n in range(max_count + 1):
        radius = COUNT_BASE + COUNT_STEP * hit_n
        ax.plot(angles, [radius] * len(angles), color="#A9BAC4" if hit_n == 0 else GRID,
                lw=0.85 if hit_n == 0 else 0.66, zorder=1)
        ax.text(radians(90), radius, str(hit_n), ha="center", va="center",
                fontsize=7.3, color=MUTED, bbox=tick_style, zorder=7)
    ax.text(radians(90), count_top + 7, "hits", ha="left", va="center",
            fontsize=7.2, color=MUTED, zorder=7)


def main(tables, output):
    cells = load_cells(tables)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "text.color": INK,
        "figure.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "hatch.linewidth": 0.55,
    })
    fig = plt.figure(figsize=(8.25, 8.25))
    ax = fig.add_axes(AXES, projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, RADIAL_LIMIT)
    ax.set_axis_off()

    fig.text(0.5, 0.967, "Historical TF recovery in the top 10",
             ha="center", va="center", fontsize=20, weight="bold")
    fig.text(0.5, 0.93, "Outer bars count hits. Hatch marks flag input QC limits.",
             ha="center", va="center", fontsize=9.9, color=MUTED)
    range_key = Line2D([], [], color=INK, label="Rank range")
    median_key = Line2D([], [], lw=2.05, color=INK, label="Median rank")
    dot_key = Line2D([], [], marker="o", ms=4, mfc="white", mec=INK, linestyle="none", label="Individual hit rank")
    qc_key = Patch(**QC_HATCH_STYLE, label="Input QC limitation")
    fig.legend(
        handles=[range_key, median_key, dot_key, qc_key], loc="upper center", bbox_to_anchor=(0.5, 0.905),
        handler_map={range_key: RankRangeLegend()}, ncol=4, frameon=False,
        columnspacing=1.75, handletextpad=0.60, handlelength=1.5, fontsize=8.45,
    )

    draw_axes(ax, max(cell["count"] for cell in cells.values()))
    for target, (label, fill, _, angle, _) in TARGETS.items():
        sector = ax.bar(radians(angle), HUB_R, width=radians(120), bottom=0,
                        color=fill, edgecolor="white", linewidth=1.2, zorder=3)[0]
        sector.set_gid(f"hub-sector-{target}")
        label_radius = 17.8 if target == "cardiomyocyte" else 18.4
        text = ax.text(radians(angle), label_radius, label, ha="center", va="center",
                       fontsize=7.85, weight="bold", zorder=4)
        text.set_gid(f"hub-label-{target}")

    background = {"facecolor": "white", "edgecolor": "none", "alpha": 0.91, "pad": 0.25}
    for index, (method, title) in enumerate(METHODS.items()):
        method_theta = index * SECTOR
        method_id = method.lower().replace(" ", "_")
        ax.text(radians(method_theta), RADIAL_LIMIT - 7, title, ha="center", va="center",
                fontsize=11.75, weight="bold", linespacing=1.04, clip_on=False, zorder=7)

        for delta, (target, (_, fill, edge, _, _)) in zip(BAR_OFFSETS, TARGETS.items()):
            angle = radians(method_theta + delta)
            cell = cells[method, target]
            key = f"{method_id}--{target}"
            count = cell["count"]
            if cell["input_qc"] == "conditional":
                qc_strip = ax.bar(angle, QC_TRACK_TOP - QC_TRACK_BOTTOM, bottom=QC_TRACK_BOTTOM,
                                  width=radians(BAR_WIDTH), zorder=3, **QC_HATCH_STYLE)[0]
                qc_strip.set_gid(f"qc-hatch-{key}")
            if count:
                glyph = ax.bar(angle, COUNT_STEP * count, bottom=COUNT_BASE, width=radians(BAR_WIDTH),
                               color=fill, edgecolor=edge, linewidth=0.92, zorder=3)[0]
            else:
                glyph = arc(ax, angle, COUNT_BASE, radians(BAR_WIDTH * 0.48), color=edge,
                            lw=1.8, solid_capstyle="round", zorder=3)
            glyph.set_gid(f"count-bar-{key}")
            label_r = COUNT_BASE + COUNT_STEP * count + 5
            count_text = ax.text(angle, label_r, str(count), ha="center", va="center",
                                 fontsize=8.6, weight="bold" if count else "normal",
                                 bbox=background, zorder=7)
            count_text.set_gid(f"count-text-{key}")
            draw_ranks(ax, key, angle, edge, cell["ranks"])

    fig.text(0.5, 0.073, "Inner dashed guides mark ranks 1–10. Rank 1 is the highest priority.",
             ha="center", fontsize=8.05, color=MUTED)
    fig.text(0.5, 0.051, "† Shared MEF ATAC signal limits.    ‡ Low cardiomyocyte RNA mapping and ATAC signal.",
             ha="center", fontsize=7.55, color=MUTED)

    write_sources(output, tables, cells)
    fig.savefig(output / f"{STEM}.svg", metadata={"Date": None})
    fig.savefig(output / f"{STEM}.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(output / f"{STEM}.png", dpi=320)
    plt.close(fig)
    return cells


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("locked_tables", type=Path)
    p.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = p.parse_args()
    main(args.locked_tables, args.output)
