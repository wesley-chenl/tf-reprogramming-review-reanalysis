"""Rebuild the gated-yield rank panel in bold, larger type from frozen results."""
import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]

from figure_fonts import system_font
from plot_section5_2_budget10 import D, EDGE, GRID, INK, SECONDARY, digest, emit_tsv, get_data, palette

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.text import Text
from matplotlib.ticker import FormatStrFormatter


STEM = "Figure_4"
WIDTH_MM, HEIGHT_MM = 189, 89


def fonts():
    bold = system_font("ARIALNBI.TTF")
    if bold:
        fontManager.addfont(bold)
        return FontProperties(fname=bold, family="Arial Narrow", style="italic", weight="bold", size=11.6)
    return FontProperties(family="DejaVu Sans", style="italic", weight="bold", size=10.7)


def bars(ax, rows, colors, prefix):
    for row in rows:
        rank = row["eval_rank"]
        color = "blue" if rank <= 10 else "green"
        if row["binary_call"]:
            mark = ax.bar(rank, float(row["q_gated_pct"]), width=.56, color=colors[color],
                          edgecolor=EDGE[color], lw=.85, zorder=4)[0]
            kind = "bar"
        else:
            mark, = ax.plot(rank, 0, "o", ms=4.2, mfc="white", mec=EDGE[color], mew=.94,
                            clip_on=False, zorder=5)
            kind = "zero"
        mark.set_gid(f"{prefix}-rank-{rank}-{kind}")


def plot(output, ranks, values, colors):
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "axes.titleweight": "bold",
        "text.color": INK,
        "axes.labelcolor": INK,
        "xtick.color": SECONDARY,
        "ytick.color": SECONDARY,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.facecolor": "white",
        "axes.linewidth": .7,
    })
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=300)
    full = fig.add_axes((.089, .238, .903, .738))
    full.set(xlim=(.68, 12.61), ylim=(0, 1.83), yticks=(0, .4, .8, 1.2, 1.6))
    full.set_ylabel("Gated yield (%)", fontsize=13.5, labelpad=4.5)
    full.set_xticks(range(1, 13), [str(i) for i in range(1, 13)])
    full.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    full.tick_params(axis="x", labelsize=11.3, length=2.8, pad=2.9)
    full.tick_params(axis="y", labelsize=11.45, length=2.5, pad=2.8)
    for side in ("top", "right"):
        full.spines[side].set_visible(False)
    full.axvspan(10.5, 12.61, color=colors["green"], alpha=.23, zorder=0).set_gid("adjacent-outside-budget")
    full.grid(axis="y", color=GRID, linewidth=.52, zorder=1)
    full.axvline(10.5, color="#667D70", ls=(0, (3, 3)), lw=.93, zorder=2).set_gid("budget-boundary")
    full.text(10.30, 1.70, "Budget = 10", fontsize=11.35, color="#365649", ha="right").set_gid("budget-label")

    bars(full, ranks, colors, "full")
    gene_font = fonts()
    gene_artists = {"bhlh": [], "partner": []}
    for row in ranks:
        rank = row["eval_rank"]
        for y, key, role in ((-.120, "bhlh_mouse_symbol", "bhlh"),
                             (-.205, "partner_mouse_symbol", "partner")):
            artist = full.text(rank, y, row[key], transform=full.get_xaxis_transform(), ha="center",
                               va="top", fontproperties=gene_font, clip_on=False)
            artist.set_gid(f"rank-{rank}-{role}")
            gene_artists[role].append(artist)

    y11 = float(values["rank11"])
    full.text(11, y11 + .061, f'{values["rank11"]:.2f}%', color=EDGE["green"], ha="center",
              va="bottom", fontsize=13.25).set_gid("rank11-yield-label")
    full.hlines(y11, 11.28, 12.61, color=EDGE["green"], linestyles=(0, (2, 2)), linewidth=.72, zorder=3)
    refs = [
        ("mean", 11.63, "yellow", (0, (3, 2))),
        ("maximum", 12.24, "purple", (0, (1, 1.5))),
    ]
    for tag, x, color, linestyle in refs:
        y = float(values[tag])
        full.hlines(y, x - .22, x + .20, color=EDGE[color], linestyles=linestyle, linewidth=.78, zorder=6)
        arrow = full.annotate("", xy=(x, y11), xytext=(x, y), arrowprops={
            "arrowstyle": "<->", "color": EDGE[color], "lw": .82, "shrinkA": 0, "shrinkB": 0,
            "mutation_scale": 7,
        }, zorder=6)
        arrow.arrow_patch.set_gid(f"local-delta-{tag}")
        delta = values["rank11"] - values[tag]
        full.annotate(f"{delta:.3f} pp", (x, (y11 + y) / 2), xytext=(4.25, 0),
                      textcoords="offset points", rotation=90, ha="left", va="center", fontsize=11.65,
                      color=EDGE[color], zorder=7).set_gid(f"rank11-delta-{tag}")

    zoom = full.inset_axes((.124, .218, .621, .603), zorder=8)
    zoom.set_facecolor("#FCFEFE")
    zoom.set(xlim=(.48, 10.53), ylim=(0, .08), xticks=range(1, 11), yticks=(0, .04, .08))
    zoom.set_title("Enlarged view · ranks 1–10", loc="left", fontsize=12.3, pad=3.9)
    zoom.set_ylabel("Yield (%)", fontsize=11.25, labelpad=2.8)
    zoom.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    zoom.tick_params(axis="both", labelsize=10.3, direction="out", length=2, pad=2.3)
    for spine in zoom.spines.values():
        spine.set_color("#9BACB4")
        spine.set_linewidth(.65)
    for tag, _, color, linestyle in refs:
        zoom.axhline(float(values[tag]), color=EDGE[color], ls=linestyle, lw=.92, zorder=2).set_gid(f"zoom-{tag}-rule")
    bars(zoom, ranks[:10], colors, "zoom")
    for row in ranks[:10]:
        if row["binary_call"]:
            rank = row["eval_rank"]
            y = D(row["q_gated_pct"])
            zoom.text(rank, float(y) + .0027, f"{y:.2f}%", color=EDGE["blue"], fontsize=10.45,
                      ha="center", va="bottom", zorder=7).set_gid(f"zoom-yield-{rank}")
    for x, y, tag, name, color in ((5.9, .0127, "mean", "Mean", "yellow"),
                                    (5.9, .0720, "maximum", "Maximum", "purple")):
        zoom.text(x, y, f"{name} {values[tag]:.3f}%", va="bottom", ha="center", fontsize=10.85,
                  color=EDGE[color], zorder=7, bbox={"facecolor": "#FCFEFE", "edgecolor": "none", "pad": .15}
                  ).set_gid(f"zoom-{tag}-label")

    labels = [text for text in fig.findobj(Text) if text.get_visible() and text.get_text()]
    for label in labels:
        label.set_fontweight("bold")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    if gene_font.get_file() is None:
        # A wider fallback font needs to fit the same rank slots.
        width_limit = full.transData.transform((2, 0))[0] - full.transData.transform((1, 0))[0]
        width_limit -= 2.0 * fig.dpi / 72
        widest = max(
            (left.width + right.width) / 2
            for artists in gene_artists.values()
            for boxes in [[artist.get_window_extent(renderer) for artist in artists]]
            for left, right in zip(boxes, boxes[1:])
        )
        if widest > width_limit:
            gene_font.set_size(gene_font.get_size_in_points() * width_limit / widest)
            for artists in gene_artists.values():
                for artist in artists:
                    artist.set_fontproperties(gene_font)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
    minimum_gaps = {}
    for role, artists in gene_artists.items():
        boxes = [artist.get_window_extent(renderer) for artist in artists]
        minimum_gap = min(right.x0 - left.x1 for left, right in zip(boxes, boxes[1:])) * 72 / fig.dpi
        if minimum_gap < 1:
            raise ValueError(f"{role} labels have only {minimum_gap:.2f} pt separation")
        minimum_gaps[role] = minimum_gap
    label_boxes = []
    for label in labels:
        box = label.get_window_extent(renderer)
        if box.x0 < 0 or box.x1 > fig.bbox.x1 or box.y0 < 0 or box.y1 > fig.bbox.y1:
            raise ValueError(f"Text extends outside the canvas: {label.get_text()}")
        for other_label, other_box in label_boxes:
            if box.overlaps(other_box):
                raise ValueError(f"Overlapping labels: {other_label!r}, {label.get_text()!r}")
        label_boxes.append((label.get_text(), box))
    for extension in ("svg", "pdf", "png"):
        options = {"dpi": 300} if extension == "png" else {
            "metadata": {"Date": None} if extension == "svg" else {"CreationDate": None, "ModDate": None}
        }
        fig.savefig(output / f"{STEM}.{extension}", **options)
    plt.close(fig)
    return minimum_gaps, gene_font.get_name(), gene_font.get_size_in_points()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=REPO / "data/pairs")
    parser.add_argument("--output", type=Path, default=REPO / "results/figures")
    args = parser.parse_args(argv)
    ranks, values, provenance = get_data(args.source.resolve())
    colors, _ = palette(None)
    args.output.mkdir(parents=True, exist_ok=True)
    original_data = REPO / "data/figures/SECTION5_2_BUDGET10_TOP12.tsv"
    plotted_data = args.output / f"{STEM}_TOP12.tsv"
    emit_tsv(plotted_data, ranks)
    if digest(original_data.read_bytes()) != digest(plotted_data.read_bytes()):
        raise ValueError("The rebuilt rank table differs from the frozen figure input.")
    gaps, gene_face, gene_pt = plot(args.output, ranks, values, colors)
    receipt = dict(provenance, paired_rank_count=len(ranks), original_rank_tsv_sha256=digest(original_data.read_bytes()),
                   figure_mm=[WIDTH_MM, HEIGHT_MM], gene_font=gene_face, gene_font_pt=gene_pt,
                   minimum_gene_gap_pt=gaps, blue_positive_ranks=[r["eval_rank"] for r in ranks[:10] if r["binary_call"]],
                   exact_top_ten_mean_pct=str(values["mean"]), exact_rank11_pct=str(values["rank11"]),
                   output_sha256={ext: digest((args.output / f"{STEM}.{ext}").read_bytes()) for ext in ("svg", "pdf", "png")})
    (args.output / "Figure_4_render.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return ranks, values


if __name__ == "__main__":
    main()
