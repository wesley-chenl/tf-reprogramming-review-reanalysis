"""Rebuild the §5.2 budget-ten figure from the published Tsunemoto result files."""
import argparse
from collections import Counter
from contextlib import contextmanager
import csv
from decimal import Decimal
from hashlib import sha256
import io
import json
from pathlib import Path
from zipfile import ZipFile

D = Decimal
REPO = Path(__file__).resolve().parents[2]
PANEL = {"endpoint": "GSE53951", "variant": "main", "panel_type": "native", "method": "expression"}
LOCATORS = (
    "evaluation/top_pairs.tsv",
    "evaluation/utility_results.tsv",
    "sources_private/outcomes_private.tsv",
    "sources_private/design_pairs.tsv",
    "scores/all_pair_scores.tsv",
)
ARCHIVE = "TSUNEMOTO_PAIR_UTILITY_LIGHT.zip"
ARCHIVE_SHA256 = "c1e1d9ece936b9ea7dab15d190d6864f84f389a8ca646a51a69e584b809a3f79"
PUBLIC_INPUTS = {
    "evaluation/top_pairs.tsv": "utility/evaluation/top_pairs.tsv",
    "evaluation/utility_results.tsv": "utility/evaluation/utility_results.tsv",
    "sources_private/outcomes_private.tsv": "labels/outcomes_private.tsv",
    "sources_private/design_pairs.tsv": "labels/design_pairs.tsv",
    "scores/all_pair_scores.tsv": "utility/scores/all_pair_scores.tsv",
}
STEM = "SECTION5_2_BUDGET10"
INK = "#182932"
SECONDARY = "#435967"
GRID = "#DFE7EA"
PALETTE_SAMPLES = {"blue": (2725, 2959), "yellow": (451, 2594), "purple": (3603, 2431), "green": (4381, 3059)}
PALETTE_VALUES = {"blue": "#CAE4EA", "yellow": "#F0E2B6", "purple": "#AD90D0", "green": "#BDDFB6"}
EDGE = {"blue": "#3F7688", "yellow": "#8B691F", "purple": "#735397", "green": "#457A50"}


def tsv(data):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig")), delimiter="\t"))


def emit_tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest(data):
    return sha256(data).hexdigest()


@contextmanager
def frozen_reader(path):
    if all((path / name).is_file() for name in PUBLIC_INPUTS.values()):
        manifest = REPO / "config/figures/fig4_source_hashes.tsv"

        def read_public(name):
            return manifest.read_bytes() if name == "RESULT_MANIFEST.tsv" else (path / PUBLIC_INPUTS[name]).read_bytes()

        yield read_public, "repository data/pairs"
        return
    if path.is_dir() and all((path / name).is_file() for name in ("RESULT_MANIFEST.tsv", *LOCATORS)):
        yield lambda name: (path / name).read_bytes(), str(path)
        return
    archive = path if path.is_file() else (path / ARCHIVE if (path / ARCHIVE).is_file() else path.parent / ARCHIVE)
    if digest(archive.read_bytes()) != ARCHIVE_SHA256:
        raise ValueError("Archive SHA-256 differs from the frozen deliverable.")
    with ZipFile(archive) as source:
        yield source.read, str(archive)


def select(rows, **filters):
    return [row for row in rows if all(row.get(k) == v for k, v in filters.items())]


def exactly_one(rows):
    if len(rows) != 1:
        raise ValueError(f"Expected one result row, found {len(rows)}.")
    return rows[0]


def agree(actual, expected, label, tol=D("0")):
    if isinstance(actual, Decimal) and expected is not None:
        matches = abs(actual - D(expected)) <= tol
    else:
        matches = actual == expected
    if not matches:
        raise ValueError(f"{label} mismatch: calculated {actual}, comparison {expected}")


def get_data(source):
    with frozen_reader(source) as (read_bytes, locator):
        manifest = {r["path"]: r for r in tsv(read_bytes("RESULT_MANIFEST.tsv"))}
        data = {}
        hashes = {}
        for path in LOCATORS:
            raw = read_bytes(path)
            hashes[path] = digest(raw)
            if hashes[path] != manifest[path]["sha256"] or len(raw) != int(manifest[path]["bytes"]):
                raise ValueError(f"Source does not match the published manifest: {path}")
            data[path] = tsv(raw)

    top_all = data[LOCATORS[0]]
    b20 = sorted(select(top_all, **PANEL, utility="gated", B="20"), key=lambda r: int(r["eval_rank"]))
    agree(len(b20), 20, "B=20 result size")
    agree([int(r["eval_rank"]) for r in b20], list(range(1, 21)), "Unique evaluation ranks")
    agree(len({r["pair_id"] for r in b20}), 20, "Unique ranked pairs")
    for utility in ("binary", "gated"):
        short = sorted(select(top_all, **PANEL, utility=utility, B="10"), key=lambda r: int(r["eval_rank"]))
        agree([r["pair_id"] for r in short], [r["pair_id"] for r in b20[:10]], f"{utility} budget prefix")
    binary_stat = exactly_one(select(data[LOCATORS[1]], **PANEL, utility="binary", B="10"))
    gated_stat = exactly_one(select(data[LOCATORS[1]], **PANEL, utility="gated", B="10"))

    labels = {}
    for row in select(data[LOCATORS[2]], source_version="nature2018"):
        slot = row["pair_slot_id"]
        if slot in labels:
            raise ValueError(f"Duplicate measured slot {slot}")
        labels[slot] = row
    identity = {r["pair_id"]: r for r in select(data[LOCATORS[3]], variant=PANEL["variant"])}
    old_scores = select(data[LOCATORS[4]], endpoint=PANEL["endpoint"], variant=PANEL["variant"], method=PANEL["method"], scoreable="1")
    eligible = [r for r in old_scores if identity[r["pair_id"]]["eligible"] == "1"
                and labels[identity[r["pair_id"]]["pair_slot_id"]]["q_gated"] not in ("", "NA")]
    eligible.sort(key=lambda r: int(r["rank"]))
    agree(len({r["pair_id"] for r in eligible}), len(eligible), "Unique native nominations")
    agree([r["pair_id"] for r in eligible[:20]], [r["pair_id"] for r in b20], "Budget order from independent score/mask/label join")
    agree([int(r["rank"]) for r in eligible[:20]], [int(r["global_score_rank"]) for r in b20], "Global scores traced separately")
    pool = [labels[identity[r["pair_id"]]["pair_slot_id"]] for r in eligible]
    n = len(pool)
    k = sum(int(D(r["binary_call"])) for r in pool)
    for row in (binary_stat, gated_stat):
        agree((n, k), (int(row["N"]), int(row["K"])), "Native N and K")
    agree((n,k), (533,69), "Protocol N and K")
    uniform_h = D(10 * k) / n
    uniform_gated = sum(D(r["q_gated"]) for r in pool) / n
    eps = D("1e-13")
    agree(uniform_h, D(binary_stat["null_mean"]) * 10, "Expected binary hits", eps)
    agree(uniform_gated, gated_stat["null_mean"], "Expected mean gated yield", eps)

    ranks = []
    for row in b20[:12]:
        ident = identity[row["pair_id"]]
        y = labels[ident["pair_slot_id"]]
        agree(ident["pair_slot_id"], row["pair_slot_id"], "Physical slot")
        agree((ident["bhlh_symbol"], ident["partner_symbol"]), (row["bhlh"], row["partner"]), "Native gene identities")
        call = int(D(y["binary_call"]))
        gated = D(y["q_gated"])
        raw = None if row["raw_s1_pct"] in ("", "NA") else D(row["raw_s1_pct"])
        use_swap = bool(ident.get("raw_s1_column_if_mapped", "") and y.get("raw_s1_cell_s1_thesis", "").startswith(ident["raw_s1_column_if_mapped"]))
        raw_value = y["raw_s1_pct_s1_thesis"] if use_swap else y["raw_s1_pct"]
        agree(raw, D(raw_value) if raw_value not in ("", "NA") else None, "S1 raw value trace")
        agree(D(row["binary_call"]), call, "Screen call")
        agree(D(row["q_gated"]), gated, "Gated yield")
        agree(D(row["q"]), gated, "Plot metric is q_gated on the gated B20 row")
        agree(gated, D(y["fig1b_pct_printed"]) if call else D(0), "Gate definition")
        agree(row["source_locator"], y["source_cell_or_bbox"], "Source locator")
        ranks.append({
            "endpoint": PANEL["endpoint"], "variant": PANEL["variant"], "panel_type": PANEL["panel_type"],
            "method": PANEL["method"], "source_B": 20, "source_utility": "gated",
            "eval_rank": int(row["eval_rank"]), "global_score_rank": int(row["global_score_rank"]),
            "pair_id": row["pair_id"], "pair_slot_id": row["pair_slot_id"],
            "bhlh_mouse_symbol": row["bhlh"], "partner_mouse_symbol": row["partner"],
            "binary_call": call, "q_gated_pct": str(gated),
            "reported_fig1b_pct": y["fig1b_pct_printed"], "raw_s1_pct": "" if raw is None else str(raw),
            "raw_s1_status": "missing" if raw is None else "measured_zero" if raw==0 else "measured_nonzero",
            "raw_s1_cell": y["raw_s1_cell_s1_thesis"] if use_swap else y["raw_s1_cell"],
            "readout_status": y["measurement_status"], "call_basis": y["call_basis"],
            "readout_day": y["day"], "yield_denominator": y["denominator"],
            "source_locator": y["source_cell_or_bbox"],
            "source_url": "https://www.nature.com/articles/s41586-018-0103-5/figures/1",
        })

    selected = ranks[:10]
    hits = sum(r["binary_call"] for r in selected)
    ys = [D(r["q_gated_pct"]) for r in selected]
    mean, maximum = sum(ys) / 10, max(ys)
    neighbor = D(ranks[10]["q_gated_pct"])
    oracle = max(D(p["q_gated"]) for p in pool)
    regret = oracle-maximum
    agree(hits, int(binary_stat["hit_n"]), "Observed B10 hits")
    agree(D(hits)/10, binary_stat["mean_utility"], "Observed binary mean")
    for name, x, field in (("mean",mean,"mean_utility"),("maximum",maximum,"best_selected"),
                           ("oracle",oracle,"oracle"),("BestRegret",regret,"best_regret")):
        agree(x, gated_stat[field], name, eps)
    agree([(r["eval_rank"],D(r["q_gated_pct"])) for r in selected if r["binary_call"]],
          [(3,D(".01")),(4,D(".03")),(9,D(".07"))], "Expected top-ten hits and yields")
    agree((ranks[10]["bhlh_mouse_symbol"],ranks[10]["partner_mouse_symbol"],neighbor), ("Ascl1","Pou3f2",D("1.60")), "Rank 11")
    agree((ranks[11]["binary_call"],D(ranks[11]["q_gated_pct"])), (0,D("0")), "Rank 12")
    agree(regret,"5.63","Formal BestRegret")
    values = {"N":n,"K":k,"B":10,"hits":hits,"uniform_hits":uniform_h,"mean":mean,"uniform_mean":uniform_gated,
              "maximum":maximum,"rank11":neighbor,"rank11_delta_mean":neighbor-mean,"rank11_delta_max":neighbor-maximum,
              "oracle":oracle,"regret":regret}
    return ranks, values, {"frozen_source":locator,"validated_source_hashes":hashes}


def palette(image_path):
    if image_path is None:
        return PALETTE_VALUES, {"source":"Cached exact swatches from provided figure 1 (2).png","colors":PALETTE_VALUES}
    from PIL import Image
    with Image.open(image_path) as image:
        img = image.convert("RGB")
        if img.size != (6519,5697):
            raise ValueError("Unexpected size for original color-reference image.")
        colors = {}
        for name,(x,y) in PALETTE_SAMPLES.items():
            patch = img.crop((x-4,y-4,x+5,y+5))
            pixels=patch.get_flattened_data() if hasattr(patch,"get_flattened_data") else list(patch.getdata())
            rgb = Counter(pixels).most_common(1)[0][0]
            colors[name] = "#{:02X}{:02X}{:02X}".format(*rgb)
    agree(colors, PALETTE_VALUES, "Exact reference palette")
    return colors, {"source":str(image_path),"sha256":digest(image_path.read_bytes()),"samples_xy":PALETTE_SAMPLES,"colors":colors,"ink_strokes":EDGE}


def save_data(output, ranks, v, provenance):
    emit_tsv(output/f"{STEM}_TOP12.tsv",ranks)
    definitions = [
        ("N","eligible_pairs","Score lock, identity mask and Nature 2018 outcomes"),
        ("K","screen_positive_pairs","binary_call on exactly the same native pool"),
        ("B","nominations","B=10"),
        ("hits","screen_positive_pairs","sum of binary_call in the first 10 eligible ranks"),
        ("uniform_hits","screen_positive_pairs","B*K/N, uniform without replacement"),
        ("mean","percent","arithmetic mean of all 10 gated yields, including zero-valued nominations"),
        ("uniform_mean","percent","native pool mean of q_gated, uniform without replacement"),
        ("maximum","percent","highest q_gated among the first 10 eligible ranks"),
        ("rank11","percent","q_gated at eval_rank 11"),
        ("rank11_delta_mean","percentage_points","rank 11 gated yield minus top-ten mean"),
        ("rank11_delta_max","percentage_points","rank 11 gated yield minus top-ten maximum"),
        ("oracle","percent","largest q_gated on the same 533-pair native pool"),
        ("regret","percentage_points","BestRegret@10, native pool oracle minus top-ten maximum"),
    ]
    emit_tsv(output/f"{STEM}_CHECKS.tsv",[
        {"metric":k,"computed":str(v[k]),"unit":unit,"definition_or_crosscheck":description,"check":"PASS"}
        for k,unit,description in definitions
    ])
    (output/f"{STEM}_PROVENANCE.json").write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    caption = (
        "**Fixed-budget evaluation of a frozen TF-pair ranking.** In the Tsunemoto mouse embryonic fibroblast (MEF) to induced-neuron screen, "
        "an expression rule averaged two role-specific percentiles estimated from the independent GSE53951 endpoint. "
        f"The conservative expression-native panel retained {v['N']} evaluated pairs, of which {v['K']} were screen-positive. "
        f"**a,** The first ten nominations recovered {v['hits']} positives, compared with {v['uniform_hits']:.2f} expected under uniform sampling without replacement from the same panel. "
        f"**b,** Mean hit-gated reported yield was {v['mean']:.4f}% for the ten nominations and {v['uniform_mean']:.4f}% for uniform sampling. "
        "**c,** Bars show hit-gated reported yield in frozen evaluation order. For screen-positive pairs, yield is the published normalized TUJ1-positive percentage, "
        "with the initial number of plated MEFs as denominator. Screen-negative pairs are assigned zero gated yield and shown as hollow circles. "
        "The accompanying data distinguish measured S1 zeros from missing uncorrected values. The inset magnifies ranks 1–10 on a linear scale. "
        f"Rank 11 (*Ascl1* + *Pou3f2*, {v['rank11']:.2f}%) is the adjacent out-of-budget example. "
        f"The two arrows compare its yield with the top-ten mean and maximum, giving local differences of {v['rank11_delta_mean']:.3f} and {v['rank11_delta_max']:.3f} percentage points, respectively. "
        f"The same-panel maximum was {v['oracle']:.2f}%, yielding a formal BestRegret@10 of {v['regret']:.2f} percentage points. "
        "[Tsunemoto et al., Figure 1b](https://www.nature.com/articles/s41586-018-0103-5/figures/1).\n"
    )
    (output/f"{STEM}_CAPTION_EN.md").write_text(caption,encoding="utf-8")


def plot(output, ranks, v, colors):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FormatStrFormatter

    matplotlib.rcParams.update({
        "font.family":"DejaVu Sans","text.color":INK,"axes.labelcolor":INK,
        "xtick.color":SECONDARY,"ytick.color":SECONDARY,"pdf.fonttype":42,"svg.fonttype":"none",
        "figure.facecolor":"white","axes.linewidth":0.65,
    })
    fig=plt.figure(figsize=(178/25.4,151/25.4))
    fig.text(.065,.961,"Fixed-budget evaluation of a TF-pair ranking",fontsize=14.0,weight="bold")
    fig.text(.065,.921,f"Mouse MEF-to-neuron screen  ·  GSE53951 expression ranking  ·  {v['N']} eligible pairs",
             fontsize=8.1,color=SECONDARY)
    legend=[
        Line2D([],[],marker="o",mfc=EDGE["blue"],mec=EDGE["blue"],ms=5.8,ls="none",label="Selected ten"),
        Line2D([],[],marker="|",color=EDGE["yellow"],mew=1.35,ms=8,ls="none",label="Uniform random expectation"),
    ]
    fig.legend(handles=legend,loc="upper left",bbox_to_anchor=(.058,.924),ncol=2,fontsize=7.65,frameon=False,
               columnspacing=1.3,handletextpad=.38,handlelength=1.25)

    def axis_a_b(box,letter,title,obs,rngobs,labels,limits,ticks,unit):
        ax=fig.add_axes(box)
        fig.text(box[0],box[1]+box[3]+.018,letter,fontsize=11,weight="bold",va="bottom")
        fig.text(box[0]+.026,box[1]+box[3]+.021,title,fontsize=9,weight="bold",va="bottom")
        ax.set_xlim(*limits)
        ax.set_ylim(0,1)
        ax.set_yticks([])
        ax.set_xticks(ticks)
        if unit=="Percent":
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        ax.set_xlabel(unit,fontsize=7.35,labelpad=1.7)
        ax.tick_params(axis="x",labelsize=7,length=2.5,pad=2)
        for side in ("left","top","right"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#97A8AE")
        ax.hlines(.48,*limits,color=GRID,lw=.95)
        obs_mark=ax.plot(float(obs),.48,marker="o",mfc=EDGE["blue"],mec=EDGE["blue"],ms=5.6,ls="none",zorder=3)[0]
        obs_mark.set_gid(f"{letter}-observed-mark")
        rand_mark=ax.vlines(float(rngobs),.32,.64,color=EDGE["yellow"],lw=1.5,zorder=3)
        rand_mark.set_gid(f"{letter}-uniform-mark")
        for xx, yy, txt, name, color in (
                (obs,.74,labels[0],"observed",EDGE["blue"]),(rngobs,.74,labels[1],"uniform",EDGE["yellow"])):
            ax.text(float(xx),yy,txt,ha="center",va="bottom",fontsize=8.25,weight="bold",color=color).set_gid(f"{letter}-{name}-label")

    axis_a_b((.088,.757,.388,.069),"a","Screen-positive pairs among ten",v["hits"],v["uniform_hits"],
             (str(v["hits"]),f"{v['uniform_hits']:.2f}"),(0,4),(0,1,2,3,4),"Pairs")
    axis_a_b((.560,.757,.394,.069),"b","Mean gated yield per nomination",v["mean"],v["uniform_mean"],
             (f"{v['mean']:.4f}%",f"{v['uniform_mean']:.4f}%"),(0,.10),(0,.025,.05,.075,.10),"Percent")
    fig.text(.065,.677,"c",fontsize=12,weight="bold")
    fig.text(.091,.679,"Gated yield in frozen rank order",fontsize=10.5,weight="bold")

    full=fig.add_axes((.087,.255,.888,.405))
    full.set(xlim=(.45,12.96),ylim=(0,1.83),yticks=(0,.4,.8,1.2,1.6))
    full.set_ylabel("Gated yield (%)",fontsize=10.2,labelpad=4)
    full.set_xticks(range(1,13),[str(i) for i in range(1,13)])
    full.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    full.tick_params(axis="x",labelsize=8.8,length=2.8,pad=2.2)
    full.tick_params(axis="y",labelsize=8.9,length=2.2,pad=2)
    for side in ("top","right"):
        full.spines[side].set_visible(False)
    full.axvspan(10.5,12.96,color=colors["green"],alpha=.23,zorder=0).set_gid("adjacent-outside-budget")
    full.grid(axis="y",color=GRID,linewidth=.52,zorder=1)
    full.axvline(10.5,color="#667D70",ls=(0,(3,3)),lw=.93,zorder=2).set_gid("budget-boundary")
    full.text(10.30,1.70,"Budget = 10",fontsize=8.8,color="#365649",ha="right",weight="bold").set_gid("budget-label")

    def bars(ax,rows,prefix):
        for r in rows:
            x=r["eval_rank"]
            y=float(r["q_gated_pct"])
            color="blue" if x<=10 else "green"
            if r["binary_call"]:
                shape=ax.bar(x,y,width=.56,color=colors[color],edgecolor=EDGE[color],lw=.85,zorder=4)[0]
            else:
                shape=ax.plot(x,0,"o",ms=4.0,mfc="white",mec=EDGE[color],mew=.88,clip_on=False,zorder=5)[0]
            shape.set_gid(f"{prefix}-rank-{x}-{'bar' if r['binary_call'] else 'zero'}")
    bars(full,ranks,"full")
    tf=full.get_xaxis_transform()
    narrow_path=next((p for p in (Path("C:/Windows/Fonts/ARIALNI.TTF"),
                                  Path("/mnt/c/Windows/Fonts/ARIALNI.TTF")) if p.is_file()),None)
    gene_font=FontProperties(fname=narrow_path, family="Arial Narrow", style="italic", size=10) if narrow_path else FontProperties(
        family="DejaVu Sans", style="italic", size=8.1)
    gene_artists={"bhlh":[],"partner":[]}
    for r in ranks:
        x=r["eval_rank"]
        for y,text,role in ((-.127,r["bhlh_mouse_symbol"],"bhlh"),(-.209,r["partner_mouse_symbol"],"partner")):
            artist=full.text(x,y,text,transform=tf,ha="center",va="top",
                             fontproperties=gene_font,clip_on=False)
            artist.set_gid(f"rank-{x}-{role}")
            gene_artists[role].append(artist)
    fig.text(.531,.127,"Frozen pair rank",fontsize=10,ha="center")

    full.text(11, float(v["rank11"])+.061, f"{v['rank11']:.2f}%", color=EDGE["green"], ha="center",
              va="bottom",fontsize=10.2,weight="bold").set_gid("rank11-yield-label")
    y11=float(v["rank11"])
    full.hlines(y11,11.28,12.61,color=EDGE["green"],linestyles=(0,(2,2)),linewidth=.72,zorder=3)
    refs=[("mean",v["mean"],11.67,"yellow",(0,(3,2)),"rank11_delta_mean"),
          ("maximum",v["maximum"],12.42,"purple",(0,(1,1.5)),"rank11_delta_max")]
    for tag,value,x,color,linestyle,key in refs:
        y=float(value)
        full.hlines(y,x-.22,x+.20,color=EDGE[color],linestyles=linestyle,linewidth=.78,zorder=6)
        arrow=full.annotate("",xy=(x,y11),xytext=(x,y),arrowprops={
            "arrowstyle":"<->","color":EDGE[color],"lw":.78,"shrinkA":0,"shrinkB":0,"mutation_scale":6.5},zorder=6)
        arrow.arrow_patch.set_gid(f"local-delta-{tag}")
        full.annotate(f"{v[key]:.3f} pp",(x,(y11+y)/2),xytext=(4.8,0),textcoords="offset points",
                      rotation=90,ha="left",va="center",fontsize=9.2,weight="bold",
                      color=EDGE[color],zorder=7).set_gid(key)

    zoom=full.inset_axes((.119,.215,.626,.603),zorder=8)
    zoom.set_facecolor("#FCFEFE")
    zoom.set(xlim=(.48,10.53),ylim=(0,.08),xticks=range(1,11),yticks=(0,.04,.08))
    zoom.set_title("Enlarged view  ·  ranks 1–10",loc="left",fontsize=9.5,pad=4.1,weight="bold")
    zoom.set_ylabel("Yield (%)",fontsize=8.5,labelpad=2.1)
    zoom.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    zoom.tick_params(axis="both",labelsize=8.0,direction="out",length=2,pad=2)
    for spine in zoom.spines.values():
        spine.set_color("#9BACB4")
        spine.set_linewidth(.62)
    for tag,value,_,color,style,_ in refs:
        zoom.axhline(float(value),color=EDGE[color],ls=style,lw=.9,zorder=2).set_gid(f"zoom-{tag}-rule")
    bars(zoom,ranks[:10],"zoom")
    for row in ranks[:10]:
        if row["binary_call"]:
            rank=row["eval_rank"]
            y=D(row["q_gated_pct"])
            zoom.text(rank,float(y)+.0027,f"{y:.2f}%",color=EDGE["blue"],fontsize=8,
                      ha="center",va="bottom",weight="bold",zorder=7).set_gid(f"zoom-yield-{rank}")
    for x,y,txt,tag in (
        (5.9,.0127,f"Mean {v['mean']:.3f}%","mean"),(5.9,.0719,f"Maximum {v['maximum']:.3f}%","maximum")):
        zoom.text(x,y,txt,va="bottom",ha="center",fontsize=8.4,
                  color=EDGE["yellow" if tag=="mean" else "purple"],weight="bold",zorder=7,
                  bbox={"facecolor":"#FCFEFE","edgecolor":"none","pad":.15}).set_gid(f"zoom-{tag}-label")

    for x,color,text in ((.087,"yellow","Rank 11 minus top-ten mean"),
                         (.50,"purple","Rank 11 minus top-ten maximum")):
        fig.add_artist(Line2D([x,x+.027],[.082,.082],transform=fig.transFigure,color=EDGE[color],lw=1.5))
        fig.text(x+.038,.082,text,va="center",fontsize=8.4)
    fig.text(.087,.042,"pp = percentage points",fontsize=7.8,color=SECONDARY)

    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for role,artists in gene_artists.items():
        edges=[a.get_window_extent(renderer) for a in artists]
        minimum_gap=min(right.x0-left.x1 for left,right in zip(edges,edges[1:]))*72/fig.dpi
        if minimum_gap < 1:
            raise ValueError(f"Gene labels in the {role} row are separated by only {minimum_gap:.2f} pt")

    for ext in ("svg","pdf","png"):
        kwargs={"dpi":300} if ext=="png" else {"metadata":{"Date":None} if ext=="svg" else {"CreationDate":None,"ModDate":None}}
        fig.savefig(output/f"{STEM}.{ext}",**kwargs)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    here=Path(__file__).resolve().parent
    p.add_argument("--source",type=Path,default=REPO/"data/pairs")
    p.add_argument("--output",type=Path,default=here)
    p.add_argument("--palette-image",type=Path)
    args=p.parse_args()
    ranks, values, metadata=get_data(args.source.resolve())
    colors, info=palette(args.palette_image)
    args.output.mkdir(parents=True,exist_ok=True)
    save_data(args.output,ranks,values,{**metadata,"palette":info,"figure_width_mm":178})
    plot(args.output,ranks,values,colors)
    print(f"Plotted 12 frozen ranks, N={values['N']}, K={values['K']}, B=10.")
