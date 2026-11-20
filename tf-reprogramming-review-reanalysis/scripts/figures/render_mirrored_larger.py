"""Render the frozen mirrored rank distributions in larger bold Arial."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
import csv
from fractions import Fraction
from math import isclose
from pathlib import Path
from statistics import median
from matplotlib import font_manager, rcParams

from figure_fonts import system_font
from render_rounded import (LABEL, MUTED, OUTLINE, PINK, PINK_TEXT, TEAL, TEAL_TEXT,
                            Text, np, plt, read, rounded_heights, style)

STEM="Figure_5"
PANELS={
    "A": (.493,.350,.912,"Shared panel","GSE53951"),
    "B": (.302,.104,.449,"Taiji-reprogram comparison","GSE43916"),
    "C": (.111,.104,.258,"TFcomb comparison","GSE53951"),
}
LEFT, RIGHT=.288,.916
CENTER=(LEFT+RIGHT)/2
AUC_X=.957
NBINS=20

def write_tsv(path, rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def row_data(summary, slots):
    n,k=int(summary["N"]),int(summary["K"])
    assert len(slots)==n and len({r["slot"] for r in slots})==n
    assert sorted(int(r["rank"]) for r in slots)==list(range(1,n+1))
    ranks={0:[],1:[]}
    exported=[]
    for r in slots:
        cls=int(r["screen_positive"])
        exact_midrank=Fraction(int(r["rank_lo"])+int(r["rank_hi"]),2)
        assert exact_midrank==Fraction(r["midrank"]) and 1<=exact_midrank<=n
        ranks[cls].append(exact_midrank)
        sign=1 if cls else -1
        exported.append(dict(panel=r["panel"],endpoint=r["endpoint"],method=r["method"],slot=r["slot"],
                             bhlh_symbol=r["bhlh_symbol"],partner_symbol=r["partner_symbol"],
                             full_panel_n=n,screen_positive=cls,frozen_rank=r["rank"],rank_lo=r["rank_lo"],
                             rank_hi=r["rank_hi"],midrank=float(exact_midrank),signed_plot_rank=float(sign*exact_midrank)))
    assert len(ranks[0])==int(summary["negative_N"])==n-k and len(ranks[1])==k
    neg=[float(x) for x in ranks[0]]
    pos=[float(x) for x in ranks[1]]
    assert median(neg)==float(summary["negative_median_midrank"])
    assert median(pos)==float(summary["positive_median_midrank"])
    auc=sum((p<q)+.5*(p==q) for p in pos for q in neg)/(k*len(neg))
    assert isclose(auc,float(summary["AUROC"]))
    return ranks,exported


def bins_for(ranks,n):
    counts=Counter(min(NBINS-1,int(r*NBINS/n)) for r in ranks)
    probabilities=np.array([counts[i]/len(ranks) for i in range(NBINS)])
    assert np.isclose(probabilities.sum(),1)
    return counts,probabilities


def main(argv=None):
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--tables",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args(argv)
    args.out.mkdir(parents=True,exist_ok=True)
    summary=read(args.tables/"PANEL_METHOD_SUMMARY.tsv")
    source=defaultdict(list)
    for r in read(args.tables/"PLOT_SLOTS.tsv"):
        source[r["panel"],r["method"]].append(r)

    style()
    font=system_font("arialbd.ttf")
    family="DejaVu Sans"
    if font:
        font_manager.fontManager.addfont(font)
        family=font_manager.FontProperties(fname=font).get_name()
    rcParams.update({"font.family":family,
                     "font.weight":"bold","axes.labelweight":"bold"})
    fig=plt.figure(figsize=(178/25.4,225/25.4),facecolor="white")
    fig.text(CENTER,.974,"Full rankings of screened TF pairs",ha="center",va="center",size=20.1,weight="bold")
    fig.text(.018,.854,"Scoring rule",ha="left",va="center",size=10.3,weight="bold")
    fig.text((LEFT+CENTER)/2,.854,"Screen negatives",ha="center",va="center",
             color=PINK_TEXT,size=11.05,weight="bold")
    fig.text((CENTER+RIGHT)/2,.854,"Screen positives",ha="center",va="center",
             color=TEAL_TEXT,size=11.05,weight="bold")
    fig.text(AUC_X,.854,"AUROC",ha="center",va="center",size=9.85,weight="bold")
    grid=np.linspace(-1,1,1001)
    edges=np.linspace(-1,1,NBINS+1)
    centers=(edges[:-1]+edges[1:])/2
    all_slots,all_bins=[],[]
    for pane,(bottom,height,title_y,title,endpoint) in PANELS.items():
        sub=[s for s in summary if s["panel"]==pane]
        n,k=int(sub[0]["N"]),int(sub[0]["K"])
        assert all((int(s["N"]),int(s["K"]),s["endpoint"])==(n,k,endpoint) for s in sub)
        fig.text(.018,title_y,pane.lower(),ha="left",va="center",size=14.35,weight="bold")
        fig.text(CENTER,title_y,f"{title} · {endpoint}",ha="center",va="center",size=14.8,weight="bold")
        fig.text(CENTER,title_y-.023,f"{n} pairs   ·   {n-k} negative   ·   {k} positive",
                 ha="center",va="center",size=10.5,color=MUTED)
        ax=fig.add_axes([LEFT,bottom,RIGHT-LEFT,height],facecolor="white")
        ax.set(xlim=(-1.045*n,1.045*n),ylim=(0,len(sub)))
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(1.02)
            spine.set_color(OUTLINE)
        ax.axvline(0,linestyle=(0,(4,3)),color="#888888",lw=.85,zorder=1)
        step=max(5,(n//15)*5)
        ticks=[-n,-2*step,-step,0,step,2*step,n]
        ax.set_xticks(ticks,[str(abs(t)) for t in ticks],fontsize=9.85,fontweight='bold')
        ax.tick_params(axis="x",length=3.3,width=.8,direction="in",pad=3)
        for index,s in enumerate(sub):
            m=s["method"]
            row_key=f"{pane}-{m}"
            ranks,exported=row_data(s,source[pane,m])
            all_slots.extend(exported)
            base=len(sub)-index-.689
            yfig=bottom+height*(base+.23)/len(sub)
            fig.text(.018,yfig,LABEL.get(m,s["display_name"]),size=11.2,ha="left",
                     va="center",weight="bold",linespacing=1.03)
            shapes={}
            for cls in (0,1):
                counts,probabilities=bins_for(ranks[cls],n)
                shapes[cls]=rounded_heights(centers,probabilities,grid)
                for i,mass in enumerate(probabilities):
                    all_bins.append(dict(panel=pane,method=m,screen_positive=cls,bin=i+1,
                                         outward_rank_lo=i*n/NBINS,outward_rank_hi=(i+1)*n/NBINS,
                                         count=counts[i],class_size=len(ranks[cls]),mass=float(mass)))
            scale=.60/max(shape.max() for shape in shapes.values())
            for cls,colour,ink,side in ((0,PINK,PINK_TEXT,"negative"),(1,TEAL,TEAL_TEXT,"positive")):
                sign=1 if cls else -1
                ridge=ax.fill_between(sign*n*(grid+1)/2,base,base+scale*shapes[cls],
                                      facecolor=colour,linewidth=0,zorder=2)
                ridge.set_gid(f"ridge-{row_key}-{side}")
                mid=float(median(ranks[cls]))
                mark,=ax.plot([sign*mid,sign*mid],[base,base+.14],color=ink,lw=1.05,zorder=4,
                              solid_capstyle="butt")
                mark.set_gid(f"median-tick-{row_key}-{side}")
                tag=ax.text(sign*mid,base-.041,f"{mid:g}",ha="center",va="top",size=9.15,color=ink)
                tag.set_gid(f"median-label-{row_key}-{side}")
            auroc=fig.text(AUC_X,yfig,f'{float(s["AUROC"]):.3f}',ha="center",va="center",size=11.5)
            auroc.set_gid(f"auroc-{row_key}")

    fig.text(CENTER,.056,"Full-panel rank",ha="center",va="center",size=12.55,weight="bold")
    fig.text(CENTER,.033,"On either side, rank 1 is nearest the center. Ticks on ridges mark the medians.",
             ha="center",va="center",size=9.65 if font else 8.15,color=MUTED)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    shown=[]
    for label in fig.findobj(Text):
        if label.get_visible() and label.get_text():
            box=label.get_window_extent(renderer)
            assert box.x0>=-1 and box.x1<=fig.bbox.x1+1 and box.y0>=-1 and box.y1<=fig.bbox.y1+1,label.get_text()
            shown.append((label.get_text(),box))
    for i,(name,box) in enumerate(shown):
        for other_name,other in shown[:i]:
            if box.overlaps(other):
                raise ValueError(f"Text boxes overlap: {name!r}, {other_name!r}")
    assert len(all_slots)==sum(int(s["N"]) for s in summary)
    write_tsv(args.out/"MIRRORED_SLOT_RANKS.tsv",all_slots)
    write_tsv(args.out/"MIRRORED_RANK_BINS.tsv",all_bins)
    for ext in ("svg","pdf","png"):
        metadata = {"Date": None} if ext == "svg" else {"CreationDate": None, "ModDate": None} if ext == "pdf" else None
        fig.savefig(args.out/f"{STEM}.{ext}",dpi=300,metadata=metadata)
    plt.close(fig)
    return summary

if __name__=="__main__":
    main()
