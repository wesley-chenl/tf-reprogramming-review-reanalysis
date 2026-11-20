"""Render overlapping, rounded ridges from the frozen rank bins."""
from argparse import ArgumentParser
import csv
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.text import Text
import numpy as np
from figure_fonts import system_font

PINK, TEAL = "#E3999A", "#7DB1BC"
PINK_TEXT, TEAL_TEXT = "#99474F", "#356B79"
DISPLAY_HALF_WIDTH = .23
OUTLINE, MUTED = "#221F1F", "#555555"
POSITIONS = {
    "A": (.445, .381, .902, "Shared panel", "GSE53951"),
    "B": (.264, .121, .428, "Taiji-reprogram comparison", "GSE43916"),
    "C": (.083, .121, .247, "TFcomb comparison", "GSE53951"),
}
LABEL = {"tfcomb":"TFcomb scoring module\n(bulk adaptation)"}

def read(file):
    with file.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f,delimiter="\t"))

def rank_x(rank,n):
    return 1-2*(float(rank)-1)/(n-1)

def style():
    paths=[system_font(p) for p in ("times.ttf","timesbd.ttf","timesi.ttf","timesbi.ttf")]
    family="DejaVu Serif"
    if all(paths):
        for path in paths:
            font_manager.fontManager.addfont(path)
        family=font_manager.FontProperties(fname=paths[0]).get_name()
    mpl.rcParams.update({"font.family":family,"font.size":9,"text.color":OUTLINE,"axes.labelcolor":OUTLINE,
                         "svg.fonttype":"none","pdf.fonttype":42,"pdf.use14corefonts":False,
                         "axes.unicode_minus":True,"savefig.facecolor":"white"})

def rounded_heights(centers, mass, grid):
    """Draw each bin's mass as a rounded, bounded contribution for display."""
    distance = (grid[:, None] - centers[None, :]) / DISPLAY_HALF_WIDTH
    support = np.abs(distance) < 1
    waves = np.where(support, .5 + .5*np.cos(np.pi*distance), 0.)
    # Close the silhouette smoothly at the physical endpoints of the rank axis.
    edge = np.clip((1-np.abs(grid))/.085, 0, 1)
    waves *= (.5-.5*np.cos(np.pi*edge))[:, None]
    areas = np.trapezoid(waves, grid, axis=0)
    assert np.all(areas > 0)
    ridge = (waves/areas) @ mass
    assert np.isclose(np.trapezoid(ridge, grid), 1) and np.min(ridge) >= 0
    return ridge


def main():
    parser=ArgumentParser(description=__doc__)
    parser.add_argument("--tables",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    summary=read(args.tables/"PANEL_METHOD_SUMMARY.tsv")
    hist=read(args.tables/"PLOT_BINS.tsv")
    style()
    fig=plt.figure(figsize=(178/25.4,221/25.4),facecolor="white")
    left,width=.244,.675
    center=left+width/2
    scale_min,scale_max=-1.49,1.49
    gutter=1.24
    to_fig=lambda x:left+width*(x-scale_min)/(scale_max-scale_min)
    col_auc=.952
    fig.text(center,.970,"Where screened TF pairs rank",ha="center",va="center",size=17,weight="bold")
    fig.text(.020,.842,"Scoring rule",ha="left",va="bottom",size=8.5,weight="bold")
    for x,name,ink in ((-gutter,"Bottom 10\nnegative recall",PINK_TEXT),(gutter,"Top 10\npositive recall",TEAL_TEXT)):
        fig.text(to_fig(x),.838,name,ha="center",va="bottom",size=8.4,weight="bold",color=ink,linespacing=.95)
    fig.text(col_auc,.842,"AUROC",ha="center",va="bottom",size=8.4,weight="bold")
    for pane,(bottom,height,title_y,title,endpoint) in POSITIONS.items():
        data=[r for r in summary if r["panel"]==pane]
        n,k=int(data[0]["N"]),int(data[0]["K"])
        assert all((int(r["N"]),int(r["K"]))==(n,k) for r in data)
        fig.text(.020,title_y,pane.lower(),ha="left",va="center",size=12,weight="bold")
        fig.text(center,title_y,f"{title} · {endpoint}",ha="center",va="center",size=12.5,weight="bold")
        fig.text(center,title_y-.024,f"{n} pairs   ·   {k} screen positive   ·   {n-k} screen negative",
                 ha="center",va="center",size=9,color=MUTED)
        ax=fig.add_axes([left,bottom,width,height],facecolor="white")
        ax.set(xlim=(scale_min,scale_max),ylim=(0,len(data)))
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_linewidth(1.05)
            side.set_color(OUTLINE)
        ax.axvline(0,linestyle=(0,(4,3)),color="#888888",lw=.85,zorder=1)
        if pane=="C":
            ax.set_xticks([-1,0,1],["Lower priority","Middle","Higher priority"],fontsize=8.7)
            ax.tick_params(axis="x",length=4.5,width=.85,direction="in",pad=5)
        else:
            ax.set_xticks([])
        row_scale=1/len(data)
        for i,r in enumerate(data):
            y=len(data)-i-.76
            row_label_y=y+.18
            yf=bottom+height*row_label_y*row_scale
            name=LABEL.get(r["method"],r["display_name"])
            fig.text(.020,yf,name,size=9.3,ha="left",va="center",weight="bold",linespacing=1.02)
            selected=[h for h in hist if h["panel"]==pane and h["method"]==r["method"]]
            assert len(selected)==40
            grid=np.linspace(-1,1,1001)
            traces=[]
            for cls,fill,shade,med in (
                (0,PINK,PINK_TEXT,"negative_median_midrank"),
                (1,TEAL,TEAL_TEXT,"positive_median_midrank"),
            ):
                h=sorted((v for v in selected if int(v["screen_positive"])==cls),key=lambda v:int(v["bin"]))
                weights=np.array([float(v["class_mass"]) for v in h])
                edges=np.array([float(v["lo"]) for v in h]+[float(h[-1]["hi"])])
                assert np.isclose(weights.sum(),1) and np.allclose(np.diff(edges),.1)
                shape=rounded_heights((edges[:-1]+edges[1:])/2,weights,grid)
                traces.append((cls,fill,shade,med,shape))
            scale=.64/max(shape.max() for _,_,_,_,shape in traces)
            for cls,fill,shade,med,shape in traces:
                curve = scale*shape
                gid=f"{pane}-{r['method']}-{'positive' if cls else 'negative'}"
                ridge=ax.fill_between(grid,y,y+curve,facecolor=fill,alpha=.76,linewidth=0,zorder=2)
                ridge.set_gid(gid)
                ax.plot(grid,np.where(curve>1e-5,y+curve,np.nan),c=shade,lw=.67,
                        ls="-" if cls else (0,(2.4,1.7)),alpha=.9,zorder=3)
                med_x=rank_x(r[med],n)
                top=y+.085+.04*cls
                ax.plot([med_x,med_x],[y,top],color="white",lw=2,zorder=4,solid_capstyle="butt")
                ax.plot([med_x,med_x],[y,top],color=shade,lw=1,zorder=5,solid_capstyle="butt")
            for x,key,col in ((-gutter,"bottom10_negative",PINK_TEXT),(gutter,"top10_positive",TEAL_TEXT)):
                recall=float(r[key+"_recall_pct"])
                star="*" if int(r[key+"_variable"]) else ""
                ax.text(x,row_label_y,f"{recall:.1f}%{star}",ha="center",va="center",size=10.1,color=col)
            fig.text(col_auc,yf,f'{float(r["AUROC"]):.3f}',ha="center",va="center",size=9.3)
    legend=[
        Patch(facecolor=PINK,edgecolor=PINK_TEXT,linewidth=.65,linestyle=(0,(2.4,1.7)),label="Screen negative"),
        Patch(facecolor=TEAL,edgecolor=TEAL_TEXT,linewidth=.5,label="Screen positive"),
    ]
    fig.legend(handles=legend,ncol=2,loc="center",bbox_to_anchor=(center,.037),frameon=False,
               fontsize=8.6,handlelength=1.3,handleheight=1.05,columnspacing=2.2)
    fig.text(center,.016,"Rounded rank-bin profiles. Ticks mark medians. * Cutoff count varies within a tie.",
             ha="center",va="center",size=7.6,color=MUTED)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for label in fig.findobj(Text):
        if label.get_visible() and label.get_text():
            box=label.get_window_extent(renderer)
            assert box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.x1+1 and box.y1<=fig.bbox.y1+1,label.get_text()
    for ext in ("svg","pdf","png"):
        fig.savefig(args.out/f"PAIR_RANK_DISTRIBUTIONS_ROUNDED.{ext}",dpi=300)
    plt.close(fig)

if __name__=="__main__":
    main()
