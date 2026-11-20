"""Recompute the summaries printed in the two paired-screen figures. Requires Python 3."""
from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader, DictWriter
from pathlib import Path
from statistics import mean, median
from math import isclose
import gzip

ROOT = Path(__file__).resolve().parents[2]

def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        yield from DictReader(stream,delimiter='\t')

def export(path, rows):
    with path.open('w',newline='',encoding='utf8') as stream:
        writer=DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def ranks(source, output):
    panels=defaultdict(list)
    for row in table(source):
        panels[row['panel'],row['method']].append(row)
    rows=[]
    for (panel,method),group in panels.items():
        n=len(group)
        assert len({x['slot'] for x in group})==n
        assert set(int(x['frozen_rank']) for x in group)==set(range(1,n+1))
        assert all(int(x['full_panel_n'])==n for x in group)
        pos=[x for x in group if int(x['screen_positive'])==1]
        neg=[x for x in group if int(x['screen_positive'])==0]
        assert all(float(x['signed_plot_rank'])==float(x['midrank']) for x in pos)
        assert all(float(x['signed_plot_rank'])==-float(x['midrank']) for x in neg)
        pairs=[(x,y) for x in pos for y in neg]
        within=[(x,y) for x,y in pairs if x['bhlh_symbol']==y['bhlh_symbol']]
        def concordance(subset):
            return mean((float(x['midrank'])<float(y['midrank']))+.5*(float(x['midrank'])==float(y['midrank'])) for x,y in subset)
        def quota(budget):
            means=defaultdict(list)
            for x in group: means[x['bhlh_symbol']].append(int(x['screen_positive']))
            return sum(mean(means[x['bhlh_symbol']]) for x in group if int(x['frozen_rank'])<=budget)
        rows.append(dict(panel=panel,method=method,N=n,K=len(pos),neg_median=median([float(x['midrank']) for x in neg]),
                         pos_median=median([float(x['midrank']) for x in pos]),auroc=concordance(pairs),
                         within_bhlh=concordance(within),within_contrasts=len(within),
                         hits10=sum(int(x['frozen_rank'])<=10 for x in pos),hits20=sum(int(x['frozen_rank'])<=20 for x in pos),
                         negative_bottom10=sum(int(x['frozen_rank'])>n-10 for x in neg),quota10=quota(10),quota20=quota(20)))
    export(output / 'mirrored_panel_metrics.tsv',rows)
    return rows

def budget(source, output):
    grp=sorted((r for r in table(source) if r['variant']=='main' and r['endpoint']=='GSE53951'
                and r['method']=='expression' and r['native_rank']),key=lambda r:float(r['native_rank']))
    n=len(grp)
    k=sum(float(x['binary']) for x in grp)
    assert (n,k)==(533,69)
    assert len({x['slot'] for x in grp})==n
    yield_pct=[float(x['gated']) for x in grp]
    sel=grp[:10]
    results=dict(eligible=n,positives=int(k),selected_positive=int(sum(float(r['binary']) for r in sel)),
                 random_expected_hits=10*k/n,mean_yield_pct=mean(yield_pct[:10]),best_selected_pct=max(yield_pct[:10]),
                 uniform_mean_yield_pct=mean(yield_pct),rank11_yield_pct=yield_pct[10],
                 rank11_minus_mean_pp=yield_pct[10]-mean(yield_pct[:10]),rank11_minus_best_pp=yield_pct[10]-max(yield_pct[:10]),
                 panel_max_pct=max(yield_pct),regret_pp=max(yield_pct)-max(yield_pct[:10]))
    assert isclose(results['regret_pp'],5.63)
    export(output / 'budget10_figure_summary.tsv',[results])
    export(output / 'budget10_first12.tsv',[{z:r[z] for z in ('slot','bhlh_symbol','partner_symbol','native_rank','binary','gated')} for r in grp[:12]])
    return results

if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--mirrored-ranks', type=Path, default=ROOT / 'results/figures/MIRRORED_SLOT_RANKS.tsv')
    parser.add_argument('--initial-ranks', type=Path, default=ROOT / 'data/pairs/ranks/pair_ranks_initial.tsv.gz')
    parser.add_argument('--output', type=Path, default=ROOT / 'results/figure_checks')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    r=ranks(args.mirrored_ranks, args.output)
    b=budget(args.initial_ranks, args.output)
    assert len(r)==10 and sum(x['N'] for x in r)==1570 and b['selected_positive']==3
    print('Verified 1,570 figure observations and the 533-pair budget-ten example.')
