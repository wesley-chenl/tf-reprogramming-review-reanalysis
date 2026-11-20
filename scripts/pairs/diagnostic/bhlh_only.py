"""Recompute the matched 65-pair bHLH ablations without refitting or new data.
Run with Python 3 and NumPy. The score function never takes screen labels.
"""
from pathlib import Path
from collections import defaultdict
from fractions import Fraction
import csv
import gzip
import numpy as np

ROOT = Path(__file__).resolve().parent

def read(root, name):
    with gzip.open(root / name, 'rt', encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(out, name, rows):
    with (out / name).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

def rank_blind(inputs=ROOT):
    panel = read(inputs, 'panel_inputs.tsv.gz')
    effects = read(inputs, 'gene_effect.tsv.gz')
    edges = read(inputs, 'direct_edges.tsv.gz')
    index = {s['stable_id']: j for j, s in enumerate(effects)}
    effect = np.array([float(s['delta']) for s in effects])
    columns = defaultdict(lambda: np.zeros(len(effect)))
    for e in edges:
        columns[e['tf_gene_id']][index[e['target_gene_id']]] = int(e['sign'])
    by_method = defaultdict(list)
    for r in panel:
        b = columns[r['bhlh_id']] * float(r['bhlh_coef'])
        p = columns[r['partner_id']] * float(r['partner_coef'])
        complete = np.corrcoef(effect, b+p)[0,1]
        assert round(complete,12) == round(float(r['tfcomb_pair_score']),12)
        scores = {'Expression: both': Fraction(r['expression_pair_fraction']),
                  'Expression: bHLH only': Fraction(r['bhlh_fraction']),
                  'TFcomb: both': round(complete, 12),
                  'TFcomb: bHLH only': round(np.corrcoef(effect,b)[0,1],12)}
        for name,score in scores.items():
            by_method[name].append(dict(slot=r['slot'], bhlh=r['bhlh_symbol'], partner=r['partner_symbol'],
                                        score=score, tie_key=r['tie_key']))
    result = []
    for method, data in by_method.items():
        rows = sorted(data, key=lambda a: (-a['score'], a['tie_key']))
        for i,r in enumerate(rows,1):
            others=[j for j,a in enumerate(rows,1) if a['score']==r['score']]
            result.append(dict(method=method, slot=r['slot'], bhlh=r['bhlh'],partner=r['partner'],
                               score=r['score'],rank=i,midrank=sum(others)/len(others)))
    return result

def evaluate(ranks, labels=ROOT):
    annotations = {r['slot']:r for r in read(labels, 'screen_outcomes.tsv.gz')}
    summary = []
    for name in dict.fromkeys(r['method'] for r in ranks):
        rows=[dict(r, positive=int(annotations[r['slot']]['positive']),
                   gated=float(annotations[r['slot']]['gated_yield_pct'])) for r in ranks if r['method']==name]
        assert len(rows)==65 and sum(r['positive'] for r in rows)==18
        match=[(p,n) for p in rows if p['positive'] for n in rows if not n['positive']]
        within=[(p,n) for p,n in match if p['bhlh']==n['bhlh']]
        def auc(pairs):
            return sum((p['midrank']<n['midrank']) + .5*(p['midrank']==n['midrank']) for p,n in pairs)/len(pairs)
        for budget in (10,20):
            picks=rows[:budget]
            def within_quota(r):
                group=[v for v in rows if v['bhlh']==r['bhlh']]
                return sum(v['positive'] for v in group)/len(group)
            summary.append(dict(method=name,budget=budget,N=len(rows),positive=sum(r['positive'] for r in rows),
                hits=sum(r['positive'] for r in picks),random_hits=budget*18/65,
                quota_expected_hits=sum(map(within_quota,picks)),
                negative_in_bottom_budget=sum(1-r['positive'] for r in rows[-budget:]),
                AUROC=auc(match),within_bhlh_concordance=auc(within),within_comparisons=len(within),
                mean_gated_yield_pct=np.mean([r['gated'] for r in picks]),
                best_regret_pp=max(r['gated'] for r in rows)-max(r['gated'] for r in picks)))
    return summary

if __name__=='__main__':
    from argparse import ArgumentParser
    from hashlib import sha256
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--labels', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    ranks = rank_blind(args.inputs)
    rank_name = 'tfcomb65_blind_ranks.tsv'
    write(args.out, rank_name, ranks)
    (args.out / 'TFComb_DIAGNOSTIC_SCORE.sha256').write_text(
        sha256((args.out / rank_name).read_bytes()).hexdigest() + '  ' + rank_name + '\n',
        encoding='ascii',
    )
    write(args.out, 'tfcomb65_diagnostic.tsv', evaluate(ranks, args.labels))
