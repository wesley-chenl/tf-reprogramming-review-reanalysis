"""Recompute §5.2 tables from frozen ranks and evaluate the two manuscript ablations."""
from argparse import ArgumentParser
from collections import Counter, defaultdict
from fractions import Fraction
import json
from math import isclose
from pathlib import Path
from statistics import mean
from types import SimpleNamespace

import rebuild_rank_panels as panels
from diagnostic import bhlh_only, expression_bhlh_210


DATA = Path(__file__).resolve().parents[2] / "data/pairs"
KEY = ("endpoint", "variant", "panel_type", "method", "utility", "B")
MEMBERS = ("bhlh_id", "partner_id")
EP = ("GSE53951", "GSE114809")


def equal(result, expected):
    for k, v in result.items():
        target = expected[k]
        if isinstance(v, (int, float, Fraction)):
            if not isclose(float(v), float(target), rel_tol=1e-10, abs_tol=1e-12):
                raise ValueError((k, v, target))
        elif str(v) != str(target):
            raise ValueError((k, v, target))


def frozen_utility(data, out):
    # Use the positional identity and outcome criteria of the original evaluator.
    design = {(r["variant"], r["pair_id"]): r for r in panels.rows(data / "labels/design_pairs.tsv")}
    ids = Counter((v, *(r[x] for x in MEMBERS)) for (v, _), r in design.items())
    outcomes = {r["pair_slot_id"]: r for r in panels.rows(data / "labels/outcomes_private.tsv")
                if r["source_version"] == "nature2018"}
    scores: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for r in panels.rows(data / "utility/scores/all_pair_scores.tsv"):
        ep, v, m = r["endpoint"], r["variant"], r["method"]
        z = design[v, r["pair_id"]]
        assert all(r[x] == z[x] for x in MEMBERS)
        if r["scoreable"] != "1" or z["eligible"] not in {"1", "True", "true"}:
            continue
        if ids[v, *(r[x] for x in MEMBERS)] != 1:
            continue
        slot = z["pair_slot_id"]
        readout = outcomes[slot]
        for u, value_key in (("binary", "binary_call"), ("gated", "q_gated")):
            if readout[value_key] not in {"", "NA", "unknown", "None"}:
                scores[ep, v, m, u].append(dict(
                    pair_id=r["pair_id"], slot=slot, bhlh=r["bhlh_symbol"], partner=r["partner_symbol"],
                    bhlh_id=r["bhlh_id"], global_score_rank=int(r["rank"]), score=r["score"],
                    binary=int(readout["binary_call"]), q=Fraction(readout[value_key]),
                    gated=Fraction(readout["q_gated"]), raw_s1_pct=readout["raw_s1_pct"],
                    source_locator=readout["source_cell_or_bbox"], measurement_status=readout["measurement_status"],
                ))
    common = defaultdict(list)
    for (ep, v, m, u), ranked in scores.items():
        ranked.sort(key=lambda r: r["global_score_rank"])
        common[ep, v, u].append({r["pair_id"] for r in ranked})
    null = {tuple(r[x] for x in KEY): r for r in panels.rows(data / "utility/evaluation/partner_null.tsv")}
    metrics, nominations, checks = [], [], []
    for archived in panels.rows(data / "utility/evaluation/utility_results.tsv"):
        ep, v, scope, method, utility, b = [archived[k] for k in KEY]
        budget = int(b)
        native = scores[ep, v, method, utility]
        match = set.intersection(*common[ep, v, utility])
        pool = native if scope == "native" else [r for r in native if r["pair_id"] in match]
        ids = {r["pair_id"] for r in pool}
        assert len(ids) == len(pool) == int(archived["N"])
        if len(pool) < budget or archived["mean_utility"] == "":
            continue
        top = pool[:budget]
        by_bhlh = defaultdict(list)
        for z in pool:
            by_bhlh[z["bhlh_id"]].append(z["q"])
        exact_null_mean = mean(r["q"] for r in pool)
        fixed_bhlh_mean = mean(mean(by_bhlh[r["bhlh_id"]]) for r in top)
        quantities = dict(
            N=len(pool), mean_utility=mean(r["q"] for r in top), null_mean=exact_null_mean,
            best_selected=max(r["q"] for r in top), oracle=max(r["q"] for r in pool),
            best_regret=max(r["q"] for r in pool) - max(r["q"] for r in top),
        )
        if utility == "binary":
            quantities.update(K=sum(r["q"] for r in pool), hit_n=sum(r["q"] for r in top))
        equal(quantities, archived)
        equal(dict(partner_null_mean=fixed_bhlh_mean), null[ep, v, scope, method, utility, b])
        k = {field: archived[field] for field in KEY}
        metrics.append(dict(**k, **{x: float(y) if isinstance(y, Fraction) else y for x, y in quantities.items()},
                            fixed_bhlh_mean=float(fixed_bhlh_mean),
                            frozen_positive=sum(r["binary"] for r in top),
                            uniform_expected_hits=float(Fraction(budget * sum(r["binary"] for r in pool), len(pool)))))
        for rank, row in enumerate(top, 1):
            nominations.append(dict(**k, eval_rank=rank, **{x: float(y) if isinstance(y, Fraction) else y
                                                           for x, y in row.items()}))
        if (ep, v, scope, method, utility, b) == ("GSE53951", "main", "native", "expression", "gated", "10"):
            checks = [dict(eval_rank=i, **r) for i, r in enumerate(pool[:12], 1)]
            n, k0 = len(pool), sum(r["binary"] for r in pool)
            assert (n, k0, sum(r["binary"] for r in top)) == (533, 69, 3)
            rank11 = checks[10]["gated"]
            assert quantities["mean_utility"] == Fraction("0.011") and rank11 == Fraction("1.60")
            figure4 = [dict(eligible=n, positives=k0, selected_positive=3, random_expected_hits=10*k0/n,
                            mean_yield_pct=float(quantities["mean_utility"]), best_selected_pct=float(quantities["best_selected"]),
                            uniform_mean_yield_pct=float(exact_null_mean), rank11_yield_pct=float(rank11),
                            rank11_minus_mean_pp=float(rank11-quantities["mean_utility"]),
                            rank11_minus_best_pp=float(rank11-quantities["best_selected"]),
                            panel_max_pct=float(quantities["oracle"]), regret_pp=float(quantities["best_regret"]))]
            equal(figure4[0], next(panels.rows(data / "manuscript/budget10_figure_summary.tsv")))
            panels.write(out / "FIG4_RECOMPUTED.tsv", figure4)
    assert len(checks) == 12
    panels.write(out / "FIXED_BUDGET_UTILITY.tsv", metrics)
    panels.write(out / "FIXED_BUDGET_NOMINATIONS.tsv", nominations)
    panels.write(out / "FIG4_RANKED_FIRST12.tsv", [dict(r, q=float(r["q"]), gated=float(r["gated"])) for r in checks])


def matched_paper_metrics(data, out, old, followups, full):
    summary = {(r["panel"], r["method"]): r for r in panels.rows(out / "PANEL_METHOD_SUMMARY.tsv")}
    for archived in panels.rows(data / "manuscript/mirrored_panel_metrics.tsv"):
        p, m = archived["panel"], archived["method"]
        s = summary[p, m]
        assert (int(s["N"]), int(s["K"])) == {"A": (210, 37), "B": (90, 7), "C": (65, 18)}[p]
        assert len(full[p, m]) == int(s["N"])
        equal(dict(N=int(s["N"]), K=int(s["K"]), neg_median=float(s["negative_median_midrank"]),
                   pos_median=float(s["positive_median_midrank"]), auroc=float(s["AUROC"]),
                   within_bhlh=float(s["bhlh_weighted_concordance"]), within_contrasts=int(s["bhlh_positive_negative_comparisons"]),
                   hits10=int(s["top10_positive_frozen_count"]), hits20=int(s["top20_positive_frozen_count"]),
                   negative_bottom10=int(s["bottom10_negative_frozen_count"]),
                   quota10=float(s["top10_positive_fixed_bhlh_expected_count"]),
                   quota20=float(s["top20_positive_fixed_bhlh_expected_count"])), archived)
    panel_b = full["B", "taiji_pair"].keys()
    endpoint_rows = []
    for ep in EP:
        rec = panels.rerank(old[ep, "expression"], panel_b)
        stats, *_ = panels.metric(dict(panel="B", endpoint=ep, method="expression"), rec, detail=False)
        expected = {"GSE53951": 0.682, "GSE114809": 0.654}[ep]
        assert round(stats["AUROC"], 3) == expected
        endpoint_rows.append(stats)
    panels.write(out / "TAIJI_FIXED_SLOTS_EXPRESSION_ENDPOINTS.tsv", endpoint_rows)

    a_slots = full["A", "expression"].keys()
    expected_auc = {
        "GSE53951": (0.706, 0.673, 0.571, 0.613, 0.637, 0.734),
        "GSE114809": (0.662, 0.678, 0.635, 0.650, 0.645, 0.758),
    }
    endpoint_sensitivity = []
    for ep in EP:
        for method, auc in zip(panels.A, expected_auc[ep]):
            group = old[ep, method]
            match = panels.valid(group, "reference")
            rec = panels.rerank(group, match)
            stats, *_ = panels.metric(dict(endpoint=ep, method=method, scope="endpoint_matched"), rec, detail=False)
            assert (stats["N"], stats["K"]) == ((210, 37) if ep == EP[0] else (220, 37))
            assert round(stats["AUROC"], 3) == auc
            endpoint_sensitivity.append(stats)
            if ep == EP[1] and method in {"two_union", "imcell_p010"}:
                rec = panels.rerank(group, a_slots)
                stats, *_ = panels.metric(dict(endpoint=ep, method=method, scope="GSE53951_matched_slots"), rec, detail=False)
                assert round(stats["AUROC"], 3) == {"two_union": 0.624, "imcell_p010": 0.641}[method]
                endpoint_sensitivity.append(stats)
    panels.write(out / "ENDPOINT_SENSITIVITY.tsv", endpoint_sensitivity)

    pair_groups = []
    for ep in EP:
        for scheme in ("p010", "in_degree", "mean"):
            for mode in ("joint", "additive"):
                method = f"imcell_weight_{scheme}_{mode}"
                group = followups[ep, method]
                match = panels.valid(old[ep, "expression"], "reference")
                rec = panels.rerank(group, match)
                assert len(rec) == len(match)
                top = rec[:10]
                pair_groups.append(dict(endpoint=ep, scheme=scheme, mode=mode, N=len(rec),
                                        K=sum(r["screen_positive"] for r in rec),
                                        hits10=sum(r["screen_positive"] for r in top),
                                        Myc_pairs10=sum(r["bhlh_symbol"] == "Myc" for r in top),
                                        mean_gated_yield_pct=float(mean(Fraction(r["q_gated"]) for r in top))))
    panels.write(out / "IMCELL_EDGE_PROBABILITY_SENSITIVITY.tsv", pair_groups)

    priors = []
    for ep in EP:
        for method in ("network_union", "two_union"):
            base = old[ep, method]
            for label in ("modern", "through2015"):
                alt = followups[ep, f"collectri_{label}_{method}"]
                a, b = panels.valid(base, "native"), panels.valid(alt, "native")
                common = a & b
                before, _ = panels.pair_auc(panels.rerank(base, common))
                after, _ = panels.pair_auc(panels.rerank(alt, common))
                priors.append(dict(endpoint=ep, method=method, prior=label, TRRUST_N=len(a), CollecTRI_N=len(b),
                                   shared_N=len(common), TRRUST_AUROC=float(before), CollecTRI_AUROC=float(after)))
    panels.write(out / "SIGNED_PRIOR_SENSITIVITY.tsv", priors)


def tfcomb_flow(data, out, old, full):
    ep = "GSE53951"
    adapter = {r["pair_id"]: r for r in panels.rows(data / "initial/tfcomb_primary_adapter.tsv")}
    pool = old[ep, "expression"]
    scorer_slots = full["C", "tfcomb"].keys()
    detail = []
    for slot, r in sorted(pool.items()):
        record = adapter[r["pair_id"]]
        expressed = bool(r["native_rank"])
        fitted = json.loads(record["details"])
        coefficients = [fitted["bhlh_coef"], fitted["partner_coef"]]
        connects = expressed and all(v is not None for v in coefficients)
        positive = connects and all(v > 0 for v in coefficients)
        assert positive == (slot in scorer_slots)
        detail.append(dict(slot=slot, pair_id=r["pair_id"], screen_positive=r["screen_positive"],
                           expression_pass=int(expressed), graph_pass=int(connects), both_positive=int(positive),
                           bhlh_coef=coefficients[0], partner_coef=coefficients[1], tfcomb_reason=record["reason"]))
    result = [dict(stage=key, pairs=sum(r[key] for r in detail),
                   positives=sum(r[key]*r["screen_positive"] for r in detail))
              for key in ("expression_pass", "graph_pass", "both_positive")]
    assert [(r["pairs"], r["positives"]) for r in result] == [(533, 69), (230, 45), (65, 18)]
    panels.write(out / "TFCOMB_ELIGIBILITY.tsv", result)
    panels.write(out / "TFCOMB_ELIGIBILITY_SLOTS.tsv", detail)


def diagnostics(data, out, full):
    # Artifact generation below reads separately stored blind inputs for each ablation.
    hold = data / "manuscript/diagnostic"
    path = out / "diagnostic"
    path.mkdir(exist_ok=True)
    blind = bhlh_only.rank_blind(data / "diagnostic/score_inputs")
    name = "tfcomb65_blind_ranks.tsv"
    bhlh_only.write(path, name, blind)
    assert (path / name).read_bytes() == (hold / name).read_bytes()
    (path / "TFComb_DIAGNOSTIC_SCORE.sha256").write_text(panels.digest(path / name) + "  " + name + "\n", encoding="ascii")
    results = bhlh_only.evaluate(blind, data / "diagnostic/labels")
    name = "tfcomb65_diagnostic.tsv"
    for actual, archive in zip(results, panels.rows(hold / name), strict=True):
        equal(actual, archive)
    bhlh_only.write(path, name, results)
    members = path / "expression_shared_210_slot_membership.tsv"
    panels.write(members, [dict(slot=r["slot"], pair_id=r["pair_id"]) for r in full["A", "expression"].values()])
    ranked = expression_bhlh_210.score_pairs(members, data / "utility/scores/all_pair_scores.tsv")
    name = "expression_bhlh_210_ranks.tsv"
    expression_bhlh_210.write_table(path / name, ranked)
    assert (path / name).read_bytes() == (hold / name).read_bytes()
    (path / "EXPRESSION_DIAGNOSTIC_SCORE.sha256").write_text(panels.digest(path / name) + "  " + name + "\n", encoding="ascii")
    computed = expression_bhlh_210.summarize(ranked, out / "PLOT_SLOTS.tsv")
    name = "expression_bhlh_210_summary.tsv"
    equal(computed, next(panels.rows(hold / name)))
    expression_bhlh_210.write_table(path / name, [computed])


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    frozen_utility(DATA, out)
    old, followups, full = panels.main(SimpleNamespace(
        old=DATA / "ranks/pair_ranks_initial.tsv.gz", new=DATA / "ranks/pair_ranks_followup.tsv.gz", out=out))
    matched_paper_metrics(DATA, out, old, followups, full)
    tfcomb_flow(DATA, out, old, full)
    diagnostics(DATA, out, full)
    print("Pair results verified: Figure 4, Figure 5, fixed-budget tables, TFcomb eligibility and manuscript diagnostics.")


if __name__ == "__main__":
    main()
