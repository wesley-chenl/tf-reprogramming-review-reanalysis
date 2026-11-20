"""Reuse the frozen direct and two-hop scoring definitions with a new graph."""

from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import sys

from historical_score import read, network, make_program, score_row, stable_id, values, write
from core_score import path_members

ENDPOINTS = ("GSE53951", "GSE114809")
DIRECT = ("network_add", "network_union")
TWO_HOP = ("two_add", "two_union")
SIDES = ("bhlh_id", "partner_id")


def details_for(row):
    keys = ("bhlh_measured_degree", "partner_measured_degree", "bhlh_hits_n", "partner_hits_n", "intersection_n")
    return {k: row[k] for k in keys}


def output_row(pair, endpoint, prior, base, val, reason, details):
    return dict(endpoint=endpoint, variant=pair["variant"], pair_id=pair["pair_id"], method=f"{prior}_{base}",
                score=format(float(val), ".17g") if val is not None else "", eligible=int(val is not None),
                exclusion_reason=reason, score_num=val.numerator if val is not None else "",
                score_den=val.denominator if val is not None else "", prior=prior, base_method=base,
                details=json.dumps(details, sort_keys=True, separators=(",", ":")))


def path_check(edges, tf, allowed):
    """Independent explicit join of first edges to a lookup by middle node."""
    outgoing = defaultdict(list)
    for a, b, s in edges:
        outgoing[a].append((b, s))
    events = [(b, s) for b, s in outgoing[tf] if b != tf]
    events += [(target, x * y) for m, x in events for target, y in outgoing[m] if target not in {m, tf}]
    pos = {g for g, s in events if s == 1} & allowed
    neg = {g for g, s in events if s == -1} & allowed
    return pos - neg, neg - pos


def independent_checks(design, edges, allowed, program, scores, ep):
    lookup = {(r["variant"], r["pair_id"], r["base_method"]): r for r in scores if r["endpoint"] == ep}
    by_id = {stable_id(g[0]): s for s, g in allowed.items()}
    direction = {g: program[s] for g, s in by_id.items() if s in program}
    direct = defaultdict(dict)
    for a, b, sign in edges:
        if b in by_id:
            direct[a][b] = sign
    candidates = [r for r in design if r["variant"] == "main"]
    candidates.sort(key=lambda r: sha256(("prior_coverage_replay|" + r["pair_id"]).encode()).hexdigest())
    outputs = []
    for pair in candidates[:24]:
        a, b = (stable_id(pair[side]) for side in SIDES)
        valid = a.isdigit() and b.isdigit()
        ad = {g for g, s in direct[a].items() if direction.get(g) == s}
        bd = {g for g, s in direct[b].items() if direction.get(g) == s}
        distinct_direct = bool(direct[a] and direct[b])
        left, right = (path_check(edges, t, by_id.keys()) for t in (a, b))
        ac = {g for g in left[0] if direction.get(g) == 1} | {g for g in left[1] if direction.get(g) == -1}
        bc = {g for g in right[0] if direction.get(g) == 1} | {g for g in right[1] if direction.get(g) == -1}
        distinct_hop = a != b and any(left) and any(right)
        for base, hits, eligible in (
            ("network_add", len(ad) + len(bd), distinct_direct),
            ("network_union", len(ad | bd), distinct_direct),
            ("two_add", len(ac) + len(bc), distinct_hop),
            ("two_union", len(ac | bc), distinct_hop),
        ):
            derived = Fraction(hits, len(program)) if valid and eligible else None
            r = lookup[pair["variant"], pair["pair_id"], base]
            saved = Fraction(r["score_num"], r["score_den"]) if r["eligible"] else None
            if saved != derived:
                raise AssertionError((r, saved, derived))
            outputs.append(dict(endpoint=ep, pair_id=pair["pair_id"], method=base, score=str(derived) if derived is not None else "",
                                check="same"))
    return outputs


def main(prior):
    design = list(read("CLEAN_PAIRS.tsv"))
    policy = json.loads(Path("SCORING_POLICY.json").read_text(encoding="utf-8"))
    program_rows = list(read("PROGRAM.tsv"))
    assert len(design) == 598 * 9 and len({(r["variant"], r["pair_id"]) for r in design}) == len(design)
    panel_ids = {stable_id(r[key]): r[sym] for r in design
                 for key, sym in zip(SIDES, ("bhlh_symbol", "partner_symbol")) if r[key]}
    panel = set(panel_ids)
    edges = [(r["tf_gene_id"], r["target_gene_id"], int(r["sign"])) for r in read("GRAPH.tsv")]
    assert len(edges) == len({(a, b) for a, b, s in edges}) and all(s in (-1, 1) for a, b, s in edges)
    outputs, qa = [], []
    for ep in ENDPOINTS:
        fx = values(ep + "_EFFECT.tsv")
        symbols_by_id = {stable_id(g): k for k, (g, _) in fx.items()}
        excluded = set(panel_ids.values()) | {symbols_by_id[g] for g in panel if g in symbols_by_id}
        allowed, programme = make_program(fx, excluded, policy["program_each_direction"])
        expected = {r["symbol"]: int(r["direction"]) for r in program_rows if r["endpoint"] == ep}
        assert programme == expected and len(programme) == 500
        measured, matches, ready = network(Path("GRAPH.tsv"), allowed, programme, symbols_by_id, panel_ids)
        by_gene = path_members(edges, panel, symbols_by_id, {fx[s][0]: d for s, d in programme.items()})
        ep_rows = []
        for pair in design:
            for name in DIRECT:
                val, row = score_row(pair, name, ep, ({}, {}), measured, matches, ready, len(programme), policy["tie_salt"])
                ep_rows.append(output_row(pair, ep, prior, name, val, row["reason"], details_for(row)))
            a, b = (stable_id(pair[key]) for key in SIDES)
            if a in by_gene and b in by_gene and a != b:
                x, y = by_gene[a], by_gene[b]
                p, q = x["hits"], y["hits"]
                reason = "" if x["degree"] and y["degree"] else "member_has_no_unambiguous_measured_reachable_target"
                details = dict(member_degree=[x["degree"], y["degree"]], member_hits=[sorted(p), sorted(q)],
                               intersection_ids=sorted(p & q), conflicting_measured=[x["conflicting"], y["conflicting"]],
                               program_n=len(programme))
                for name, n in (("two_add", len(p) + len(q)), ("two_union", len(p | q))):
                    ep_rows.append(output_row(pair, ep, prior, name, None if reason else Fraction(n, len(programme)), reason, details))
            else:
                for name in TWO_HOP:
                    ep_rows.append(output_row(pair, ep, prior, name, None, "unresolved_or_equal_gene_id", {}))
        qa.extend(independent_checks(design, edges, {k: fx[k] for k in allowed}, programme, ep_rows, ep))
        outputs.extend(ep_rows)
    write("COVERAGE_SCORES.tsv", outputs, list(outputs[0]))
    write("COVERAGE_INDEPENDENT.tsv", qa, list(qa[0]))
    print(f"{prior}: {len(edges)} edges, {len(outputs)} physical pair score records, {len(qa)} independent rule comparisons")


if __name__ == "__main__":
    main(sys.argv[1])
