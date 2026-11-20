"""Deterministic scoring in the offline, result-free scoring directory."""
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

POLICY = "SCORING_POLICY.json"
DESIGN = "CLEAN_PAIRS.tsv"
PRIOR = "TRRUST_MOUSE_SIGNED.tsv"
ENDPOINTS = ("GSE53951", "GSE114809")
MODELS = ("expression", "network_add", "network_union")
ALLOWED = {POLICY, DESIGN, PRIOR, *(s + "_EFFECT.tsv" for s in ENDPOINTS)}
READS = set()


def restrict(event, args):
    if event.startswith(("socket.", "subprocess.")):
        raise PermissionError("Score process may not use network or spawn processes")
    if event == "open":
        path, mode, _ = args
        if not isinstance(path, (str, bytes)) or not isinstance(mode, str):
            return
        path = Path(path).absolute()
        if path.parent == Path("/input") and path.name in ALLOWED and mode == "r":
            READS.add(path.name)
        elif path.parent == Path("/out") and "w" in mode:
            pass
        elif str(path).startswith("/usr/") and mode == "r":
            pass
        else:
            raise PermissionError(f"Score process attempted undeclared access: {path}")


def read(name):
    with Path(name).open(encoding="utf-8-sig", newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def write(name, rows, fields):
    with Path("/out", name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def values(path):
    effects = {}
    for row in read(path):
        symbol = row["symbol"]
        identity = row["stable_id"]
        if not symbol or not identity or row["d_g"] in {"", "NA", "NaN", "Inf", "-Inf"}:
            continue
        if symbol in effects:
            raise ValueError(f"Nonunique measured symbol {path}/{symbol}")
        effects[symbol] = identity, Fraction(row["d_g"])
    return effects


def stable_id(raw):
    return raw.removeprefix("GeneID:")


def percentiles(members, fx):
    ordered = sorted((fx[s][1], s) for s in members if s in fx)
    result = {}
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][0] == ordered[start][0]:
            end += 1
        mid = Fraction(start + end, 2 * len(ordered))
        for _, symbol in ordered[start:end]:
            result[symbol] = mid
        start = end
    return result


def make_program(fx, excluded, n):
    eligible = {s for s in fx if s not in excluded}
    up = sorted((s for s in eligible if fx[s][1] > 0), key=lambda s: (-fx[s][1], fx[s][0]))
    down = sorted((s for s in eligible if fx[s][1] < 0), key=lambda s: (fx[s][1], fx[s][0]))
    return eligible, {**dict.fromkeys(up[:n], 1), **dict.fromkeys(down[:n], -1)}


def network(path, eligible, program, by_id, panel_ids):
    use = {}
    match = {}
    if not Path(path).is_file():
        return use, match, False
    for row in read(path):
        tf = panel_ids.get(stable_id(row["tf_gene_id"]), row["tf_symbol"])
        target = by_id.get(stable_id(row["target_gene_id"]))
        sign = int(row["sign"])
        if sign not in {-1, 1}:
            raise ValueError("Non-signed network input")
        if target in eligible:
            if target in use.setdefault(tf, set()):
                raise ValueError(f"Duplicate measured edge {tf}/{target}")
            use[tf].add(target)
            if program.get(target) == sign:
                match.setdefault(tf, set()).add(target)
    return use, match, True


def score_row(pair, method, endpoint, role_p, measured, matches, net_ready, n_program, salt):
    i, j = pair["bhlh_symbol"], pair["partner_symbol"]
    aid, bid = pair["bhlh_id"], pair["partner_id"]
    p_i, p_j = role_p[0].get(i), role_p[1].get(j)
    a, b = matches.get(i, set()), matches.get(j, set())
    dx, dy = len(measured.get(i, ())), len(measured.get(j, ()))
    valid_identity = all(s and s not in {"NA", "unresolved"} for s in (i, j, aid, bid))
    value, reason = None, ""
    if not valid_identity:
        reason = "identity_unresolved"
    elif method == "expression":
        if p_i is None or p_j is None:
            reason = "missing_member_expression"
        else:
            value = (p_i + p_j) / 2
    elif not net_ready:
        reason = "network_resource_unavailable"
    elif not n_program:
        reason = "empty_target_program"
    elif not dx or not dy:
        reason = "no_usable_signed_edge_for_either_or_both_members"
    else:
        count = len(a) + len(b) if method == "network_add" else len(a | b)
        value = Fraction(count, n_program)
    tie = hashlib.sha256("|".join((salt, aid, bid)).encode()).hexdigest()
    return value, {
        **pair, "method": method, "endpoint": endpoint, "scoreable": int(value is not None),
        "reason": reason, "score": format(float(value), ".17g") if value is not None else "",
        "num": value.numerator if value is not None else "", "den": value.denominator if value is not None else "",
        "program_n": n_program, "bhlh_role_measured_n": len(role_p[0]), "partner_role_measured_n": len(role_p[1]),
        "bhlh_p": str(p_i) if p_i is not None else "", "partner_p": str(p_j) if p_j is not None else "",
        "bhlh_measured_degree": dx, "partner_measured_degree": dy, "bhlh_hits_n": len(a), "partner_hits_n": len(b),
        "intersection_n": len(a & b), "tie_key": tie, "rank": "", "rank_min": "", "rank_max": "",
    }


PAIR = ("variant", "pair_id", "bhlh_id", "partner_id", "bhlh_symbol", "partner_symbol")
FIELDS = (*PAIR, "method", "endpoint", "scoreable", "reason", "score", "num", "den", "program_n",
          "bhlh_role_measured_n", "partner_role_measured_n", "bhlh_p", "partner_p", "bhlh_measured_degree",
          "partner_measured_degree", "bhlh_hits_n", "partner_hits_n", "intersection_n", "tie_key",
          "rank", "rank_min", "rank_max")


def ordered(rows):
    available = sorted((x for x in rows if x[0] is not None),
                       key=lambda x: (-x[0], x[1]["tie_key"], x[1]["pair_id"]))
    start = 0
    while start < len(available):
        stop = start + 1
        while stop < len(available) and available[stop][0] == available[start][0]:
            stop += 1
        for index in range(start, stop):
            available[index][1].update(rank=index + 1, rank_min=start + 1, rank_max=stop)
        start = stop
    return [r for _, r in available] + sorted((r for v, r in rows if v is None), key=lambda r: r["pair_id"])


def main():
    sys.addaudithook(restrict)
    policy = json.loads(Path(POLICY).read_text(encoding="utf-8"))
    design = list(read(DESIGN))
    if not design or set(design[0]) != set(PAIR):
        raise ValueError("Unexpected clean design shape")
    unique = {(r["variant"], r["pair_id"]) for r in design}
    if len(unique) != len(design):
        raise ValueError("Duplicate public pair ID")
    variants = sorted({r["variant"] for r in design})
    separated = {v: [r for r in design if r["variant"] == v] for v in variants}
    panel_ids = {stable_id(r[side+"_id"]):r[side+"_symbol"] for r in design for side in ("bhlh","partner")
                 if r[side+"_id"] not in {"", "NA", "unresolved"}}
    results, targets, frozen_top, input_report = [], [], [], []
    for endpoint in ENDPOINTS:
        endpoint_file = endpoint + "_EFFECT.tsv"
        if not Path(endpoint_file).is_file():
            input_report.append({"endpoint": endpoint, "state": "blocked_missing_endpoint", "genes": 0, "network": False})
            continue
        fx = values(endpoint_file)
        by_id = {stable_id(identity):symbol for symbol,(identity, _) in fx.items()}
        if len(by_id) != len(fx):
            raise ValueError(f"Measured stable gene IDs are not unique: {endpoint}")
        excluded = set(panel_ids.values()) | {by_id[gene] for gene in panel_ids if gene in by_id}
        role_fx = {symbol:fx[by_id[gene]] for gene,symbol in panel_ids.items() if gene in by_id}
        eligible, program = make_program(fx, excluded, policy["program_each_direction"])
        measured, matches, usable = network(PRIOR, eligible, program, by_id, panel_ids)
        input_report.append({"endpoint": endpoint, "state": "scored", "genes": len(fx), "network": usable})
        for gene, sign in program.items():
            targets.append({"endpoint": endpoint, "stable_id": fx[gene][0], "symbol": gene, "direction": sign, "effect": str(fx[gene][1])})
        for rows in separated.values():
            role_p = tuple(percentiles({r[symbol] for r in rows if r[identity] not in {"", "NA", "unresolved"}}, role_fx)
                           for symbol, identity in (("bhlh_symbol", "bhlh_id"), ("partner_symbol", "partner_id")))
            for method in MODELS:
                partial = ordered([score_row(pair, method, endpoint, role_p, measured, matches, usable,
                                             len(program), policy["tie_salt"]) for pair in rows])
                results.extend(partial)
                for row in partial[:20]:
                    if row["scoreable"]:
                        frozen_top.append(row)
    write("all_pair_scores.tsv", results, FIELDS)
    write("top_blind.tsv", frozen_top, FIELDS)
    write("program.tsv", targets, ("endpoint", "stable_id", "symbol", "direction", "effect"))
    with Path("/out/READ_ACCESS.json").open("w", encoding="utf-8") as stream:
        json.dump({"opened_files": sorted(READS), "endpoints": input_report}, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
