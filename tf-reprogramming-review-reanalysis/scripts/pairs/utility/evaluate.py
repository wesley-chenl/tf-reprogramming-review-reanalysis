"""Join locked screens and scores. Compute same-panel and conditional nulls."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from math import comb, sqrt
import os
from pathlib import Path

import numpy as np
from common import check_lock, load, now, read, record, sha, table

MODELS = ("expression", "network_add", "network_union")
BUDGETS = (10, 20)
SEED = 20261119
UNIFORM_DRAWS = PARTNER_DRAWS = 100000
TIE_DRAWS = 10000
UNKNOWN = {"", "NA", "unknown", "None"}
TRUE = {"1", "true", "True"}


def number(row, key):
    value = row.get(key, "")
    if value in UNKNOWN:
        return None
    if value.endswith("%"):
        value = value[:-1]
    return float(value)


def rng(key):
    digest = sha256(repr(key).encode()).digest()
    return np.random.default_rng(np.random.SeedSequence([SEED, int.from_bytes(digest[:8], "big")]))


def sample(n, size, times, generator):
    """Sample each row uniformly without replacement. All quasi-random work uses a private seed."""
    if size == n:
        return np.broadcast_to(np.arange(n), (times, n))
    if n < size * 3:
        keys = generator.random((times, n))
        return np.argpartition(keys, size - 1, axis=1)[:, :size]
    result = np.empty((times, size), dtype=np.int32)
    for j in range(size):
        ballot = generator.integers(0, n, size=times, dtype=np.int32)
        collision = np.any(ballot[:, None] == result[:, :j], axis=1)
        while collision.any():
            ballot[collision] = generator.integers(0, n, size=collision.sum(), dtype=np.int32)
            collision = np.any(ballot[:, None] == result[:, :j], axis=1)
        result[:, j] = ballot
    return result


def quantile(x):
    return np.quantile(x, [0.025, 0.975], method="inverted_cdf").tolist()


def max_distribution(values, size):
    ordered = sorted(values)
    denom = comb(len(ordered), size)
    masses, last = [], 0
    for index in range(size - 1, len(ordered)):
        stops = comb(index + 1, size)
        masses.append((ordered[index], (stops - last) / denom, stops / denom))
        last = stops
    return masses


def uniform(values, size, key, binary):
    n = len(values)
    array = np.asarray(values)
    average = float(array.mean())
    oracle = float(array.max())
    masses = max_distribution(values, size)
    regret_mean = sum((oracle - val) * mass for val, mass, _ in masses)
    def regret_quantile(p):
        running = 0.0
        for v, weight, _ in reversed(masses):
            running += weight
            if running >= p:
                return oracle - v
        return oracle - masses[0][0]
    sampled = array[sample(n, size, UNIFORM_DRAWS, rng(("uniform", key, size)))]
    mean_s = sampled.mean(axis=1)
    out = {
        "null_mean": average, "randomization_mean_q025": quantile(mean_s)[0], "randomization_mean_q975": quantile(mean_s)[1],
        "null_regret_mean": regret_mean, "randomization_regret_q025": regret_quantile(0.025),
        "randomization_regret_q975": regret_quantile(0.975), "random_zero_regret_probability": sum(m for v, m, _ in masses if v == oracle),
        "uniform_mc_mean": float(mean_s.mean()), "uniform_mc_zero_regret_probability": float((sampled.max(axis=1) == oracle).mean()),
        "top_b_oracle_mean": float(np.mean(sorted(values)[-size:])), "uniform_draws": UNIFORM_DRAWS,
    }
    if binary:
        k = int(array.sum())
        weights = [((comb(k, h) * comb(n-k, size-h)) / comb(n, size), h)
                   for h in range(max(0, size - (n-k)), min(size, k) + 1)]
        def hyperquant(p):
            total = 0.0
            for mass, h in weights:
                total += mass
                if total >= p:
                    return h / size
            return weights[-1][1] / size
        zero_prob = 1.0 if k == 0 or n - k < size else 1 - comb(n-k, size) / comb(n, size)
        tol = 6 * sqrt(zero_prob * (1-zero_prob) / UNIFORM_DRAWS) + 1 / UNIFORM_DRAWS
        if abs(out["uniform_mc_zero_regret_probability"] - zero_prob) > tol or abs(out["random_zero_regret_probability"] - zero_prob) > 1e-12:
            raise ValueError("Binary analytic and sampled zero-regret checks disagree")
        out.update(randomization_mean_q025=hyperquant(0.025), randomization_mean_q975=hyperquant(0.975),
                   random_zero_regret_probability=zero_prob, null_distribution="exact_hypergeometric")
    else:
        out["null_distribution"] = "uniform_100000_without_replacement_mean_exact_order_statistic_regret"
    return out


def partner_null(pool, selected, utility, size, key):
    quotas = Counter(row["bhlh_id"] for row in selected)
    groups = defaultdict(list)
    for row in pool:
        if row["bhlh_id"] in quotas:
            groups[row["bhlh_id"]].append(row[utility])
    mean_total, oracle, sampled_rows, variable = 0.0, max(r[utility] for r in pool), [], 0
    generator = rng(("partner", key, size, sorted(quotas.items())))
    for role, b in sorted(quotas.items()):
        values = np.asarray(groups[role])
        mean_total += b * float(values.mean())
        variable += int(len(values) > b)
        sampled_rows.append(values[sample(len(values), b, PARTNER_DRAWS, generator)])
    mean = sum(a.sum(axis=1) for a in sampled_rows) / size
    regret = oracle - np.maximum.reduce([a.max(axis=1) for a in sampled_rows])
    bounds, reg_bounds = quantile(mean), quantile(regret)
    observed = sum(row[utility] for row in selected) / size
    return {
        "observed_mean": observed, "partner_null_mean": mean_total / size, "delta_over_fixed_bhlh": observed - mean_total / size,
        "partner_randomization_mean_q025": bounds[0], "partner_randomization_mean_q975": bounds[1],
        "partner_randomization_regret_q025": reg_bounds[0], "partner_randomization_regret_q975": reg_bounds[1],
        "partner_mc_mean": float(mean.mean()), "partner_draws": PARTNER_DRAWS,
        "variable_bhlh_rows": variable, "nominated_bhlh_rows": len(quotas),
        "bhlh_quotas": ",".join(f"{role}:{size}" for role, size in sorted(quotas.items())),
        "bhlh_pool_sizes": ",".join(f"{role}:{len(values)}" for role, values in sorted(groups.items())),
    }


def boundary(pool, size, utility, oracle, key):
    score = pool[size-1]["num"], pool[size-1]["den"]
    block = [i for i, r in enumerate(pool) if (r["num"], r["den"]) == score]
    choice = size - block[0]
    summary = {"boundary_tie_n": len(block), "boundary_places": choice, "boundary_fraction": choice / len(block),
               "score_unique_n": len({(r["num"], r["den"]) for r in pool}), "tie_draws": 0}
    if choice < len(block):
        fixed = [row[utility] for row in pool[:block[0]]]
        in_block = np.array([pool[j][utility] for j in block])
        draws = in_block[sample(len(in_block), choice, TIE_DRAWS, rng(("tie", key, size)))]
        mean = (sum(fixed) + draws.sum(axis=1)) / size
        regret = oracle - np.maximum(max(fixed, default=0), draws.max(axis=1))
        qmean, qregret = quantile(mean), quantile(regret)
        summary.update(tie_draws=TIE_DRAWS, tie_mean_q025=qmean[0], tie_mean_q975=qmean[1],
                       tie_regret_q025=qregret[0], tie_regret_q975=qregret[1], tie_mean_min=float(mean.min()), tie_mean_max=float(mean.max()))
    return summary


def evaluate_one(task):
    key, method, endpoint, variant, panel_kind, utility, pool, size, global_oracle, null = task
    panel_n = len(pool)
    values = [r[utility] for r in pool]
    status = {"endpoint": endpoint, "variant": variant, "panel_type": panel_kind,
              "method": method, "utility": utility, "B": size, "N": panel_n,
              "K": sum(bool(row["binary"]) for row in pool),
              "status": "insufficient_n" if panel_n < size else "no_choice" if panel_n == size else "evaluated"}
    if panel_n < size:
        return status, None, []
    top = pool[:size]
    oracle = max(values)
    observed = sum(r[utility] for r in top) / size
    best = max(r[utility] for r in top)
    details = boundary(pool, size, utility, oracle, (key, method))
    status.update(
        hit_n=sum(bool(r["binary"]) for r in top), mean_utility=observed,
        best_selected=best, oracle=oracle, best_regret=oracle-best, full_position_oracle=global_oracle,
        full_position_regret=global_oracle-best if global_oracle is not None else None,
        delta_over_uniform=observed-null["null_mean"],
        raw_s1_selected_coverage=sum(r["raw_s1"] is not None for r in top),
        all_equal=int(details["score_unique_n"] == 1),
        utility_uninformative=int(len(set(values)) <= 1),
        **details, **null
    )
    partner = {k: status[k] for k in ("endpoint", "variant", "panel_type", "method", "utility", "B", "N", "K")}
    partner.update(partner_null(pool, top, utility, size, key))
    rows = [{
        **partner_key(status), "eval_rank": i, "global_score_rank": row["rank"], "pair_id": row["pair_id"],
        "pair_slot_id": row["pair_slot_id"], "bhlh": row["bhlh_symbol"], "partner": row["partner_symbol"],
        "score": row["score"], "exact_num": row["num"], "exact_den": row["den"], "q": row[utility],
        "binary_call": row["binary"], "q_gated": row["gated"], "raw_s1_pct": row["raw_s1"],
        "source_locator": row["locator"], "call_basis": row["call_basis"],
        "measurement_status": row["measurement_status"], "design_exclusion_reason": row["design_exclusion_reason"],
    } for i, row in enumerate(top, 1)]
    return status, partner, rows


def partner_key(row):
    return {k: row[k] for k in ("endpoint", "variant", "panel_type", "method", "utility", "B")}


def main(root, workers):
    scoring = check_lock(root, "SCORE")
    labels = check_lock(root, "LABEL")
    output = root / "evaluation"
    output.mkdir(parents=True, exist_ok=True)
    join = output / "JOIN_STARTED.json"
    digests = {"score_lock_sha256": sha(root / "locks/SCORE_LOCK.json"),
               "label_lock_sha256": sha(root / "locks/LABEL_LOCK.json")}
    if not join.is_file():
        record(join, {"at": now(), "score_sealed": scoring["sealed_utc"], "labels_sealed": labels["sealed_utc"], **digests})
    elif any(load(join)[k] != val for k,val in digests.items()):
        raise ValueError("Cannot change locks after first join")
    private = read(root / "sources_private/design_pairs.tsv")
    design = {(r["variant"], r["pair_id"]): r for r in private}
    if len(design) != len(private):
        raise ValueError("Duplicate private design join key")
    records = [r for r in read(root / "sources_private/outcomes_private.tsv") if r["source_version"] == "nature2018"]
    outcomes = {r["pair_slot_id"]:r for r in records}
    if len(outcomes) != len(records):
        raise ValueError("Duplicate Nature position in outcomes")
    variants = sorted({v for v, _ in design})
    identities = Counter((r["variant"], r["bhlh_id"], r["partner_id"]) for r in private
                         if r["bhlh_id"] not in UNKNOWN and r["partner_id"] not in UNKNOWN)
    outcome_values = {utility: [number(row, col) for row in outcomes.values()]
                      for utility, col in (("binary", "binary_call"), ("gated", "q_gated"))}
    full = {k:max(vals) if len(vals) == 598 and None not in vals else None for k, vals in outcome_values.items()}
    scored = defaultdict(list)
    for row in read(root / "scores/all_pair_scores.tsv"):
        v = row["variant"]
        d = design[v, row["pair_id"]]
        result = outcomes.get(d["pair_slot_id"], {})
        if (d["bhlh_id"],d["partner_id"],d["bhlh_symbol"],d["partner_symbol"]) != (row["bhlh_id"],row["partner_id"],row["bhlh_symbol"],row["partner_symbol"]):
            raise ValueError("Design identity mismatch across lock")
        both = all(d[k] not in UNKNOWN for k in ("bhlh_id","partner_id"))
        unique = both and identities[v, d["bhlh_id"], d["partner_id"]] == 1
        mapped_column = d.get("raw_s1_column_if_mapped", "")
        raw_field = "raw_s1_pct_s1_thesis" if mapped_column and result.get("raw_s1_cell_s1_thesis","").startswith(mapped_column) else "raw_s1_pct"
        row.update(
            pair_slot_id=d["pair_slot_id"], identity_unique=unique, allowed=d["eligible"] in TRUE and unique,
            design_exclusion_reason=d.get("exclusion_reason",""),
            binary=number(result, "binary_call"), gated=number(result, "q_gated"),
            raw_s1=number(result, raw_field),
            locator=result.get("source_cell_or_bbox", ""), call_basis=result.get("call_basis", ""),
            measurement_status=result.get("measurement_status", ""),
        )
        scored[row["endpoint"], v, row["method"]].append(row)

    endpoints = sorted({e for e, _, _ in scored})
    panels, task_definitions, null_specs, full_panel_order = [], [], {}, {}
    for endpoint in endpoints:
        for variant in variants:
            mrows = {m:scored[endpoint, variant, m] for m in MODELS}
            for utility in ("binary", "gated"):
                def filter_rows(rows):
                    return sorted([r for r in rows if r[utility] is not None and r["allowed"] and r["scoreable"] == "1"],
                                  key=lambda r: int(r["rank"]))
                eligible = {m:filter_rows(mrows[m]) for m in MODELS}
                indices = {m:{r["pair_id"] for r in rows} for m, rows in eligible.items()}
                if indices["network_add"] != indices["network_union"]:
                    raise ValueError("Network eligibility is not shared")
                common = set.intersection(*indices.values())
                any_rows = mrows["expression"]
                panels.append({
                    "endpoint": endpoint, "variant": variant, "utility": utility,
                    "design_slots": len(any_rows), "identity_unique_n": sum(r["identity_unique"] for r in any_rows),
                    "identity_allowed_n": sum(r["allowed"] for r in any_rows), "position_label_n": sum(r[utility] is not None for r in any_rows),
                    "allowed_and_label_n": sum(r["allowed"] and r[utility] is not None for r in any_rows),
                    "expression_panel_n": len(indices["expression"]), "network_panel_n": len(indices["network_add"]),
                    "common_panel_n": len(common),
                })
                for model in MODELS:
                    for panel_type in ("native", "common"):
                        pool = eligible[model] if panel_type == "native" else [r for r in eligible[model] if r["pair_id"] in common]
                        full_panel_order[endpoint,variant,utility,panel_type,model] = [r["pair_id"] for r in pool]
                        by_key = sorted((r["pair_id"], r[utility]) for r in pool)
                        cache_key = sha256(repr((utility, by_key)).encode()).hexdigest()
                        for size in BUDGETS:
                            nd = (cache_key, size)
                            if len(pool) >= size:
                                null_specs[nd] = ([q for _, q in by_key], size, cache_key, utility=="binary")
                            task_definitions.append((nd, model, endpoint, variant, panel_type, utility, pool, size, full[utility]))
    results, partners, tops = [], [], []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        computed = {key:executor.submit(uniform, *defn) for key, defn in null_specs.items()}
        nulls = {k:f.result() for k,f in computed.items()}
        jobs = [executor.submit(evaluate_one, (*d, nulls.get(d[0]))) for d in task_definitions]
        for job in jobs:
            metrics, partner, picks = job.result()
            results.append(metrics)
            if partner:
                partners.append(partner)
            tops.extend(picks)
    counts = Counter(row["binary_call"] for row in records)
    gated = Counter("missing" if number(r,"q_gated") is None else "defined" for r in records)
    comparison = []
    for std in [r for r in results if r["method"] == "network_add"]:
        group = partner_key(std)
        alt = next(r for r in results if partner_key(r) == {**group, "method":"network_union"})
        ids = [r["pair_id"] for r in tops if partner_key(r)==group]
        other = [r["pair_id"] for r in tops if partner_key(r)=={**group, "method":"network_union"}]
        comparison.append({**group, "comparison": "network_add_vs_union", "N":std["N"], "B":std["B"],
                           "identical_order":int(ids == other) if ids else None,
                           "identical_complete_eligible_order":int(full_panel_order[std["endpoint"],std["variant"],std["utility"],std["panel_type"],"network_add"]
                                                                  == full_panel_order[std["endpoint"],std["variant"],std["utility"],std["panel_type"],"network_union"]),
                           "overlap":len(set(ids) & set(other)) if ids else None, "delta_union_over_add":alt["mean_utility"]-std["mean_utility"] if ids else None})
    table(output / "panel_flow.tsv", panels)
    table(output / "utility_results.tsv", results)
    table(output / "partner_null.tsv", partners)
    table(output / "top_pairs.tsv", tops)
    table(output / "sensitivity.tsv", results)
    table(output / "network_comparison.tsv", comparison)
    record(output / "EVALUATION_META.json", {"evaluated_utc":now(), "workers":workers, "position_binary_counts":counts,
        "position_gated_counts":gated, "full_position_oracle":full, "raw_s1_position_n":sum(number(r,"raw_s1_pct") is not None for r in records),
        "utility_rows":len(results),"eligible_metric_rows":len(partners),"uniform_cache_groups":len(null_specs),
        "seed":SEED, "numpy":np.__version__})
    print(f"{len(results)} evaluation rows, {len(partners)} complete budget rows, {len(tops)} selected-pair receipts")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("--workers", type=int, default=min(32, os.cpu_count() or 1))
    args = p.parse_args()
    main(args.root, args.workers)
