"""Build two signed mouse graph snapshots using public evidence and frozen IDs."""

import argparse
import csv
import gzip
import hashlib
import json
import re
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

MD5 = {"CollecTRI_regulons.csv": "cee4a3943c059e6dd8796ced6dec44f6",
       "CollecTRI_source.tsv": "c8622a46a49e21e955182f000cf1c067"}
ZENO = "https://zenodo.org/records/8192729/files/"
EP = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CUTOFF = "2015-12-31"
SCREEN_PMID = "29743677"
HCOP_MIN_SUPPORT = 3
PMID = re.compile(r"(?<!\d)\d{5,9}(?!\d)")
SOURCE_EVIDENCE = {
    "TRRUST": ("Regulation", {"Activation": 1, "Repression": -1}),
    "TFactS": ("Sign", {"UP": 1, "DOWN": -1}),
    "CytReg": ("Activation/Repression", {"Activation": 1, "Repression": -1}),
    "NTNU Curated": ("Sign", {"+": 1, "-": -1}),
    "Pavlidis2021": ("Mode of action", {"activation": 1, "repression": -1}),
    "SIGNOR": ("Effect", {"up-regulates quantity by expression": 1,
                          "down-regulates quantity by repression": -1}),
}


def digest(path, alg="sha256"):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, alg).hexdigest()


def rows(path, delim="\t", skip=0):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        for _ in range(skip):
            next(stream)
        yield from csv.DictReader(stream, delimiter=delim)


def table(path, records, fields=None):
    records = list(records)
    columns = fields or list(records[0])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def jsonfile(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def fetch(path, url, retries=4):
    if path.exists():
        return path.read_bytes()
    for attempt in range(retries):
        try:
            with urlopen(url, timeout=35) as res:
                data = res.read()
            path.write_bytes(data)
            return data
        except OSError as exc:
            if attempt == retries - 1:
                raise
            print(f"Retry {attempt + 1} {path.name} {exc}", flush=True)
            time.sleep(2 ** attempt)


def orthologs(hcop, gene_info):
    """Require unique forward and reverse edges among the high-support proposals."""
    mm = {r["GeneID"]: r for r in rows(gene_info) if r["#tax_id"] == "10090"}
    by_symbol, forward, reverse = defaultdict(set), defaultdict(set), defaultdict(set)
    good_pairs = {}
    for r in rows(hcop):
        hs, hi, mi = (r[k] for k in ("human_symbol", "human_entrez_gene", "mouse_entrez_gene"))
        support = r["support"].split(",")
        if not (hs != "-" and hi.isdigit() and mi.isdigit() and len(set(support)) >= HCOP_MIN_SUPPORT):
            continue
        by_symbol[hs].add(hi)
        forward[hi].add(mi)
        reverse[mi].add(hi)
        key = hi, mi
        if key not in good_pairs or len(good_pairs[key]["support"].split(",")) < len(support):
            good_pairs[key] = r
    assert good_pairs, "HCOP parsing yielded no supported pairs"
    return mm, by_symbol, forward, reverse, good_pairs


def map_symbol(symbol, mapping):
    mm, by_symbol, forward, reverse, good_pairs = mapping
    hh = by_symbol[symbol]
    hi = next(iter(hh)) if len(hh) == 1 else ""
    dests = forward[hi] if hi else set()
    mi = next(iter(dests)) if len(dests) == 1 else ""
    item = good_pairs.get((hi, mi), {})
    record = dict(human_symbol=symbol, human_entrez=";".join(sorted(hh)),
                  human_ensembl=item.get("human_ensembl_gene", ""),
                  source_hgnc=item.get("hgnc_id", ""),
                  mouse_entrez=";".join(sorted(dests)), mouse_ensembl=item.get("mouse_ensembl_gene", ""),
                  mouse_hcop_symbol=item.get("mouse_symbol", ""),
                  mouse_symbol=mm.get(mi, {}).get("Symbol", ""),
                  HCOP_support=item.get("support", ""), map_status="")
    if not hh:
        record["map_status"] = "no_high_support_HCOP_human_symbol"
    elif not hi:
        record["map_status"] = "ambiguous_human_symbol"
    elif not mi:
        record["map_status"] = "ambiguous_HCOP_forward"
    elif reverse[mi] != {hi}:
        record["map_status"] = "ambiguous_HCOP_reverse"
    elif mi not in mm:
        record["map_status"] = "mouse_entrez_absent_in_frozen_gene_info"
    else:
        # Stable NCBI IDs are checked with at least one name or cross-reference.
        m = mm[mi]
        ok = (m["Symbol"] == record["mouse_hcop_symbol"]
              or "Ensembl:" + record["mouse_ensembl"] in m["dbXrefs"].split("|")
              or "MGI:" + item["mgi_id"] in m["dbXrefs"].split("|"))
        record["map_status"] = "one_to_one" if ok else "mouse_identity_disagrees"
    return record


def pmids(txt):
    return set(PMID.findall(txt))


def direction_provenance(source_row):
    """Only paired citation/sign tokens in the corresponding signed field."""
    for resource, (direction_field, signs) in SOURCE_EVIDENCE.items():
        directions = source_row[f"[{resource}] {direction_field}"].split("|")
        refs = source_row[f"[{resource}] PMID"].split("|")
        if len(directions) != len(refs):
            continue
        for d, pp in zip(directions, refs):
            sign = signs.get(d)
            if sign:
                for p in sorted(pmids(pp)):
                    yield resource, d, sign, p


def get_dates(folder, required, archived_dates):
    dates = {r["pmid"]: (r["earliest_interval_end"], "frozen_TRRUST_date_" + r["date_basis"])
             for r in rows(archived_dates) if re.fullmatch(r"\d{4}-\d{2}-\d{2}", r["earliest_interval_end"])}
    remaining = sorted(required - dates.keys(), key=int)
    print(f"Directional PMIDs {len(required)}, frozen dated {len(required & dates.keys())}, "
          f"Europe PMC lookup {len(remaining)}", flush=True)
    date_cache = folder / "europepmc_dates"
    date_cache.mkdir(exist_ok=True)

    def date_chunk(ids):
        query = " OR ".join(f"EXT_ID:{i}" for i in ids)
        url = EP + "/search?" + urlencode(dict(query=f"({query}) AND SRC:MED",
                                               format="json", resultType="core", pageSize=len(ids)))
        key = hashlib.sha256(query.encode()).hexdigest()[:16]
        try:
            data = json.loads(fetch(date_cache / f"{key}.json", url))["resultList"]["result"]
            return {r["id"]: (r["firstPublicationDate"], "EuropePMC_firstPublicationDate")
                    for r in data if r["id"] in ids and r.get("firstPublicationDate")}
        except (OSError, ValueError, KeyError) as exc:
            print(f"Date chunk {key} unavailable: {exc}", flush=True)
            return {}

    chunks = [remaining[i:i + 80] for i in range(0, len(remaining), 80)]
    with ThreadPoolExecutor(max_workers=3) as workers:
        for result in workers.map(date_chunk, chunks):
            dates.update(result)
    table(folder / "PMID_DATES.tsv", (
        dict(pmid=p, first_public_by=dates.get(p, ("", ""))[0], date_source=dates.get(p, ("", ""))[1],
             cutoff=("on_or_before" if dates[p][0] <= CUTOFF else "after") if p in dates else "date_unresolved")
        for p in sorted(required, key=int)))
    return dates


def write_edges(scoring, prior, candidates):
    kept, decisions = {}, {}
    for key, es in sorted(candidates.items()):
        possible = {e["sign"] for e in es}
        if len(possible) > 1:
            decisions[key] = "conflicting_projected_signs"
            continue
        e = es[0]
        kept[key] = dict(tf_symbol=e["tf_symbol"], target_symbol=e["target_symbol"],
                         sign=e["sign"], tf_gene_id=key[0], target_gene_id=key[1],
                         source_release_records=";".join(str(i["release_line"]) for i in es),
                         supporting_pmids=";".join(sorted(set().union(*(i["pmids"] for i in es)), key=int)),
                         source_context=";".join(sorted(set().union(*(i["evidence_resources"] for i in es)))))
        decisions[key] = "retained" if len(es) == 1 else "retained_deduplicated"
    names = ("tf_symbol", "target_symbol", "sign", "tf_gene_id", "target_gene_id",
             "source_release_records", "supporting_pmids", "source_context")
    table(scoring / "inputs" / (prior + ".tsv"), kept.values(), names)
    return kept, decisions


def main(scoring, old, utility, audit, modern_only=False):
    source = scoring / "sources"
    hcop = audit / "sources/annotation/human_mouse_hcop_fifteen_column.txt.gz"
    mm = utility / "expression_raw/Mus_musculus.gene_info.gz"
    dates_path = utility / "inputs_scoring/TRRUST_PUBMED_DATES.tsv"
    assets = []
    for name, expected in MD5.items():
        f = source / name
        fetch(f, ZENO + name)
        actual = digest(f, "md5")
        if actual != expected:
            raise ValueError(f"Official MD5 mismatch {f} {actual}")
        assets.append(dict(path=str(f), source=ZENO + name, md5=actual, sha256=digest(f),
                           bytes=f.stat().st_size, version="Zenodo 8192729 / v2.0"))
    citation_path = source / "Tsunemoto_citing_EuropePMC.json"
    citation_url = f"{EP}/MED/{SCREEN_PMID}/citations?format=json&page=1&pageSize=1000"
    citation_data = json.loads(fetch(citation_path, citation_url))
    cited_by = {r["id"]: r for r in citation_data["citationList"]["citation"] if r.get("source") == "MED"}
    assert len(citation_data["citationList"]["citation"]) == citation_data["hitCount"]
    for f, version, ref in (
        (citation_path, "retrieved snapshot", citation_url),
        (hcop, "frozen HCOP 2026-10-01", "https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz"),
        (mm, "existing Tsunemoto frozen Entrez snapshot", "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz"),
        (dates_path, "old verified publication dates", str(dates_path)),
    ):
        assets.append(dict(path=str(f), source=ref, md5=digest(f, "md5"), sha256=digest(f),
                           bytes=f.stat().st_size, version=version))
    table(scoring / "PRIOR_SOURCES.tsv", assets)
    network = list(rows(source / "CollecTRI_regulons.csv", delim=","))
    source_rows = {(r["Transcription Factor (Associated Gene Name)"], r["Target Gene (Associated Gene Name)"]): (i, r)
                   for i, r in enumerate(rows(source / "CollecTRI_source.tsv", skip=1), start=3)}
    symbols = {r[k] for r in network for k in ("source", "target")}
    orthology = orthologs(hcop, mm)
    mapped = {s: map_symbol(s, orthology) for s in sorted(symbols)}
    table(source / "ORTHOLOGY_DECISIONS.tsv", mapped.values())
    hits = set()
    evidence, pending = [], []
    date_needed = set()
    for line, r in enumerate(network, start=2):
        a, b = r["source"], r["target"]
        i, evidence_table = source_rows.get((a, b), ("", {}))
        refs = pmids(r["references"])
        source_refs = {p for k, v in evidence_table.items() if "PMID" in k for p in pmids(v)}
        found = (refs | source_refs) & ({SCREEN_PMID} | cited_by.keys())
        hits.update(found)
        x, y = mapped[a], mapped[b]
        mstatus = "pass" if x["map_status"] == y["map_status"] == "one_to_one" else "mapping_excluded"
        status = "screen_direct_citation_removed" if SCREEN_PMID in found else mstatus
        if found - {SCREEN_PMID}:
            status = "citing_work_edge_removed_pending_independence"
        weight = float(r["weight"])
        assert weight in (-1.0, 1.0)
        signed = sorted(set(direction_provenance(evidence_table))) if status == "pass" and evidence_table else []
        date_needed.update(sign_tuple[3] for sign_tuple in signed)
        pending.append((dict(release_line=line, source_line=i, human_tf=a, human_tg=b,
                             human_tf_entrez=x["human_entrez"], human_tg_entrez=y["human_entrez"],
                             source_hgnc=x["source_hgnc"], target_hgnc=y["source_hgnc"],
                             tf_mapping=x["map_status"], target_mapping=y["map_status"],
                             tf_gene_id=x["mouse_entrez"], target_gene_id=y["mouse_entrez"],
                             tf_symbol=x["mouse_symbol"], target_symbol=y["mouse_symbol"],
                             tf_hcop_support=x["HCOP_support"], target_hcop_support=y["HCOP_support"],
                             weight=weight, sign=int(weight), release_sign_basis=r["sign_decision"],
                             resources=r["resources"], references=r["references"], contaminated_pmids=";".join(sorted(found)),
                             resource_pmids=";".join(sorted(source_refs, key=int)), mapping_and_exposure_status=status),
                        refs | source_refs, signed))
    jsonfile(source / "SCREEN_PROVENANCE_CHECK.json", dict(
        direct_screen_pmid=SCREEN_PMID, indexed_citing_records=citation_data["hitCount"],
        matched_pmids=sorted(hits), matched_citing_metadata=[cited_by[p] for p in sorted(hits) if p in cited_by],
        action="Conservative entire edge exclusion for attributable direct citation or citation-index overlap",
        limitation="Non-indexed citations and indirect reuse without citation remain unresolved"))
    if modern_only:
        snapshot = defaultdict(list)
        for r, refs, _ in pending:
            if r["mapping_and_exposure_status"] == "pass":
                snapshot[r["tf_gene_id"], r["target_gene_id"]].append(
                    {**r, "pmids": refs, "evidence_resources": set(r["resources"].split(";"))})
        graph, _ = write_edges(scoring, "collectri_modern", snapshot)
        print(f"Modern signed graph ready: {len(graph)} edges", flush=True)
        return
    dates = get_dates(source, date_needed, dates_path)
    layers = {"collectri_modern": defaultdict(list), "collectri_through2015": defaultdict(list)}
    for record, refs, signed in pending:
        key = record["tf_gene_id"], record["target_gene_id"]
        status = record["mapping_and_exposure_status"]
        early = []
        for resource, mode, sign, p in signed:
            pub, basis = dates.get(p, ("", ""))
            decision = "date_unresolved" if not pub else "early" if pub <= CUTOFF else "after_cutoff"
            if decision == "early":
                early.append((resource, sign, p))
            evidence.append(dict(release_line=record["release_line"], source_line=record["source_line"],
                                 human_tf=record["human_tf"], human_tg=record["human_tg"], tf_gene_id=key[0],
                                 target_gene_id=key[1], resource=resource, source_mode=mode, sign=sign,
                                 pmid=p, pubdate=pub, date_basis=basis, decision=decision))
        signs = {x[1] for x in early}
        historic = status if status != "pass" else "no_linked_explicit_signed_source_citation" if not signed else (
            "no_dated_early_direction" if not early else "early_sign_conflict" if len(signs) > 1 else "pass")
        record.update(historical_eligibility=historic,
                      historic_signed_pmids=";".join(sorted({x[2] for x in early}, key=int)),
                      historic_sign=next(iter(signs)) if len(signs) == 1 else "")
        if status == "pass":
            layers["collectri_modern"][key].append({**record, "pmids": refs,
                                                    "evidence_resources": set(record["resources"].split(";"))})
        if historic == "pass":
            layers["collectri_through2015"][key].append(
                {**record, "sign": next(iter(signs)), "pmids": {x[2] for x in early},
                 "evidence_resources": {x[0] for x in early}})
    trrust = {(r["tf_gene_id"], r["target_gene_id"]): r for r in rows(old / "inputs/TRRUST_MOUSE_SIGNED.tsv")}
    kept_decisions = {p: write_edges(scoring, p, es) for p, es in layers.items()}
    for record, _, _ in pending:
        key = record["tf_gene_id"], record["target_gene_id"]
        for p, (_, d) in kept_decisions.items():
            state = record["mapping_and_exposure_status"] if p == "collectri_modern" else record["historical_eligibility"]
            record[p + "_decision"] = d.get(key, "projection_missing") if state == "pass" else state
    table(scoring / "PRIOR_PROVENANCE.tsv", (r for r, _, _ in pending))
    table(source / "DIRECTIONAL_EVIDENCE.tsv", evidence)

    summary = []
    overlaps = []
    for p, (edges, _) in kept_decisions.items():
        common = edges.keys() & trrust.keys()
        overlap_agree = sum(edges[k]["sign"] == int(trrust[k]["sign"]) for k in common)
        summary.append(dict(prior=p, edges=len(edges), nodes=len({x for k in edges for x in k}),
                            activating=sum(r["sign"] == 1 for r in edges.values()), shared=len(common),
                            signed_agree=overlap_agree, signed_disagree=len(common) - overlap_agree,
                            trrust_only=len(trrust) - len(common), collectri_only=len(edges) - len(common),
                            signed_citation_pmids=len(date_needed), dated_pmids=len(date_needed & dates.keys())))
        for a, b in sorted(edges.keys() | trrust.keys()):
            x, y = trrust.get((a, b), {}), edges.get((a, b), {})
            state = ("shared_sign_agree" if int(x["sign"]) == int(y["sign"]) else "shared_sign_disagree") if x and y else (
                "collectri_only" if y else "trrust_only")
            overlaps.append(dict(prior=p, tf_gene_id=a, target_gene_id=b, trrust_sign=x.get("sign", ""),
                                 collectri_sign=y.get("sign", ""), relation=state))
    table(scoring / "PRIOR_GRAPH_SUMMARY.tsv", summary)
    table(scoring / "PRIOR_EDGE_OVERLAP.tsv", overlaps)
    print(Counter(r["historical_eligibility"] for r, _, _ in pending), flush=True)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for name in ("scoring", "old", "utility", "audit"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--modern-only", action="store_true")
    args = parser.parse_args()
    main(args.scoring, args.old, args.utility, args.audit, modern_only=args.modern_only)
