"""Download frozen mouse TRRUST and PubMed records using only official endpoints."""

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import time
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from xml.etree import ElementTree as ET


RELEASE_URL = "https://www.grnpedia.org/trrust/downloadnetwork.php"
MOUSE_URL = "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv"
BIOC_URL = "https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.bioc.xml"
ARTICLE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC5753191/?report=classic"
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CUTOFF = date(2015, 12, 31)
SIGNED = {"Activation": 1, "Repression": -1}
DATE_STATUSES = {"epublish", "ppublish", "aheadofprint"}
MONTHS = {v.lower(): i for i, v in enumerate((
    "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
), 1)}
# Deliberately adverse-evidence rule: one source row in BioC can associate a
# list of PMIDs with a shared MoR, without separately assigning that sign
# to each citation. Mixed-date rows need specific qualifying evidence.
CACHE_BATCH_SIZE = 150


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fetch(url, destination, timeout=90):
    if destination.exists():
        return destination.read_bytes()
    tmp = destination.with_suffix(destination.suffix + ".partial")
    subprocess.run([
        "curl", "-sSL", "--fail", "--retry", "2", "--retry-delay", "1",
        "--max-time", str(timeout), url, "-o", str(tmp),
    ], check=True, stdout=subprocess.DEVNULL)
    data = tmp.read_bytes()
    if not data:
        raise ValueError(f"Empty HTTP response at {url}")
    tmp.replace(destination)
    return data


def raw_tsv(data):
    rows = [line.split("\t") for line in data.decode("utf-8-sig").splitlines() if line]
    if not rows or any(len(row) != 4 or row[2] not in SIGNED | {"Unknown": 0}
                       or not re.fullmatch(r"[0-9]+(?:;[0-9]+)*", row[3]) for row in rows):
        raise ValueError("Unexpected mouse TSV format")
    return [(*row[:3], tuple(row[3].split(";"))) for row in rows]


def read_original(folder):
    html = fetch(RELEASE_URL, folder / "trrust_release_20180416.html")
    if b"16/04/2018" not in html or b"6,490" not in html or b"mouse.tsv" not in html:
        raise ValueError("Unexpected release page")
    tsv = fetch(MOUSE_URL, folder / "trrust_rawdata.mouse.tsv")
    bioc = fetch(BIOC_URL, folder / "trrust_rawdata.mouse.bioc.xml")
    rows = raw_tsv(tsv)
    root = ET.fromstring(bioc)
    if root.tag != "collection" or root.findtext("key") != "BioCTRRUST":
        raise ValueError("Unexpected BioC XML format")
    passages = {}
    for passage in root.findall(".//passage"):
        all_pmids = tuple(passage.findtext("infon[@key='article-id_pmid']").split(";"))
        g1 = passage.find("annotation[@id='G1']")
        g2 = passage.find("annotation[@id='G2']")
        mode = passage.findtext("annotation[@id='E1']/text")
        key = g1.findtext("text"), g2.findtext("text"), mode, all_pmids
        ids = g1.findtext("infon[@key='GeneID']"), g2.findtext("infon[@key='GeneID']")
        if key in passages and passages[key] != ids:
            raise ValueError(f"Conflicting BioC gene IDs for {key}")
        passages[key] = ids
    if set(passages) - set(rows):
        raise ValueError("Official BioC XML contains unexpected rows outside TSV")
    return rows, passages


def fetch_pubmed(pmids, folder):
    records = {}
    manifests = []
    for start in range(0, len(pmids), CACHE_BATCH_SIZE):
        batch = pmids[start:start + CACHE_BATCH_SIZE]
        name = f"pubmed_{start // CACHE_BATCH_SIZE + 1:03d}.xml"
        params = {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
        url = NCBI + "?" + urlencode(params, safe=",")
        compressed = folder / (name + ".gz")
        url_record = folder / (name + ".url")
        if compressed.exists():
            data = gzip.decompress(compressed.read_bytes())
            # The first acquisition attempt used percent-escaped commas.
            url = url_record.read_text().strip() if url_record.exists() else NCBI + "?" + urlencode(
                {**params, "tool": "tsunemoto_historical_prior_audit"})
        else:
            path = folder / name
            data = fetch(url, path)
            root = ET.fromstring(data)
            if root.tag != "PubmedArticleSet" or root.find("./ERROR") is not None:
                raise ValueError(f"Unexpected EFetch payload: {name}")
            zip_tmp = compressed.with_suffix(".gz.partial")
            zip_tmp.write_bytes(gzip.compress(data, mtime=0))
            zip_tmp.replace(compressed)
            url_record.write_text(url + "\n", encoding="ascii")
            path.unlink()
            time.sleep(0.45)
        root = ET.fromstring(data)
        if root.tag != "PubmedArticleSet" or root.find("./ERROR") is not None:
            raise ValueError(f"Unexpected EFetch payload: {name}")
        for article in root:
            pmid = article.findtext("MedlineCitation/PMID") or article.findtext("BookDocument/PMID")
            if pmid in batch:
                records[pmid] = article
        manifests.append((url, compressed, sha(data), len(batch), sum(p in records for p in batch)))
        print(f"{name}: {len(records)}/{len(pmids)} articles available", flush=True)
    return records, manifests


def element_text(node):
    if node is None:
        return ""
    return " ".join("".join(node.itertext()).split())


def date_interval(node):
    """Return a conservative interval without inventing a missing month/day."""
    year = node.findtext("Year")
    if not year:
        raw = node.findtext("MedlineDate") or ""
        years = [int(y) for y in re.findall(r"\b(?:18|19|20)\d{2}\b", raw)]
        if not years:
            return None
        return date(min(years), 1, 1), date(max(years), 12, 31), raw, "medline_year_or_range"
    try:
        year_int = int(year)
        month = node.findtext("Month") or ""
        day = node.findtext("Day") or ""
        month_int = int(month) if month.isdigit() else MONTHS.get(month[:3].lower())
        if month_int is None:
            return date(year_int, 1, 1), date(year_int, 12, 31), year, "year"
        if not day:
            next_month = date(year_int + (month_int == 12), 1 if month_int == 12 else month_int + 1, 1)
            first = date(year_int, month_int, 1)
            return first, date.fromordinal(next_month.toordinal() - 1), first.strftime("%Y-%m"), "month"
        exact = date(year_int, month_int, int(day))
        return exact, exact, exact.isoformat(), "day"
    except ValueError:
        return None


def publication(article):
    article_prefix = "MedlineCitation/Article/"
    dates = [
        (node, f"ArticleDate:{node.get('DateType', '')}")
        for node in article.findall(article_prefix + "ArticleDate")
        if node.get("DateType", "").lower() in {"electronic", "print"}
    ]
    dates.extend((n, "JournalIssue/PubDate") for n in article.findall(article_prefix + "Journal/JournalIssue/PubDate"))
    dates.extend((n, "Book/PubDate") for n in article.findall("BookDocument/Book/PubDate"))
    dates.extend(
        (n, f"PubMedHistory:{n.get('PubStatus')}")
        for n in article.findall("PubmedData/History/PubMedPubDate")
        if n.get("PubStatus") in DATE_STATUSES
    )
    known = [(*interval, kind) for node, kind in dates if (interval := date_interval(node))]
    if not known:
        return ("", "", "", "", "", "date_unresolved")
    # The earliest reconstructed lower bound need not be exact. A qualifying
    # date requires that the entire interval for at least one source field
    # lie before the threshold. Prefer the earliest verifiable upper bound.
    lower, upper, reported, precision, basis = min(known, key=lambda value: (value[1], value[0], value[4]))
    state = "at_or_before_cutoff" if upper <= CUTOFF else (
        "after_cutoff" if all(value[0] > CUTOFF for value in known) else "cutoff_unresolved")
    return (reported, precision, lower.isoformat(), upper.isoformat(), basis, state)


def publication_records(pmids, records):
    out = []
    for pmid in pmids:
        article = records.get(pmid)
        if article is None:
            out.append({"pmid": pmid, "public_date_reported": "", "date_precision": "",
                        "earliest_interval_start": "", "earliest_interval_end": "", "date_basis": "",
                        "cutoff_status": "pubmed_unavailable", "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
            continue
        fields = ("public_date_reported", "date_precision", "earliest_interval_start",
                  "earliest_interval_end", "date_basis", "cutoff_status")
        out.append({"pmid": pmid, **dict(zip(fields, publication(article))),
                    "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
    return out


def put_tsv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def get_file_info(path):
    data = path.read_bytes()
    return {
        "path": str(path.resolve()), "filename": path.name, "bytes": len(data), "sha256": sha(data),
        "download_time_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
    }


def main(root, fetch_only):
    outputs = root / "inputs_scoring"
    source = outputs / "prior_sources"
    source.mkdir(parents=True, exist_ok=True)
    rows, passages = read_original(source)
    article_filename = "PMC5753191_trrust_methods.html"
    article_html = fetch(ARTICLE_URL, source / article_filename)
    if b"gkx1013" not in article_html or b"23 577" not in article_html:
        raise ValueError("Unexpected primary TRRUST article HTML")
    pmids = sorted({pmid for _, _, _, group in rows for pmid in group}, key=int)
    records, batches = fetch_pubmed(pmids, source)
    pmid_rows = publication_records(pmids, records)
    put_tsv(outputs / "TRRUST_PUBMED_DATES.tsv", pmid_rows, list(pmid_rows[0]))
    by_pmid = {r["pmid"]: r for r in pmid_rows}
    out_manifest = []
    for url, filename, role in (
        (RELEASE_URL, "trrust_release_20180416.html", "release"),
        (MOUSE_URL, "trrust_rawdata.mouse.tsv", "mouse_raw_rows"),
        (BIOC_URL, "trrust_rawdata.mouse.bioc.xml", "mouse_BioC_Entrez_and_sign_groups"),
    ):
        out_manifest.append({"resource": role, "url": url, **get_file_info(source / filename),
                             "version": "TRRUST v2 release 2018-04-16", "license_or_access": "public; non-commercial research",
                             "content_sha256": "", "requested_records": "", "verified_records": ""})
    out_manifest.append({"resource": "original_TRRUST_v2_curation_and_early_counts", "url": ARTICLE_URL,
                         **get_file_info(source / article_filename),
                         "version": "PMC5753191; doi:10.1093/nar/gkx1013; original 2017 journal article",
                         "license_or_access": "open access; CC BY-NC 4.0",
                         "content_sha256": "", "requested_records": "", "verified_records": ""})
    for url, path, content, requested, verified in batches:
        out_manifest.append({"resource": "pubmed_official_efetch_batch_gzip", "url": url, **get_file_info(path),
                             "version": "PubMed XML DTD 2025-01-01; see payload",
                             "license_or_access": "public PubMed bibliographic metadata; abstracts not redistributed",
                             "content_sha256": content, "requested_records": requested, "verified_records": verified})
    put_tsv(root / "SOURCE_TRRUST.tsv", out_manifest, list(out_manifest[0]))
    mixed = []
    for tf, target, sign, group in rows:
        if sign not in SIGNED:
            continue
        qualifiers = [p for p in group if by_pmid[p]["cutoff_status"] == "at_or_before_cutoff"]
        if qualifiers and len(qualifiers) != len(group):
            for pmid in qualifiers:
                article = records[pmid]
                mixed.append({
                    "tf_symbol": tf, "target_symbol": target, "mode": sign, "pmid": pmid,
                    "all_pmids": ";".join(group), "pmid_url": by_pmid[pmid]["source_url"],
                    "title": element_text(article.find("MedlineCitation/Article/ArticleTitle")),
                    "abstract": element_text(article.find("MedlineCitation/Article/Abstract")),
                })
    put_tsv(source / "mixed_evidence_abstracts.tsv", mixed,
            ["tf_symbol", "target_symbol", "mode", "pmid", "all_pmids", "pmid_url", "title", "abstract"])
    manual_file = outputs / "TRRUST_MIXED_SIGN_REVIEW.tsv"
    if not mixed and not manual_file.exists():
        put_tsv(manual_file, [], ["tf_symbol", "target_symbol", "mode", "pmid", "verdict", "note", "source_url"])
    print(json.dumps({"raw_rows": len(rows), "distinct_raw_edges": len({r[:2] for r in rows}),
                      "signed_raw_rows": sum(row[2] in SIGNED for row in rows), "pubmed_requested": len(pmids),
                      "pubmed_fetched": len(records), "date_funnel": Counter(r["cutoff_status"] for r in pmid_rows),
                      "mixed_sign_groups": len({(r["tf_symbol"], r["target_symbol"], r["mode"], r["all_pmids"]) for r in mixed}),
                      "mixed_early_pmid_links": len(mixed)}, ensure_ascii=False))
    if not fetch_only:
        build_prior(root, rows, passages, by_pmid)


def build_prior(root, rows, passages, by_pmid):
    outputs = root / "inputs_scoring"
    manual_file = outputs / "TRRUST_MIXED_SIGN_REVIEW.tsv"
    with manual_file.open(encoding="utf-8-sig", newline="") as handle:
        checks = {
            (row["tf_symbol"], row["target_symbol"], row["mode"], row["pmid"]): row
            for row in csv.DictReader(handle, delimiter="\t")
        }
    mode_rows = []
    candidates = defaultdict(list)
    consulted_reviews = set()
    for index, (tf, target, mode, group) in enumerate(rows, 1):
        ids = passages.get((tf, target, mode, group), ("", ""))
        date_eligible = [p for p in group if by_pmid[p]["cutoff_status"] == "at_or_before_cutoff"]
        manual = {}
        for pmid in date_eligible:
            key = tf, target, mode, pmid
            if key in checks:
                manual[pmid] = checks[key]
                consulted_reviews.add(key)
        item = {
            "source_tsv_row": index, "tf_symbol": tf, "target_symbol": target, "source_mode": mode,
            "tf_gene_id": ids[0], "target_gene_id": ids[1],
            "raw_pmids": ";".join(group),
            "precutoff_pmids": ";".join(date_eligible),
            "publication_exclusions": ";".join(f'{p}:{by_pmid[p]["cutoff_status"]}' for p in group if p not in date_eligible),
            "mixed_review_pmids": ";".join(manual),
            "verified_pmids": "", "status": "", "basis": "",
        }
        eligible_pmids = []
        if mode not in SIGNED:
            item["status"] = "excluded_unsigned"
            item["basis"] = "TRRUST source mode Unknown"
        elif not all(value and value.isdigit() for value in ids):
            item["status"] = "excluded_bioc_gene_id_unresolved"
            item["basis"] = "Official BioC lacks an exact signed passage with both gene IDs"
        elif not date_eligible:
            item["status"] = "excluded_no_pre_cutoff_support"
            item["basis"] = "No source PMID has an auditable pre-cutoff publication date"
        elif len(date_eligible) == len(group):
            eligible_pmids = date_eligible
            item["basis"] = "All official PMIDs in signed TSV/BioC passage meet cutoff"
        else:
            if set(manual) != set(date_eligible):
                raise ValueError(f"Early primary evidence review is incomplete for {(tf, target, mode, group)}")
            eligible_pmids = [p for p in date_eligible if manual[p]["verdict"] == "early_sign_supported"]
            if eligible_pmids:
                item["basis"] = "Direct support verified in pre-cutoff primary article"
            else:
                item["status"] = "excluded_mixed_source_sign_unresolved"
                item["basis"] = "Neither BioC nor inspected early primary evidence assigns the signed mode unambiguously"
        if eligible_pmids:
            item["verified_pmids"] = ";".join(eligible_pmids)
            item["status"] = "eligible_for_edge_deduplication"
            candidates[tf, target].append((item, SIGNED[mode]))
        mode_rows.append(item)
    if set(checks) != consulted_reviews:
        raise ValueError("Manual review file contains unrelated rows")

    output = []
    for (tf, target), supporting in sorted(candidates.items()):
        signs = {sign for _, sign in supporting}
        ids = {(item["tf_gene_id"], item["target_gene_id"]) for item, _ in supporting}
        if len(signs) != 1 or len(ids) != 1:
            status = "excluded_conflicting_sign" if len(signs) != 1 else "excluded_conflicting_gene_id"
            for item, _ in supporting:
                item["status"] = status
            continue
        (tf_id, target_id), = ids
        sign, = signs
        pmids = {pmid for item, _ in supporting for pmid in item["verified_pmids"].split(";")}
        output.append({
            "tf_symbol": tf, "target_symbol": target, "sign": sign,
            "tf_gene_id": tf_id, "target_gene_id": target_id,
            "verified_pmids": ";".join(sorted(pmids, key=int)),
            "source_tsv_rows": ";".join(str(item["source_tsv_row"]) for item, _ in supporting),
        })
        for idx, (item, _) in enumerate(supporting):
            item["status"] = "retained" if idx == 0 else "retained_collapsed_support"
    if not output:
        raise ValueError("No unambiguous signed mouse edges remain")
    put_tsv(outputs / "TRRUST_MOUSE_SIGNED.tsv", output, list(output[0]))
    put_tsv(outputs / "TRRUST_EDGE_AUDIT.tsv", mode_rows, list(mode_rows[0]))
    tfs = {r["tf_symbol"] for r in output}
    tgs = {r["target_symbol"] for r in output}
    raw_tfs = {r[0] for r in rows}
    raw_tgs = {r[1] for r in rows}
    signed_counts = Counter(mode for _, _, mode, _ in rows)
    degrees = Counter(row["tf_symbol"] for row in output)
    version = {
        "research_article": "Han et al. doi:10.1093/nar/gkx1013",
        "article_original_mouse_TFs_2017": 828, "article_original_mouse_edges_2017": 6552,
        "release_date": "2018-04-16", "release_listed_mouse_TFs": 827, "release_listed_mouse_edges": 6490,
        "observed_export_rows": len(rows), "observed_export_pairs": len({r[:2] for r in rows}),
        "observed_export_TFs": len(raw_tfs), "observed_export_non_TF_targets": len(raw_tgs - raw_tfs),
        "observed_export_mode_counts": signed_counts,
        "observed_tsv_rows_missing_exact_bioc": sum(row not in passages for row in rows),
        "cutoff": CUTOFF.isoformat(),
        "cutoff_date_fields": ["electronic and print ArticleDate", "JournalIssue PubDate or MedlineDate",
                               "PubMed History epublish, ppublish, aheadofprint", "Book/PubDate"],
        "cutoff_precision_rule": "An exact release or the full month/year interval must end by cutoff. Acceptance, indexing, and revised dates are excluded. Earliest verifiable upper endpoint is retained as evidence.",
        "PMID_counts": Counter(r["cutoff_status"] for r in by_pmid.values()),
        "row_counts_by_decision": Counter(row["status"] for row in mode_rows),
        "mixed_early_source_reviews": len(checks),
        "mixed_early_source_supported": sum(r["verdict"] == "early_sign_supported" for r in checks.values()),
        "frozen_signed_edges": len(output), "frozen_signed_nodes": len(tfs | tgs),
        "frozen_TFs": len(tfs), "frozen_targets_including_TFs": len(tgs),
        "frozen_sign_counts": Counter(row["sign"] for row in output),
        "frozen_PMID_support_count": len({p for row in output for p in row["verified_pmids"].split(";")}),
        "frozen_top_10_source_TF_outdegree_edge_fraction": sum(degrees.values() if len(degrees) < 10 else
            (count for _, count in degrees.most_common(10))) / len(output),
        "signed_file_sha256": sha((outputs / "TRRUST_MOUSE_SIGNED.tsv").read_bytes()),
        "source_manifest": "SOURCE_TRRUST.tsv",
    }
    (outputs / "TRRUST_AUDIT_SUMMARY.json").write_text(json.dumps(version, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(version, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--fetch-only", action="store_true")
    args = parser.parse_args()
    main(args.root, args.fetch_only)
