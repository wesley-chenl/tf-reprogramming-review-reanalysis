"""Reconstruct the registered Tsunemoto screen independently of the scorer.

Record the literal Fig. 1b grid, the raw quantitative S1 cells, source-corroborated
identities, and the prespecified ambiguity axes. The clean handoff contains only
design identities and opaque keys. A color-only extraction is independent of the
manual numeric transcription and is repeated on the NIHMS figure rendering.
"""
from collections import Counter
from pathlib import Path
import csv
import os
import gzip
import hashlib
import itertools
import json
import xml.etree.ElementTree as ET

from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
from PIL import Image
import numpy as np

ROOT = Path(os.environ["PAIR_STAGE"]).resolve()
PRIVATE = ROOT / "sources_private"
PUBLIC = ROOT / "design_public"
FIG_COLUMNS = ["Ascl1", "Ascl2", "Ascl3", "Ascl5", "Ngn1", "Ngn2", "Ngn3", "ND1", "ND2", "ND4", "ND6", "Atoh1", "Atoh7", "Atoh8", "Myf5", "Ptf1a"]
FIG_SYMBOLS = ["Ascl1", "Ascl2", "Ascl3", "Ascl5", "Neurog1", "Neurog2", "Neurog3", "Neurod1", "Neurod2", "Neurod4", "Neurod6", "Atoh1", "Atoh7", "Atoh8", "Myf5", "Ptf1a"]
FIG_ROWS = ["Brn3c", "Brn3a", "Brn3b", "Brn2", "Brn4", "Brn1", "Oct6", "Pit1", "Oct4", "Oct1", "RPF", "Brn5", "Nurr1"]
PARTNERS = dict(zip(FIG_ROWS, ["Pou4f3", "Pou4f1", "Pou4f2", "Pou3f2", "Pou3f4", "Pou3f3", "Pou3f1", "Pou1f1", "Pou5f1", "Pou2f1", "Pou6f2", "Pou6f1", "Nr4a2"]))
# Source Data Fig. 1 is panel d. The following is a literal human reading of Fig. 1b JPG.
FIG_VALUES = [
    "5.7 3.0 0.02 0.95 0.42 0.29 2.0 0.16 0.31 - 0.02 0.30 0.04 - 0.05 -",
    "3.2 1.8 0.02 0.29 0.33 0.32 0.82 0.03 0.63 0.06 0.02 0.26 0.11 0.01 - -",
    "2.8 2.7 - - 0.21 0.28 0.97 0.20 0.33 - - 0.05 - - 0.04 0.01",
    "1.6 1.4 - - 0.76 0.56 1.8 0.07 0.02 - - 0.01 - - - -",
    "0.20 0.15 - - 0.36 0.13 3.0 0.04 - - - - - - - -",
    "- - - - 0.66 0.09 0.15 0.02 - - - - - - - -",
    "0.03 - - - 0.13 0.02 0.11 0.01 - - - - - - - -",
    "0.13 0.25 - - 0.02 0.03 0.34 0.06 - - - - - - - -",
    "1.1 0.55 - - 0.06 0.05 0.84 0.07 0.01 - - 0.02 - - - -",
    "- - - - - - - - - - - - - - - -",
    "- - - - - - - - - - - - - - - -",
    "- - - - - - - - - - - - - - - -",
    "0.05 0.03 - - - - - - - - - - - - - -",
]
RAW_ROWS = dict(zip(FIG_ROWS, [23, 24, 25, 26, 27, 33, 34, 29, 28, 36, 35, 30, 38]))
PARTNER_ROSTER_CELLS = dict(zip(FIG_ROWS, ["H10","H8","H9","G9","G10","G8","G11","H17","G14","G17","H14","H13","I7"]))
ALIAS = {
    "Hen1": "Nhlh1", "Hen2": "Nhlh2", "cMyc": "Myc",
    "MycL": "Mycl", "MycN": "Mycn", "ArntL": "Bmal1",
    "Mdx3": "", "Ascl3": "Ascl3", "Ascl4": "Ascl4",
}
CONFLICT_AXES = {"Ascl3", "Myf5", "Ptf1a", "Mdx3"}
# The alternatives are the cross of independent hypotheses, defined before release to scoring.
# The unverified Mxd3 interpretation is explicitly labeled as a sensitivity.
VARIANTS = {
    "main": ("", False, False),
    "figure_literal": ("Ascl3", False, False),
    "s1_thesis": ("Ascl4", True, False),
    "s1_named_axes": ("Ascl4", False, False),
    "figure_ascl3_swapped_axes": ("Ascl3", True, False),
    "figure_literal_mxd3_hypothesis": ("Ascl3", False, True),
    "s1_thesis_mxd3_hypothesis": ("Ascl4", True, True),
    "s1_named_axes_mxd3_hypothesis": ("Ascl4", False, True),
    "figure_ascl3_swapped_axes_mxd3_hypothesis": ("Ascl3", True, True),
}

def tsv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def occupancy(image, x_edges, y_edges):
    source = np.asarray(Image.open(image).convert("RGB"))
    result = {}
    for r, partner in enumerate(FIG_ROWS):
        for c, backbone in enumerate(FIG_COLUMNS):
            x0, x1 = x_edges[c:c+2]
            y0, y1 = y_edges[r:r+2]
            # Sample the four corners away from grid borders and printed digits.
            pixels = []
            for a,b in ((0.08,0.23),(0.77,0.92)):
                for d,e in ((0.11,0.29),(0.71,0.89)):
                    pixels.extend(source[round(y0+(y1-y0)*d):round(y0+(y1-y0)*e),
                                         round(x0+(x1-x0)*a):round(x0+(x1-x0)*b)].reshape(-1,3))
            red, green, blue = np.median(pixels, axis=0)
            call = "1" if green-red > 5 and green-blue > 7 else "0" if abs(green-red) <= 4 and abs(green-blue) <= 4 else ""
            result[backbone, partner] = (call, [round(x0),round(y0),round(x1),round(y1)], int(red),int(green),int(blue))
    return result

def mouse_genes(symbols):
    with gzip.open(PRIVATE / "mouse_gene_info.gz", "rt", encoding="utf-8", newline="") as stream:
        selected = {r["Symbol"].casefold(): r for r in csv.DictReader(stream, delimiter="\t")
                    if r["Symbol"].casefold() in symbols}
    if missing := symbols - selected.keys():
        raise ValueError(f"Official GeneInfo symbol unresolved: {sorted(missing)}")
    return selected

def construct(backbone, base_symbol, variant):
    ascl, swapped, mxd = VARIANTS[variant]
    if backbone == "Ascl3":
        return ascl
    if backbone == "Mdx3":
        return "Mxd3" if mxd else ""
    if swapped and backbone in {"Myf5", "Ptf1a"}:
        return {"Myf5":"Ptf1a", "Ptf1a":"Myf5"}[backbone]
    return base_symbol

def proof():
    # The NIH author manuscript provides aligned methods/caption and a separately rendered bitmap.
    root = ET.parse(PRIVATE / "pmc_efetch.xml")
    fig = next(f for f in root.iter("fig") if f.get("id")=="F1")
    caption = " ".join("".join(fig.find("caption").itertext()).split())
    assert "30 bHLH factors" in caption and "522" in caption
    assert (PRIVATE / "fig1.jpg").read_bytes()[:3] == b"\xff\xd8\xff"
    assert (PRIVATE / "nihms_fig1.jpg").read_bytes()[:3] == b"\xff\xd8\xff"
    return caption

def main():
    caption = proof()
    ws = load_workbook(PRIVATE / "MOESM3_ESM.xlsx", read_only=False, data_only=True).active
    roster = {}
    for row in ws.iter_rows(min_row=7,max_row=18,min_col=2,max_col=5):
        for cell in row:
            if cell.value:
                roster.setdefault(str(cell.value).strip(),[]).append(cell.coordinate)
    assert len(roster)==45 and roster["Usf2"]==["D10","C18"] and ws["O22"].value=="Atoh1"
    roster["Atoh1"]=["O22"]
    visible_roster = set(FIG_SYMBOLS) - {"Ascl3"} | {"Ascl4"}
    omitted = []
    roster_to_backbone = {}
    for literal in roster:
        sym = ALIAS.get(literal, literal)
        visible_matches = [v for v in visible_roster if sym.casefold()==v.casefold()]
        if visible_matches:
            biological = visible_matches[0]
            axis = "Ascl3" if biological=="Ascl4" else FIG_COLUMNS[FIG_SYMBOLS.index(biological)]
        else:
            axis = literal
            omitted.append((axis, sym, literal))
        assert axis not in roster_to_backbone
        roster_to_backbone[axis] = literal

    backbone_axes = [(name, symbol, roster_to_backbone[name]) for name, symbol in zip(FIG_COLUMNS, FIG_SYMBOLS)] + omitted
    assert len(backbone_axes)==46 and len(omitted)==30
    assert len(PARTNERS)==13
    symbol_set = set(PARTNERS.values()) | {construct(a,s,v) for a,s,_ in backbone_axes for v in VARIANTS}
    genes = mouse_genes({s.casefold() for s in symbol_set if s})
    def get(g):
        return genes.get(g.casefold(),{})
    def id_(g):
        return ("GeneID:" + get(g)["GeneID"]) if g else ""
    def symbol(g):
        return get(g).get("Symbol","")

    metadata = []
    for axis, default, literal in backbone_axes:
        original = (axis if axis in FIG_COLUMNS else "")
        for version in VARIANTS:
            g = construct(axis, default, version)
            candidate = get(g)
            construct_species = "Homo sapiens" if literal=="Hen2" else "Mus musculus"
            metadata.append(dict(variant=version, construct_slot=axis, role="bhlh", s1_literal=literal,
                    s1_location=",".join(roster[literal]), fig1b_literal=original,
                    discordance=axis in CONFLICT_AXES, construct_species=construct_species, recipient_species="Mus musculus",
                    mouse_symbol=symbol(g), mouse_id=id_(g), db_xrefs=candidate.get("dbXrefs",""),
                    identity_status="hypothesis" if axis in CONFLICT_AXES and g and version!="main" else
                        "unresolved" if not g or axis in CONFLICT_AXES else "resolved",
                    identity_evidence="NCBI GeneInfo official symbol and S1 roster; Fig1b header evidence where shown",
                    orthology_location="human NHLH2 to mouse Nhlh2; author methods plus NCBI orthology" if literal=="Hen2" else "",
                    notes="Exact clone accession not recovered; role/class used to disambiguate Hen1 methyltransferase alias" if literal=="Hen1" else ""))
    for axis,g in PARTNERS.items():
        for version in VARIANTS:
            metadata.append(dict(variant=version, construct_slot=axis, role="partner",
                    s1_literal=ws[PARTNER_ROSTER_CELLS[axis]].value, s1_location=PARTNER_ROSTER_CELLS[axis],
                    fig1b_literal=axis, discordance=False,
                    construct_species="Homo sapiens" if axis=="Brn3a" else "Mus musculus", recipient_species="Mus musculus",
                    mouse_symbol=symbol(g), mouse_id=id_(g), db_xrefs=get(g)["dbXrefs"],identity_status="resolved",
                    identity_evidence="S1 roster gene in parentheses, matrix row and Fig.1b",
                    orthology_location="human POU4F1 to mouse Pou4f1; author methods plus NCBI orthology" if axis=="Brn3a" else "", notes=""))
    tsv(PRIVATE / "constructs.tsv", metadata)

    private_pairs, public_pairs = [], []
    raw, figure_out, cell_audit, discrepancies = [], [], [], []
    full = occupancy(PRIVATE / "fig1.jpg", [66 + i*(1047-66)/16 for i in range(17)], [293 + i*(755-293)/13 for i in range(14)])
    nihms = occupancy(PRIVATE / "nihms_fig1.jpg", [37 + i*(598-37)/16 for i in range(17)], [167 + i*(432-167)/13 for i in range(14)])
    words = [r.split() for r in FIG_VALUES]
    assert all(len(r)==16 for r in words)
    for ib, (axis, g, _) in enumerate(backbone_axes):
        for ip, (partner, p) in enumerate(PARTNERS.items()):
            slot = f"B{ib+1:02d}-P{ip+1:02d}"
            opaque = hashlib.sha256(("MEF-neuron-design-key-v1|"+slot).encode()).hexdigest()
            main_ok = axis not in CONFLICT_AXES
            for version, (_, swapped, _) in VARIANTS.items():
                gene = construct(axis, g, version)
                row = dict(variant=version, pair_id=opaque, bhlh_id=id_(gene), partner_id=id_(p), bhlh_symbol=symbol(gene), partner_symbol=symbol(p))
                public_pairs.append(row)
                reason = "" if bool(gene) and (main_ok or version!="main") else ("disputed_axis_"+axis if not main_ok else "unmapped")
                private_pairs.append(dict(**row, pair_slot_id=slot, bhlh_axis=axis, partner_axis=partner,
                    main_eligible=int(main_ok),eligible=int(not reason),exclusion_reason=reason,
                    raw_s1_column_if_mapped=get_column_letter((33-ib) if swapped and ib in (14,15) else ib+4)
                        if ib<16 else ""))
            if ib<16:
                expected = words[ip][ib]
                digit_call = "0" if expected=="-" else "1"
                call_a,bbox,*rgb_a = full[axis,partner]
                call_b,bbox_b,*rgb_b = nihms[axis,partner]
                agrees = call_a == call_b == digit_call
                if not agrees:
                    discrepancies.append((axis,partner,expected,call_a,call_b,rgb_a,rgb_b))
                call = digit_call if agrees else ""
                basis = ("positive_fig1b_green_and_printed" if call=="1" else "negative_fig1b_blank" if call=="0" else "unresolved_color")
                val = expected if call=="1" else ""
                s1 = ws.cell(RAW_ROWS[partner],ib+4)
                s1_alt = ws.cell(RAW_ROWS[partner],(33-ib) if ib in (14,15) else ib+4)
                assert s1.value is not None and s1_alt.value is not None
                status = "printed_corrected_and_raw_s1" if call=="1" else "empty_figure_with_explicit_raw_zero" if s1.value==0 else "empty_figure_with_nonzero_raw_s1"
                src = "Fig1b row=" + partner + " col=" + axis + " bbox=" + ",".join(map(str,bbox))
                raw.append(dict(pair_slot_id=slot, source_version="nature2018_s1_raw", bhlh_axis_s1=ws.cell(22,ib+4).value,
                    partner_axis_s1=ws.cell(RAW_ROWS[partner],3).value, s1_cell=s1.coordinate, raw_s1_pct=s1.value,
                    measurement_status="explicit_raw_zero" if s1.value==0 else "positive_raw_readout",
                    denominator="initial plated MEFs",day=14))
                cell_audit.append(dict(pair_slot_id=slot, figure_row=partner, figure_column=axis,
                    springer_bbox=",".join(map(str,bbox)), nihms_bbox=",".join(map(str,bbox_b)),
                    manual_printed_value=expected, manual_call=digit_call,
                    springer_color_call=call_a, nihms_color_call=call_b, springer_corner_rgb=",".join(map(str,rgb_a)),
                    nihms_corner_rgb=",".join(map(str,rgb_b)), binary_calls_agree=agrees))
                raw_v, raw_alt = s1.value, s1_alt.value
                raw_c, raw_alt_c = s1.coordinate, s1_alt.coordinate
            else:
                call,val,raw_v,raw_alt,raw_c,raw_alt_c = "0","","","","",""
                status,basis = "not_shown_raw_unknown","group_negative_figure_caption"
                src = f"Fig1b caption 30 omitted bHLH; S1 roster {','.join(roster[roster_to_backbone[axis]])}"
            figure_out.append(dict(pair_slot_id=slot,source_version="nature2018",binary_call=call,
                q_gated=val if call=="1" else "0" if call=="0" else "",
                raw_s1_pct=raw_v,raw_s1_pct_s1_thesis=raw_alt,fig1b_pct_printed=val,
                call_basis=basis,measurement_status=status,denominator="initial plated MEFs",day=14,
                source_cell_or_bbox=src,raw_s1_cell=raw_c,raw_s1_cell_s1_thesis=raw_alt_c,
                notes="Identity ambiguity is kept in design_pairs. Negative zero is gate-defined." if not main_ok else ""))

    tsv(PRIVATE / "design_pairs.tsv", private_pairs)
    tsv(PUBLIC / "CLEAN_PAIRS.tsv", sorted(public_pairs,key=lambda r:(r["pair_id"],r["variant"])))
    tsv(PRIVATE / "figure_grid_audit.tsv", cell_audit)
    tsv(PRIVATE / "s1_raw_outcomes.tsv", raw)
    tsv(PRIVATE / "outcomes_private.tsv", figure_out)
    controls = []
    for row in list(RAW_ROWS.values()) + [31,37]:
        for col in ([20] if row not in (31,37) else range(4,21)):
            c = ws.cell(row,col)
            controls.append(dict(cell=c.coordinate, row_label=ws.cell(row,3).value, column_label=ws.cell(22,col).value, raw_s1_pct=c.value,
                role="partner-alone" if col==20 and row not in (31,37) else "control-background-row", included_as_pair=0))
    tsv(PRIVATE / "s1_control_cells.tsv", controls)
    assert len(figure_out)==598 and len(raw)==208 and len(controls)==47
    assert Counter(r["measurement_status"] for r in raw)=={"explicit_raw_zero":129,"positive_raw_readout":79}
    counts = Counter(r["binary_call"] or "unknown" for r in figure_out)
    report = {
        "design_slots": len(figure_out), "main_eligible_slots": sum(r["eligible"] for r in private_pairs if r["variant"]=="main"),
        "binary_fig": counts, "raw_measurement": Counter(r["measurement_status"] for r in raw),
        "variant_count": len(VARIANTS), "colors_agree_manual": len(discrepancies)==0,
        "color_discrepancies": discrepancies, "figure_caption": caption,
        "correction": "Printed Nature values retained; no repeated control correction.",
        "precisions": sorted(set(len(v.split(".")[1]) for row in words for v in row if "." in v)),
    }
    (PRIVATE / "OUTCOME_AUDIT.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
