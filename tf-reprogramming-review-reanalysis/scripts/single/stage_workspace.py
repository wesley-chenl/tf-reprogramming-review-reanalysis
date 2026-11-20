"""Stage a new single-TF score workspace or attach references after its score lock."""
import argparse
import csv
import gzip
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path
from zipfile import ZipFile


REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data/single"
CONFIG = REPO / "config/single"
CODE = REPO / "scripts/single"
TARGETS = ("liver", "cardiomyocyte", "betapanc")
GROUPS = {"source", *TARGETS}
R2_CFG = "FROZEN_CONFIG_ROUND2.yaml"


def tsv(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def put_tsv(path, rows):
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def jread(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def jwrite(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def fingerprint(path):
    return {"path": str(path), "bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns, "sha256": digest(path)}


def install_workspace(work):
    work.mkdir(parents=True, exist_ok=True)
    if any(work.iterdir()):
        raise ValueError("Initialize a new or empty scoring workspace")
    base, round2 = work / "base", work / "round2"
    shutil.copytree(CODE / "base", base / "code")
    shutil.copytree(CODE / "round2", round2 / "scripts")
    for root in (base, round2):
        for name in ("config", "logs", "manifests", "inputs/annotation", "environments", "validation",
                     "tmp/posix_work", "checkpoints", "sources"):
            (root / name).mkdir(parents=True, exist_ok=True)
    for name in ("FROZEN_CONFIG.yaml", "FROZEN_CONFIG_C.json"):
        shutil.copyfile(CONFIG / "base" / name, base / "config" / name)
    shutil.copyfile(CONFIG / "base/C_NATIVE_INTERFACE_POLICY.json", base / "config/C_NATIVE_INTERFACE_POLICY.json")
    for name in ("VERSIONS.json", "PREPROCESSING_TOOL_VERSIONS.json"):
        shutil.copyfile(CONFIG / "base" / name, base / "environments" / name)
    for artifact in ("CANDIDATE_UNIVERSE.tsv",):
        shutil.copyfile(DATA / artifact, base / "config" / artifact)
        (round2 / "inputs_manifest").mkdir(exist_ok=True)
        shutil.copyfile(DATA / artifact, round2 / "inputs_manifest" / artifact)
    for packed in (DATA / "raw_support").rglob("*.gz"):
        output = base / "inputs" / packed.relative_to(DATA / "raw_support").with_suffix("")
        output.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(packed, "rb") as original, output.open("wb") as destination:
            shutil.copyfileobj(original, destination)
    selected = tsv(DATA / "RAW_SOURCE_POINTERS.tsv")
    pointers = {Path(r["path"].replace(chr(92), "/")).name: r for r in selected}
    downloads = {}
    for assay in ("RNA", "ATAC"):
        rows = []
        for r in tsv(CONFIG / "base" / f"{assay}_DOWNLOAD_PLAN.tsv"):
            if r["group"] not in GROUPS:
                continue
            direct = r["group"] in {"source", "liver"}
            file = r["filename"] if direct else f'{r["run"]}.sra'
            folder = f"raw_{assay.lower()}" + ("" if direct else "_sra")
            if folder == "raw_rna_sra":
                folder = "raw_sra"
            raw = base / "inputs" / folder / file
            raw.parent.mkdir(parents=True, exist_ok=True)
            original = pointers[file]
            r["local_path"], r["reuse_cached"] = str(raw), "False"
            rows.append(r)
            stamp_dir = "file_hashes" if direct else ("sra_cloud" if assay == "RNA" else "atac_sra_cloud")
            stamp = base / "manifests" / stamp_dir / f'{file if direct else r["run"]}.json'
            downloads[str(raw)] = {"assay": assay, "sample": original["sample_id"], "path": str(raw),
                                   "url": original["source_url"], "sha256": original["sha256"],
                                   "bytes": original["bytes"], "md5": original["md5"], "stamp": str(stamp)}
        put_tsv(base / "manifests" / f"{assay}_DOWNLOAD_PLAN.tsv", rows)
        samples = [r for r in tsv(CONFIG / "base" / f"{assay}_SAMPLES_PREPROCESSING.tsv") if r["group"] in GROUPS]
        put_tsv(base / "manifests" / f"{assay}_SAMPLES_PREPROCESSING.tsv", samples)
    for resource in tsv(CONFIG / "base/PUBLIC_RESOURCE_DOWNLOADS.tsv"):
        file = base / resource["path"]
        downloads[str(file)] = {**resource, "assay": "reference", "sample": "", "path": str(file)}
    releases = jread(CONFIG / "SOURCE_RELEASES.json")
    for method in ("ANANSE", "CompForce", "IMCell"):
        release = releases[method]["archive"]
        file = round2 / "sources" / Path(release["path"]).name
        downloads[str(file)] = {"assay": "software", "sample": method, "path": str(file),
                                "url": release["URL"], "sha256": release["sha256"],
                                "bytes": release["bytes"], "md5": release["md5"] or "", "stamp": ""}
    put_tsv(work / "DOWNLOADS.tsv", downloads.values())
    run_config = jread(CONFIG / R2_CFG)
    run_config["paths"]["previous_frozen_inputs"] = str(base)
    jwrite(round2 / R2_CFG, run_config)
    shutil.copyfile(CONFIG / "SOURCE_RELEASES.json", round2 / "config/SOURCE_RELEASES.json")
    print(f"Prepared score-only workspace at {work}. Review DOWNLOADS.tsv and obtain the specified software.")


def fetch(work, select):
    for row in tsv(work / "DOWNLOADS.tsv"):
        if select != "all" and row["assay"] != select:
            continue
        path = Path(row["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        if not path.is_file():
            temporary = path.with_suffix(path.suffix + ".partial")
            request = urllib.request.Request(row["url"], headers={"User-Agent": "TF reprogramming review reproducibility"})
            with urllib.request.urlopen(request, timeout=120) as stream, temporary.open("wb") as output:
                shutil.copyfileobj(stream, output, length=8 * 1024**2)
            if digest(temporary) != row["sha256"] or temporary.stat().st_size != int(row["bytes"]):
                raise ValueError(f"Public source changed: {row['url']}")
            temporary.replace(path)
        rec = fingerprint(path)
        if rec["sha256"] != row["sha256"] or rec["bytes"] != int(row["bytes"]):
            raise ValueError(f"Source checksum mismatch: {path}")
        if row["stamp"]:
            jwrite(Path(row["stamp"]), {**rec, "md5": row["md5"], "status": "verified", "url": row["url"]})
        if row["assay"] == "software":
            unpack_release(path, path.parent / row["sample"])
        print(f"Verified {row['assay']} {path.name}", flush=True)


def unpack_release(archive, dest):
    if dest.is_dir():
        return
    temporary = dest.with_suffix(".partial")
    temporary.mkdir(parents=True)
    with ZipFile(archive) as package:
        for member in package.infolist():
            parts = Path(member.filename).parts
            if ".." in parts or Path(member.filename).is_absolute():
                raise ValueError(f"Invalid source archive member: {member.filename}")
            if member.is_dir() or len(parts) < 2:
                continue
            output = temporary.joinpath(*parts[1:])
            output.parent.mkdir(parents=True, exist_ok=True)
            with package.open(member) as stream, output.open("wb") as handle:
                shutil.copyfileobj(stream, handle)
    temporary.replace(dest)


def freeze_base(work):
    base = work / "base"
    all_inputs = {}
    for target in TARGETS:
        folder = base / "inputs/rna" / target
        names = ("TPM.tsv", "TARGET_TPM.tsv", "DESEQ2_TARGET_V_SOURCE.tsv", "SAMPLE_METADATA.tsv", "ROUNDED_COUNTS.tsv")
        task = {str(folder / name): fingerprint(folder / name) for name in names}
        all_inputs.update(task)
        jwrite(base / "config" / f"TASK_INPUT_FREEZE_{target}.json", {"task": target, "actual_input_hashes": task})
    for relative in ("config/CANDIDATE_UNIVERSE.tsv", "inputs/priors/TRIAGE_RTS_M24.tsv", "inputs/priors/COLLECTRI_MOUSE_M24.tsv"):
        f = base / relative
        all_inputs[str(f)] = fingerprint(f)
    jwrite(base / "config/FROZEN_INPUT_HASHES.json", all_inputs)
    jwrite(base / "PRE_EVAL_FREEZE.json", {"inputs": list(all_inputs.values()), "reference_files_available": False})
    gtf = base / "inputs/annotation/gencode.vM24.annotation.gtf"
    jwrite(base / "manifests/C_NATIVE_ANNOTATION_STAGE.json", {"status": "verified", "native_logical_path": str(gtf)})
    peak = base / "inputs/atac_conditions/source/source_pool_peaks.narrowPeak"
    atac = list((base / "inputs/atac").glob("*/PREPROCESSING_COMPLETE.json"))
    if len(atac) != 5:
        raise ValueError(f"Expected five complete biological ATAC inputs, found {len(atac)}")
    index = jread(base / "runs/preprocessing/TAIJI_GENOME_INDEX/RUN_RECORD.json")
    jwrite(base / "PRE_EVAL_FREEZE_C.json", {"peaks": fingerprint(peak), "atac_records": [fingerprint(f) for f in atac],
                                           "genome_index": index["native_index_sha256"],
                                           "reference_files_available": False})


def prepare_round2(work):
    base, round2 = work / "base", work / "round2"
    assembly = base / "inputs/annotation/GRCm38.p6.genome.fa"
    records = sorted((base / "inputs/atac").glob("*/PREPROCESSING_COMPLETE.json"))
    if not assembly.with_suffix(".fa.fai").is_file() or len(records) != 5:
        raise ValueError("Index the GRCm38 FASTA and complete all five biological ATAC libraries first")
    bams = []
    for record in records:
        raw = jread(record)
        bams.append({k: raw[k] for k in ("sample_id", "layout", "bam", "peaks")})
    put_tsv(round2 / "inputs_manifest/ATAC_INPUT_CHECKS.tsv", bams)
    jwrite(round2 / "inputs_manifest/VERIFIED_INPUTS.json", {"assembly": str(assembly), "reference_label_or_evaluation_accessed": False})
    cfg_hash = digest(round2 / R2_CFG)
    policy = jread(CONFIG / "ANANSE_PARTITIONED_NETWORK_POLICY.json")
    policy.update(config_sha256=cfg_hash, completed_ranked_binding_checkpoints={}, completed_chromosome_checkpoints={})
    jwrite(round2 / "config/ANANSE_PARTITIONED_NETWORK_POLICY.json", policy)
    adaptation = jread(CONFIG / "IMCELL_RUNTIME_ADAPTATION.json")
    if digest(round2 / "scripts/imcell_indexed_cpu.R") != adaptation["indexed_code_sha256"]:
        raise ValueError("IMCell adapter differs from the version that passed parity")
    jwrite(round2 / "checkpoints/IMCELL_RUNTIME_READY.json",
           {"config_sha256": cfg_hash, "indexed_code_sha256": adaptation["indexed_code_sha256"]})


def prelock(work):
    base, round2 = work / "base", work / "round2"
    lineage, source_records = [], []
    cfg_sha = digest(round2 / R2_CFG)
    for target in TARGETS:
        for method, slot in (("TRIAGE", "A"), ("VIPER", "B"), ("Taiji", "C")):
            folder = base / "runs" / target / slot
            record = jread(folder / "RUN_RECORD.json")
            score = fingerprint(folder / "FULL_U_SCORES.tsv")
            if record["status"] != "REAL_SCORING_COMPLETE" or record["output_sha256"] != score["sha256"] or record["evaluation_used"]:
                raise ValueError(f"Base method did not finish blindly: {target}/{method}")
            lineage.append({"target": target, "method": method, "original_score_file": score["path"],
                            "original_score_sha256": score["sha256"], "previous_config_sha256": record["config_hash"],
                            "status": "reuse_verified"})
            source_records.append(score)
        for native in ("ANANSE/influence_RUN.json", "CompForce/METHOD_RUN.json", "IMCell/METHOD_RUN.json"):
            f = round2 / "results" / target / native
            rec = jread(f)
            if rec["execution_status"] not in {"completed", "completed_empty", "completed_uninformative"} or rec["config_sha256"] != cfg_sha:
                raise ValueError(f"Native stage is incomplete or uses different configuration: {f}")
            source_records.append(fingerprint(f))
    put_tsv(round2 / "inputs_manifest/REUSE_LINEAGE.tsv", lineage)
    put_tsv(round2 / "SOURCE_MANIFEST.tsv", source_records)
    jwrite(round2 / "STATUS.json", {"dispatcher_complete": True, "running": [], "config_sha256": cfg_sha})


def attach_references(work):
    base, round2 = work / "base", work / "round2"
    lock = jread(round2 / "PREDICTIONS_LOCK_ROUND2.json")
    if lock["reference_files_opened"] or not lock["workflow_scoring_terminal"]:
        raise ValueError("A complete blind lock is required")
    for row in lock["outputs"] + lock["lineage_sources"]:
        path = Path(row["path"])
        if not path.resolve().is_relative_to(work) or digest(path) != row["sha256"]:
            raise ValueError(f"Blind lock checksum mismatch: {path}")
    refs = {
        "REFERENCE_TFS.tsv": "REFERENCE_TFS_FROZEN_ROUND2.tsv",
        "REFERENCE_EVIDENCE.tsv": "REFERENCE_EVIDENCE_FROZEN_ROUND2.tsv",
        "MOTIF_MAPPING.tsv": "MOTIF_MAPPING_FROZEN.tsv",
        "REFERENCE_EVIDENCE_REPORTING_VIEW.tsv": "REFERENCE_EVIDENCE_REPORTING_ROUND2.tsv",
    }
    records = []
    for name, src in refs.items():
        path = base / "references" / name
        path.parent.mkdir(exist_ok=True)
        shutil.copyfile(DATA / src, path)
        records.append({**fingerprint(path), "role": "evaluation_only"})
    put_tsv(base / "FREEZE_MANIFEST.tsv", records)
    put_tsv(base / "OUTPUT_MANIFEST.tsv", records)
    print("Blind lock verified. Historical reference tables staged for evaluation.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("init", "download", "freeze-base", "prepare-round2", "prelock", "attach-references"))
    parser.add_argument("work", type=Path)
    parser.add_argument("--assay", choices=("RNA", "ATAC", "reference", "software", "all"), default="all")
    args = parser.parse_args()
    work = args.work.resolve()
    if args.stage == "download":
        fetch(work, args.assay)
    else:
        action = {"init": install_workspace, "freeze-base": freeze_base, "prepare-round2": prepare_round2,
                  "prelock": prelock, "attach-references": attach_references}[args.stage]
        action(work)


if __name__ == "__main__":
    main()
