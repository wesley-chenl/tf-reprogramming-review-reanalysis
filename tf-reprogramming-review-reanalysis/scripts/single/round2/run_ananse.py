"""One restartable ANANSE state/task stage, with output-schema and score gates."""
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

from r2lib import ANA, ATAC_SRC, ATAC_TGT, BASE, ROOT, RUNTIME, SCRATCH, input_record, now
from r2lib import read_json, read_tsv, run_env, sha, write_json, write_tsv


def validate_binding(output, regions):
    expected = {f"{r[0]}:{r[1]}-{r[2]}" for line in regions.read_text().splitlines() if (r := line.split("\t"))}
    with pd.HDFStore(output, "r") as hdf:
        keys = set(hdf.keys())
        assert {"/_atac", "/_factor_activity", "/_index"} <= keys, keys
        indices = list(hdf["_index"].astype(str))
        assert len(indices) == len(expected) and set(indices) == expected, (len(indices), len(expected))
        factors = sorted(k[1:] for k in keys if not k.startswith("/_"))
        assert len(factors) > 0 and len(set(factors)) == len(factors)
        activity = hdf["_factor_activity"]
        assert set(activity.columns)=={"factor","activity"} and activity["factor"].is_unique
        assert set(factors)<=set(activity["factor"]) and np.isfinite(activity["activity"].astype(float)).all()
        raw_atac = hdf["_atac"]
        numeric = np.asarray(raw_atac.select_dtypes(include="number"), dtype=float)
        assert numeric.size and np.isfinite(numeric).all()
        for tf in factors:
            p = np.asarray(hdf[tf], dtype=float)
            assert len(p) == len(indices) and np.isfinite(p).all() and np.min(p) >= 0 and np.max(p) <= 1, tf
    return {"region_n": len(indices), "finite_binding_factors": factors,"finite_activity_factor_n":len(activity),
            "native_HDF5_keys": sorted(keys), "ATAC_numeric_columns": list(raw_atac.select_dtypes(include="number").columns)}


def validate_network(output, stats_file, *, expression, binding_record):
    edges = 0
    outdegree = Counter()
    nodes = set()
    expressions = read_tsv(expression)
    genes = {r["gene"]:i for i,r in enumerate(expressions)}
    assert len(genes)==len(expressions)
    regulators = set(binding_record["finite_binding_factors"]) & genes.keys()
    # Dense native expression product is complete when every possible pair appears once.
    # A bit per edge validates identity without retaining tens of millions of string keys.
    seen = {tf:bytearray((len(genes)+7)//8) for tf in regulators}
    bits = (1,2,4,8,16,32,64,128)
    with output.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream, delimiter="\t")
        names = next(reader)
        assert {"tf_target", "prob"} <= set(names), names
        i_edge, i_p = names.index("tf_target"), names.index("prob")
        for row in reader:
            tf, dest = row[i_edge].split("—", 1)
            value = float(row[i_p])
            assert tf and dest and math.isfinite(value) and 0 <= value <= 1, row
            i = genes[dest]
            ledger = seen[tf]
            byte = i>>3
            flag = bits[i&7]
            assert not ledger[byte]&flag, f"Duplicate native edge {tf}→{dest}"
            ledger[byte] |= flag
            outdegree[tf] += 1
            nodes.update((tf, dest))
            edges += 1
    assert edges==len(genes)*len(regulators)>0 and set(outdegree)==regulators and len(nodes)==len(genes), (edges,len(genes),len(regulators))
    write_tsv(stats_file, ({"symbol": tf, "out_degree": count} for tf, count in sorted(outdegree.items())))
    return {"network_edges": edges, "graph_nodes": len(nodes), "regulators": len(outdegree), "schema": names,
            "full_expected_unique_pairs_verified":True,"expected_expression_genes":len(genes),"out_degree_file": str(stats_file)}


def native_log_anomalies(log, stage):
    full = log.read_text(encoding="utf-8",errors="replace")
    tag = "Failed to communicate with scheduler during heartbeat."
    tolerated = ""
    if stage=="network" and tag in full and "Writing network" in full:
        before, after = full.split(tag,1)
        # An isolated, terminal heartbeat traceback on closure can be accepted only with a verified full output.
        endings = after.strip().splitlines()
        if ("Writing network" in before and "An error has been caught" not in after and endings and
                endings[-1].startswith("distributed.comm.core.CommClosedError:") and
                "heartbeat_worker" in endings[-1] and "Stream is closed" in endings[-1]):
            final_error = full.rfind("\n",0,full.index(tag))
            full = full[:final_error]
            tolerated = "terminal Dask heartbeat disconnect at scheduler close; exact full network cardinality and unique pairs required"
    errors = [line.strip() for line in full.splitlines() if re.search(r"Traceback \(most recent call last\)|An error has been caught|\b(CRITICAL|ERROR)\s*\|",line)]
    return errors,tolerated


def validate_influence(output, gatefile):
    rows = read_tsv(output)
    assert rows and {"factor", "influence_score"} <= set(rows[0]), rows[:1]
    names = [r["factor"] for r in rows]
    assert len(names) == len(set(names)) and all(n.strip() for n in names), "Duplicate or empty native influence identity"
    scores = [float(r["influence_score"]) for r in rows]
    tol = 8 * sys.float_info.epsilon
    assert all(math.isfinite(s) and -tol <= s <= 1 + tol for s in scores), (min(scores),max(scores))
    gate = {r["symbol"]: r for r in read_tsv(gatefile)}
    allowed = {s for s, r in gate.items() if int(r["differential_network_out_degree"]) > 0 and
               r["native_DE_present"] == "True" and float(r["native_real_fold_change"]) > 0}
    assert set(names) == allowed, (set(names) - allowed, allowed - set(names))
    return {"influence_factors": len(names), "finite_scores": len(scores), "schema": list(rows[0]),
            "native_gate_factors": len(allowed), "native_gate_genes": len(gate),
            "range_validation_FP_epsilon":tol,"out_of_exact_unit_interval_FP_n":sum(s<0 or s>1 for s in scores),
            "native_min_score":min(scores),"native_max_score":max(scores),"native_score_clipped":False}


def publish_exact(local, durable):
    if local == durable:
        return None
    durable.parent.mkdir(exist_ok=True)
    pending = durable.with_name(durable.name + ".partial")
    shutil.copyfile(local,pending)
    digest = sha(local)
    assert digest==sha(pending)
    pending.replace(durable)
    return digest


def main(override_threads=None, override_memory=None):
    task, stage, state, threads = sys.argv[1:5]
    if override_threads is not None:
        threads = str(override_threads)
    config = ROOT / "FROZEN_CONFIG_ROUND2.yaml"
    config_sha = sha(config)
    dask_policy = ROOT / "config/ANANSE_DASK_NETWORK_POLICY.json"
    if stage == "network" and dask_policy.is_file() and override_threads is None:
        strategy = read_json(dask_policy)
        assert strategy["config_sha256"] == config_sha
        threads = str(strategy["Dask_workers"])
        override_memory = strategy["memory_limit_per_worker"]
    native = ROOT / "scripts/ananse_native.py"
    working_inputs = ROOT / "inputs_manifest/ananse"
    manifest = read_json(ROOT / "inputs_manifest/VERIFIED_INPUTS.json")
    fasta = manifest["assembly"]
    staging_file = ROOT / "config/ANANSE_IO_STAGING.json"
    staged = read_json(staging_file)["copies"] if staging_file.is_file() else {}

    def identical_disk_path(original):
        if original not in staged:
            return original
        record = staged[original]
        target = Path(record["path"])
        original_record = manifest["files"].get(original) or manifest["files"][str(Path(original).resolve())]
        assert original_record["sha256"] == record["sha256"]
        if not target.is_file():
            return original
        assert target.stat().st_size == record["bytes"]
        return str(target)

    fasta = identical_disk_path(fasta)
    local = working_inputs / task
    data = ROOT / "results" / task / "ANANSE"
    data.mkdir(parents=True, exist_ok=True)
    cpu_policy = ROOT / "config/ANANSE_LOCAL_CPU_ALLOCATION.json"
    if stage == "binding" and cpu_policy.is_file():
        strategy = read_json(cpu_policy)
        assert strategy["config_sha256"] == config_sha
        # Increase native scanner slots only after all other native method jobs are terminal.
        independent = [ROOT / "results" / t / m / "METHOD_RUN.json" for t in ("liver", "cardiomyocyte", "betapanc")
                       for m in ("CompForce", "IMCell")]
        if all(p.is_file() and read_json(p)["execution_status"] in {"completed", "completed_empty", "completed_uninformative"} for p in independent):
            threads = str(strategy["binding_native_n"])
    regions = local / "regions_combined.bed"
    label = f"{stage}_{state}" if stage != "influence" else stage
    result_record = data / f"{label}_RUN.json"
    if result_record.is_file():
        prior = read_json(result_record)
        if prior.get("execution_status") == "completed" and prior.get("config_sha256") == config_sha:
            for record in prior["outputs"]:
                assert sha(record["path"]) == record["sha256"], record["path"]
            print("Validated checkpoint", task, label)
            return prior
    attempt = 1 + len(list(data.glob(f"{label}_attempt*.log")))
    work_volume = SCRATCH if staged else RUNTIME
    scratch = work_volume / "work" / "ananse" / config_sha[:12] / task / label / f"attempt_{attempt}"
    scratch.mkdir(parents=True, exist_ok=False)
    log = data / f"{label}_attempt{attempt}.log"
    all_bams = {r["sample_id"]: identical_disk_path(r["bam"]) for r in read_tsv(ROOT / "inputs_manifest/ATAC_INPUT_CHECKS.tsv")}
    env = run_env()
    env.update({"PATH": str(ANA / "bin") + ":" + env["PATH"], "TMPDIR": str(scratch),
                "PYTHONHASHSEED": "0", "R2_DASK_WORKER_MEMORY": override_memory or os.environ.get("R2_DASK_WORKER_MEMORY", "3400MiB"),
                "R2_INFLUENCE_GATE": str(data / "INFLUENCE_NATIVE_GATE.tsv"), "R2_EDGE_GATE":str(data/"INFLUENCE_NETWORK_BOUNDARY.tsv"),
                "DASK_DISTRIBUTED__DIAGNOSTICS__NVML": "False"})
    cli = [str(ANA / "bin/python"), str(native), str(ANA / "bin/ananse"), stage]
    local_output = None
    if stage == "binding":
        bams = [all_bams[s] for s in (ATAC_SRC if state == "source" else [ATAC_TGT[task]])]
        output = data / f"{state}.binding/binding.h5"
        local_output = scratch / "native.binding/binding.h5" if staged else output
        env["R2_MOTIF_CACHE_AUDIT"] = str(data / f"{label}_MOTIF_CACHE.json")
        proof = ROOT / "validation/ANANSE_MOTIF_ROW_CACHE_PARITY.json"
        cache_code_hash = sha(ROOT / "scripts/ananse_motif_cache.py")
        if proof.is_file() and read_json(proof).get("passed") and read_json(proof)["scorer_sha256"] == cache_code_hash:
            cache_identity = sha256(":".join([config_sha, sha(regions), sha(native), cache_code_hash]).encode()).hexdigest()
            env["R2_MOTIF_CACHE_PREFIX"] = str(work_volume / "work/ananse_motif_by_task" / task / cache_identity)
        cli += ["-A", *bams, "-g", fasta, "-r", str(regions), "-p", str(working_inputs / "gimme.vertebrate.v5.0.pfm"),
                "-n", threads, "-o", str(local_output.parent)]
    elif stage == "network":
        output = data / f"{state}.network.tsv"
        local_output = scratch / f"{state}.network.tsv" if staged else output
        hdf = data / f"{state}.binding/binding.h5"
        hdf_local = scratch / "binding.h5" if staged else hdf
        if hdf_local != hdf:
            shutil.copyfile(hdf,hdf_local)
            assert sha(hdf)==sha(hdf_local)
        cli += [str(hdf_local), "-e", str(local / f"{state}_mean_TPM.tsv"), "-c", "tpm", "-g", fasta,
                "-a", str(working_inputs / "M24_gene_symbols.bed"), "--include-promoter", "--include-enhancer", "-n", threads, "-o", str(local_output)]
        low_memory = ROOT/"config/ANANSE_PARTITIONED_NETWORK_POLICY.json"
        if low_memory.is_file():
            policy = read_json(low_memory)
            assert policy["config_sha256"] == config_sha
            assert policy["adapter_sha256"] == sha(ROOT/"scripts/ananse_partitioned_network.py")
            env["R2_ANANSE_PARTITIONED_NETWORK"] = "1"
            env["R2_ANANSE_CARTESIAN_REGULATORS_PER_PARTITION"] = str(policy["regulator_partition_max"])
            compact = policy.get("compact_archive_binding_rank_for_targets",[])
            if task in compact:
                assert policy["compact_rank_code_sha256"]==sha(ROOT/"scripts/ananse_compact_binding_rank.py")
                env["R2_ANANSE_COMPACT_RANK"] = "1"
                chrom_record = policy.get("completed_chromosome_checkpoints",{}).get(f"{task}/{state}")
                if chrom_record:
                    x = read_json(chrom_record)
                    assert x["config_sha256"]==config_sha and x["state_binding_sha256"]==sha(hdf)
                    assert x["region_file_sha256"]==sha(regions) and x["annotation_sha256"]==sha(working_inputs/"M24_gene_symbols.bed")
                    if all(Path(p["path"]).is_file() and sha(p["path"])==p["sha256"] for p in x["archive_tables_after_all_21_contigs"]):
                        env["R2_ANANSE_ARCHIVE_CHROMOSOME_CACHE"] = x["aggregate_chromosome_directory"]
            saved_checkpoint = policy["completed_ranked_binding_checkpoints"].get(f"{task}/{state}")
            if saved_checkpoint:
                saved = read_json(saved_checkpoint)
                cached = saved["immutable_cache"]
                assert saved["config_sha256"]==config_sha and saved["binding_h5_sha256"]==sha(hdf)
                assert saved["regions_sha256"]==sha(regions) and saved["bed_annotation_sha256"]==sha(working_inputs/"M24_gene_symbols.bed")
                if Path(cached["path"]).is_file() and sha(cached["path"])==cached["sha256"]:
                    env["R2_ANANSE_RANKED_BINDING_CSV"] = cached["path"]
    elif stage == "influence":
        output = data / "influence.tsv"
        networks = {}
        for status in ("source","target"):
            r = read_json(data/f"network_{status}_RUN.json")
            assert r["execution_status"]=="completed" and r["config_sha256"]==config_sha
            permanent = data/f"{status}.network.tsv"
            native_bytes = Path(r.get("native_writes_path",permanent))
            networks[status] = permanent
            if native_bytes!=permanent and native_bytes.is_file():
                expected = next(z["sha256"] for z in r["outputs"] if z["path"]==str(permanent))
                if sha(native_bytes)==expected:
                    networks[status] = native_bytes
        cli += ["-s", str(networks["source"]), "-t", str(networks["target"]), "-d",
                str(local / "DESEQ2_TARGET_V_SOURCE_gene_symbols.tsv"), "-a", str(BASE / "inputs/annotation/gencode.vM24.annotation.gtf.gz"),
                "-i", "500000", "-j", "0.05", "-c", "prob", "-n", threads, "-o", str(output)]
    else:
        raise ValueError(stage)
    begin = now()
    record = {"target": task, "stage": stage, "state": state, "config_sha256": config_sha,
              "wrapper_sha256": sha(Path(__file__)), "native_launcher_sha256": sha(native), "started_UTC": begin,
              "argv": cli, "run_attempt": attempt, "dask_worker_memory": env["R2_DASK_WORKER_MEMORY"],
              "native_slot_option_n": int(threads),
              "execution_status": "running", "log": str(log), "scratch": str(scratch),
              "byte_identical_WSL_staging_manifest": str(staging_file) if staged else None,
              "native_writes_path": str(local_output or output),
              "partitioned_network_export":env.get("R2_ANANSE_PARTITIONED_NETWORK") == "1",
              "compact_ranked_binding":env.get("R2_ANANSE_COMPACT_RANK") == "1",
              "reused_archived_chromosomes":env.get("R2_ANANSE_ARCHIVE_CHROMOSOME_CACHE"),
              "reused_stock_ranked_binding_path":env.get("R2_ANANSE_RANKED_BINDING_CSV")}
    write_json(result_record, record)
    try:
        with log.open("w", encoding="utf-8") as stream:
            child = subprocess.run(cli, env=env, stdout=stream, stderr=subprocess.STDOUT)
        record["return_code"] = child.returncode
        assert child.returncode == 0, f"native CLI exit {child.returncode}"
        errors,close_diagnostic = native_log_anomalies(log,stage)
        assert not errors, errors[:12]
        observed_output = local_output or output
        assert observed_output.is_file() and observed_output.stat().st_size > 0, f"Missing/nonempty native output: {observed_output}"
        outputs = [output]
        verified_digest = None
        if stage == "binding":
            detail = validate_binding(observed_output, regions)
            binding_log = log.read_text(encoding="utf-8",errors="replace")
            assert all(f"Motif activity prediction on ATAC data, run {i}/3" in binding_log for i in (1,2,3))
            assert binding_log.count("Fitting BayesianRidge")==3, "One or more native ATAC activity fits missing"
            assert len(re.findall(r"\bINFO - Done\s*$",binding_log,re.MULTILINE))>=3, "One or more native ATAC fits never reached completion"
            detail["verified_ATAC_activity_fits"]=3
            verified_digest = publish_exact(observed_output,output)
            cache_audit = Path(env["R2_MOTIF_CACHE_AUDIT"])
            if cache_audit.is_file():
                outputs.append(cache_audit)
        elif stage == "network":
            stats_file = data / f"{state}_NETWORK_OUTDEGREE.tsv"
            detail = validate_network(observed_output, stats_file,expression=local/f"{state}_mean_TPM.tsv",
                                      binding_record=read_json(data/f"binding_{state}_RUN.json"))
            verified_digest = publish_exact(observed_output,output)
            outputs.append(stats_file)
        else:
            gatefile = Path(env["R2_INFLUENCE_GATE"])
            detail = validate_influence(output, gatefile)
            outputs += [gatefile, data / "influence_diffnetwork.tsv", Path(env["R2_EDGE_GATE"])]
        record.update(detail)
        if close_diagnostic:
            record["native_shutdown_diagnostic"] = close_diagnostic
        record.update(execution_status="completed", outputs=[input_record(p, "validated native ANANSE output",
                            digest=verified_digest if p==output else None) for p in outputs])
    except Exception as exc:
        record.update(execution_status="failed_execution", reason=f"{type(exc).__name__}: {exc}")
        print(record["reason"], file=sys.stderr)
    record["ended_UTC"] = now()
    write_json(result_record, record)
    print(task, label, record["execution_status"], log)
    return record


if __name__ == "__main__":
    outcome = main()
    if outcome["stage"] == "network" and outcome["execution_status"] != "completed":
        native_log = Path(outcome["log"]).read_text(encoding="utf-8", errors="replace")
        memory_failure = re.search(r"KilledWorker|memory budget|MemoryError|worker.*memory limit|Worker process.*memory", native_log, re.IGNORECASE)
        independent_done = all((ROOT/"results"/t/m/"METHOD_RUN.json").is_file() for t in ("liver","cardiomyocyte","betapanc") for m in ("CompForce","IMCell"))
        if memory_failure and independent_done and outcome.get("native_slot_option_n",2) > 1:
            print("Native Dask memory failure. Retrying the identical network with one two-thread worker and more memory per worker.", flush=True)
            outcome = main(override_threads=1, override_memory="7350MiB")
    raise SystemExit(0 if outcome["execution_status"] == "completed" else 1)
