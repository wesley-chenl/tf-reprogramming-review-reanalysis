"""Record evidenced CPU compatibility before resuming the independent jobs."""
import difflib

from r2lib import ROOT, input_record, read_json, sg, sha, write_json


validation = ROOT / "validation"
quick = validation / "IMCELL_INDEXED_quick_PARITY.json"
full = validation / "IMCELL_INDEXED_liver_PARITY.json"
checks = [read_json(p) for p in (quick, full)]
assert all(r["pass"] for r in checks)
assert len(checks[0]["calls"]) == 15 and checks[1]["calls"][0]["niter"] == 1000
assert all(c["bit_identical_full_rank_result_and_rng"] for r in checks for c in r["calls"])
original_runner = ROOT / "scripts/imcell_run.R"
new_code = ROOT / "scripts/imcell_indexed_cpu.R"
archive_runner = ROOT / "sources/compatibility/imcell_run_at_freeze.R"
if not archive_runner.exists():
    archive_runner.write_bytes(original_runner.read_bytes())
freeze = ROOT / "FROZEN_CONFIG_ROUND2.yaml"
freeze_record = next(r for r in read_json(freeze)["versions_and_prepared_input_files"] if r["path"] == str(original_runner))
assert sha(archive_runner) == freeze_record["sha256"]
text = archive_runner.read_text()
gate_code = '''if (!preflight) {
  adapt <- fromJSON(file.path(root,"checkpoints/IMCELL_RUNTIME_READY.json"))
  code <- file.path(root,"scripts/imcell_indexed_cpu.R")
  stopifnot(adapt$config_sha256==file_sha(file.path(root,"FROZEN_CONFIG_ROUND2.yaml")),
            adapt$indexed_code_sha256==file_sha(code))
  source(code)
  r2_install_indexed_mc()
  cat("IMCell archive rank_by_spread with exact ordered integer lookups:",adapt$indexed_code_sha256,"\\n")
}

'''
anchor = 'score_job <- function(job,niter=1000L)'
assert text.count(anchor) == 1
adapted = text.replace(anchor, gate_code + anchor)
temp = original_runner.with_suffix(".partial")
temp.write_text(adapted)
patch = ROOT / "sources/compatibility/IMCell_exact_ordered_cpu.diff"
patch.write_text("".join(difflib.unified_diff(text.splitlines(True), adapted.splitlines(True),
                fromfile="sealed/imcell_run.R", tofile="working/imcell_run.R")))
revision = {"date_Singapore": sg(), "config_sha256": sha(freeze), "indexed_code_sha256": sha(new_code),
            "change": "Execution-only index lookup and integer bookkeeping inside the archived native rank_by_spread MC evaluator. Identical edge/vertex order, R uniform/permutation RNG, signs, target assignments, finished stream output, unmodified function arguments and full graph.",
            "initial_production_intervention": "17 workers exceeded memory allowance; old completed rows remain and can be reproduced from job RNG streams.",
            "proof_before_resuming": [input_record(p, "native and CPU-indexed full object/RNG equality") for p in (quick, full)],
            "further_full_parity": "Run one highest-out-degree fixed seed at niter=1000 on each remaining real program before validation completes.",
            "original_runner_sha256": sha(archive_runner), "new_runner_sha256": sha(temp), "diff_sha256": sha(patch),
            "scientific_freeze_preserved": True}
write_json(ROOT / "config/IMCELL_RUNTIME_ADAPTATION.json", revision)
temp.replace(original_runner)
write_json(ROOT / "checkpoints/IMCELL_RUNTIME_READY.json", revision)
print("Exact CPU runtime released", revision["indexed_code_sha256"], checks[1]["calls"][0])
