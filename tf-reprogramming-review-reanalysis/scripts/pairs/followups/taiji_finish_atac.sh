#!/usr/bin/env bash
set -euo pipefail
code="$(cd "$(dirname "$0")" && pwd)"
root=$(realpath "${PAIR_TAIJI_ROOT:?Set the Taiji staging root}")
assets=$(realpath "${PAIR_TAIJI_ASSETS:?Set the Taiji assets root}")
work="$root/work"
mkdir -p "$root/logs"
runsets=(SRR5822253 SRR5822288 SRR5822289 SRR5822290 SRR5822291)
# samtools and macs3 must be on PATH.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$work"
suffix=""
if [ "$#" -gt 1 ]; then suffix="_$2"; fi
exec > >(tee -a "$root/logs/TAIJI_ATAC_FINISH$suffix.log") 2>&1
set -x
callpeak() {
  local label="$1" bam="$2"
  local dest="$work/peaks/$label"
  mkdir -p "$dest"
  if [ ! -s "$dest/$label.peaks_pass" ]; then
    macs3 callpeak -t "$bam" -f BAMPE -g mm -q 0.01 --keep-dup all -n "$label" --outdir "$dest"
    test -s "$dest/${label}_peaks.narrowPeak"
    sha256sum "$dest/${label}_peaks.narrowPeak" > "$dest/$label.peaks_pass"
  fi
  sha256sum --check "$dest/$label.peaks_pass"
  cp -- "$dest/${label}_peaks.narrowPeak" "$root/taiji/receipts/${label}_peaks.narrowPeak"
}
deduplicate() {
  local label="$1"
  shift
  local paths=()
  for run in "$@"; do
    paths+=("$work/$run.names.bam")
    test -s "$root/taiji/receipts/$run.alignment.tsv"
  done
  local output_lock
  exec {output_lock}>"$work/$label.qc.lock"
  flock -n "$output_lock" || return 2
  local out="$work/$label.filtered.bam"
  if [ -s "$out" ] && [ -s "$root/taiji/receipts/$label.peak_filter.tsv" ]; then
    samtools quickcheck -v "$out"
  else
    if [ "${#paths[@]}" -gt 1 ]; then
      samtools merge -n -@ 2 -O SAM -o - "${paths[@]}"
    else
      samtools view -h "${paths[0]}"
    fi \
      | python3 "$code/taiji_atac_qc.py" filter "$root/taiji/receipts/$label.peak_filter.tsv" \
      | samtools view -u - \
      | samtools fixmate -m -u - - \
      | samtools sort -@ 3 -m 768M -T "$work/$label.coord" -u - \
      | samtools markdup -r -s -f "$root/taiji/receipts/$label.markdup.txt" - "$out.partial"
    samtools quickcheck -v "$out.partial"
    mv -- "$out.partial" "$out"
    samtools index -@ 2 "$out"
    samtools flagstat -O tsv "$out" > "$root/taiji/receipts/$label.filtered.tsv"
  fi
  callpeak "$label" "$out"
  if [ -s "$root/taiji/receipts/$label.metrics.json" ]; then exec {output_lock}>&-; return; fi
  samtools view "$out" \
    | python3 "$code/taiji_atac_qc.py" metrics \
      "$assets/TAIJI_RESOURCES/gencode.vM24.annotation.gtf" \
      "$work/peaks/$label/${label}_peaks.narrowPeak" "$root/taiji/receipts/$label"
  exec {output_lock}>&-
}
mode="${1:-all}"
case "$mode" in
  run) shift; for run in "$@"; do deduplicate "$run" "$run"; done ;;
  group)
    case "$2" in
      cortex_B1) deduplicate cortex_B1 SRR5822288 SRR5822289 ;;
      cortex_B2) deduplicate cortex_B2 SRR5822290 SRR5822291 ;;
      *) exit 2 ;;
    esac ;;
  all)
    for run in "${runsets[@]}"; do deduplicate "$run" "$run"; done
    deduplicate cortex_B1 SRR5822288 SRR5822289
    deduplicate cortex_B2 SRR5822290 SRR5822291
    target="$work/cortex_pooled.filtered.bam"
    if [ ! -s "$target" ]; then
      samtools merge -@ 4 -o "$target.partial" "$work/cortex_B1.filtered.bam" "$work/cortex_B2.filtered.bam"
      samtools quickcheck -v "$target.partial"
      mv -- "$target.partial" "$target"
      samtools index -@ 2 "$target"
    fi
    callpeak cortex_pooled "$target"
    samtools view "$target" \
      | python3 "$code/taiji_atac_qc.py" metrics \
        "$assets/TAIJI_RESOURCES/gencode.vM24.annotation.gtf" \
        "$work/peaks/cortex_pooled/cortex_pooled_peaks.narrowPeak" "$root/taiji/receipts/cortex_pooled"
    python3 "$code/taiji_atac_qc.py" summarize "$root/taiji/receipts" "$root/taiji"
    ;;
  *) exit 2 ;;
esac
