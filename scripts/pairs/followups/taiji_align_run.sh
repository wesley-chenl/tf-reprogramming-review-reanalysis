#!/usr/bin/env bash
set -euo pipefail
code="$(cd "$(dirname "$0")" && pwd)"
root=$(realpath "${PAIR_TAIJI_ROOT:?Set the Taiji staging root}")
assets=$(realpath "${PAIR_TAIJI_ASSETS:?Set the Taiji assets root}")
work="$root/work"
mkdir -p "$root/logs"
# SRA Toolkit, Bowtie2, and samtools must be on PATH.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$work" "$root/taiji/receipts"
run="$1"
case "$run" in
  SRR5822253) bytes=7454352295 ;;
  SRR5822288) bytes=4306311372 ;;
  SRR5822289) bytes=3436724569 ;;
  SRR5822290) bytes=4155730241 ;;
  SRR5822291) bytes=759994187 ;;
  *) echo "Unknown planned run $run" >&2; exit 2 ;;
esac
exec 9>"$work/$run.lock"
flock -n 9 || { echo "$run already running"; exit 2; }
if [ -s "$work/$run.names.bam" ] && [ -s "$root/taiji/receipts/$run.alignment.tsv" ]; then
  samtools quickcheck -v "$work/$run.names.bam"
  exit 0
fi
archive="$root/taiji/sra/$run.sra"
if [ ! -s "$archive" ] || [ "$(stat -c %s "$archive")" -ne "$bytes" ]; then
  echo "$run incomplete ($archive); transfer must finish before alignment" >&2
  exit 2
fi
exec > >(tee -a "$root/logs/TAIJI_ALIGN_$run.log") 2>&1
set -x
sra="$work/$run.sra"
if [ ! -s "$sra" ] || [ "$(stat -c %s "$sra")" -ne "$bytes" ]; then
  cp -- "$archive" "$sra.staging"
  mv -- "$sra.staging" "$sra"
fi
if [ ! -s "$root/taiji/receipts/$run.sra.sha256" ]; then
  vdb-validate "$sra" 2>&1 | tee "$root/logs/TAIJI_VDB_$run.log"
  sha256sum "$sra" > "$root/taiji/receipts/$run.sra.sha256"
fi
sha256sum --check "$root/taiji/receipts/$run.sra.sha256"
sra-stat --quick --xml "$sra" > "$root/taiji/receipts/$run.sra_stat.xml"
date -Is
mkdir -p "$work/$run.tmp"
fasterq-dump --split-spot --stdout --skip-technical --threads 3 --temp "$work/$run.tmp" \
  --seq-defline "@${run}.\$si/\$ri" --qual-defline '+' "$sra" \
  | python3 "$code/taiji_paired_stream.py" "$root/taiji/receipts/$run.reads.tsv" \
  | bowtie2 --very-sensitive-local --mm --seed 0 -p "${TAIJI_ALIGN_THREADS:-7}" -I 0 -X 2000 \
    --no-mixed --no-discordant -x "$assets/BOWTIE2_SHARED_INDEX/genome" --interleaved - \
  | samtools sort -n -@ 1 -m 384M -T "$work/$run.tmp/names" -o "$work/$run.partial.bam" -
python3 - "$root/taiji/receipts/$run.sra_stat.xml" "$root/taiji/receipts/$run.reads.tsv" <<'PY'
import csv
import sys
import xml.etree.ElementTree as ET
count = int(ET.parse(sys.argv[1]).getroot().attrib["spot_count"])
with open(sys.argv[2]) as file:
    aligned = sum(int(row["spots"]) for row in csv.DictReader(file,delimiter="\t"))
print(f"Streamed {aligned} paired spots (SRA spot count {count})")
assert aligned == count, (count,aligned)
PY
samtools quickcheck -v "$work/$run.partial.bam"
samtools flagstat -O tsv "$work/$run.partial.bam" > "$root/taiji/receipts/$run.alignment.tsv.partial"
mv -- "$work/$run.partial.bam" "$work/$run.names.bam"
mv -- "$root/taiji/receipts/$run.alignment.tsv.partial" "$root/taiji/receipts/$run.alignment.tsv"
date -Is
