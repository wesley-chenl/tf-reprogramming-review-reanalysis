#!/usr/bin/env bash
# Replay the native Taiji v1.3.0 invocation using explicit local assets.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../../.." && pwd)
root=$(realpath "${PAIR_TAIJI_ROOT:?Set PAIR_TAIJI_ROOT}")
assets=$(realpath "${PAIR_TAIJI_ASSETS:?Set PAIR_TAIJI_ASSETS}")
binary=$(realpath "${TAIJI_BIN:?Set TAIJI_BIN to the v1.3.0 executable}")
compat=$(realpath "${TAIJI_COMPAT_LIB:?Set TAIJI_COMPAT_LIB to the legacy library directory}")
run="$root/work/native_taiji"
config="$repo/config/pairs/taiji"
mkdir -p "$run/tmp" "$root/taiji"
cp -- "$config/input.tsv" "$config/config.yaml" "$run/"
test -s "$root/taiji/TAIJI_RUN_QC.tsv"
source_peak="$root/work/peaks/SRR5822253/SRR5822253_peaks.narrowPeak"
target_peak="$root/work/peaks/cortex_pooled/cortex_pooled_peaks.narrowPeak"
genome="$assets/GRCm38.p6.genome.fa"
annotation="$assets/TAIJI_RESOURCES/gencode.vM24.annotation.gtf"
motifs="$assets/cisBP_mouse.meme"
index="$assets/TAIJI_mm10.genome_index"
rna="$root/taiji/rna"
(
  cd "$root"
  sha256sum "$binary" "$genome" "$annotation" "$motifs" "$index" "$source_peak" "$target_peak" \
    "$rna/GSM1074345.tsv" "$rna/GSM1074346.tsv" "$rna/GSM1074357.tsv"
) > "$run/native_inputs.sha256"
/usr/bin/bwrap --die-with-parent --new-session --unshare-user --unshare-pid --unshare-net --unshare-ipc --unshare-uts \
  --cap-drop ALL --ro-bind /usr /usr --symlink usr/bin /bin --symlink usr/lib /lib --symlink usr/lib64 /lib64 \
  --dev /dev --proc /proc --tmpfs /tmp --clearenv --setenv PATH /usr/bin --setenv LC_ALL C.UTF-8 \
  --setenv TERM dumb --setenv LD_LIBRARY_PATH /compat --hostname taiji-blind \
  --ro-bind "$config/sandbox.hosts" /etc/hosts --ro-bind "$compat" /compat --ro-bind "$binary" /opt/taiji \
  --ro-bind "$genome" /input/mm10.fa --ro-bind "$annotation" /input/mm10.gtf \
  --ro-bind "$motifs" /input/cisBP_mouse.meme --ro-bind "$index" /input/mm10.genome_index \
  --ro-bind "$rna/GSM1074345.tsv" /input/GSM1074345.tsv --ro-bind "$rna/GSM1074346.tsv" /input/GSM1074346.tsv \
  --ro-bind "$rna/GSM1074357.tsv" /input/GSM1074357.tsv --ro-bind "$source_peak" /input/source.narrowPeak \
  --ro-bind "$target_peak" /input/target.narrowPeak --bind "$run" /work --chdir /work \
  /opt/taiji run --master-ip 127.0.0.1 --port 41441 --config /work/config.yaml --select Output_Ranks -n 3 +RTS -N3 -M10G -RTS
test -s "$run/output/GeneRanks.tsv"
cp -- "$run/output/GeneRanks.tsv" "$root/taiji/TAIJI_NATIVE_GENE_RANKS.tsv"
