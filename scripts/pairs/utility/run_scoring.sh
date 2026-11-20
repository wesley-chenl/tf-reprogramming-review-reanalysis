#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$1" && pwd)
out=$(mkdir -p "$2" && cd "$2" && pwd)
stage=$(mktemp -d)
trap 'rm -r -- "$stage"' EXIT
mkdir "$stage/in" "$stage/out"
install -m 444 "$root/design_public/SCORING_POLICY.json" "$root/design_public/CLEAN_PAIRS.tsv" "$stage/in/"
for name in GSE53951_EFFECT.tsv GSE114809_EFFECT.tsv TRRUST_MOUSE_SIGNED.tsv; do
  if test -f "$root/inputs_scoring/$name"; then
    install -m 444 "$root/inputs_scoring/$name" "$stage/in/$name"
  fi
done
(
  cd "$stage/in"
  sha256sum -- *
) > "$out/SCORING_INPUT_SHA256.txt"
sha256sum "$root/code/score.py" > "$out/SCORER_SHA256.txt"
test -x /usr/bin/bwrap
/usr/bin/bwrap --die-with-parent --new-session --unshare-pid --unshare-net --unshare-ipc --unshare-uts \
  --cap-drop ALL --ro-bind /usr /usr --symlink usr/bin /bin --symlink usr/lib /lib --symlink usr/lib64 /lib64 \
  --dev /dev --proc /proc --tmpfs /tmp --clearenv --setenv PATH /usr/bin --setenv LC_ALL C.UTF-8 \
  --ro-bind "$stage/in" /input --ro-bind "$root/code/score.py" /scorer.py --bind "$stage/out" /out \
  --chdir /input /usr/bin/python3 -I -B /scorer.py
install -m 644 "$stage/out/"* "$out/"
printf '%s\n' 'bubblewrap: private network and PID namespaces; isolated filesystem; score source and exact allowlist mounted read-only; Python audited reads.' > "$out/ISOLATION.txt"
