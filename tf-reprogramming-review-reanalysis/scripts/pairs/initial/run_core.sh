#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$1" && pwd)
out=$(mkdir -p "$2" && cd "$2" && pwd)
stage=$(mktemp -d)
trap 'rm -r -- "$stage"' EXIT
mkdir "$stage/in" "$stage/out"
for name in CLEAN_PAIRS.tsv SCORING_POLICY.json GSE53951_EFFECT.tsv GSE114809_EFFECT.tsv TRRUST_MOUSE_SIGNED.tsv; do
    install -m 444 "$root/inputs/$name" "$stage/in/$name"
done
expected=$(/usr/bin/python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["old_scores_sha256"])' "$root/locks/REGISTRATION.json")
( cd "$stage/in"; sha256sum -- * ) > "$out/CORE_INPUT_SHA256.tsv"
/usr/bin/bwrap --die-with-parent --new-session --unshare-user --unshare-pid --unshare-net --unshare-ipc --unshare-uts \
 --cap-drop ALL --ro-bind /usr /usr --symlink usr/bin /bin --symlink usr/lib /lib --symlink usr/lib64 /lib64 \
 --dev /dev --proc /proc --tmpfs /tmp --clearenv --setenv PATH /usr/bin --setenv LC_ALL C.UTF-8 --setenv OLD_SCORE_HASH "$expected" \
 --ro-bind "$stage/in" /input --ro-bind "$root/code/historical_score.py" /old.py \
 --ro-bind "$root/code/core_score.py" /core.py --bind "$stage/out" /out --chdir /input \
 /usr/bin/bash -c '/usr/bin/python3 -I -B /old.py && sha256sum /out/all_pair_scores.tsv | grep -qF "$OLD_SCORE_HASH" && /usr/bin/python3 -I -B /core.py'
install -m 644 "$stage/out/"* "$out/"
