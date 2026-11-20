#!/usr/bin/env bash
# Mount only score inputs, code, software and blank output in an offline namespace.
set -euo pipefail
if [ "$(id -u)" -eq 0 ] || [ "$#" -lt 6 ]; then
  echo 'Usage (ordinary Linux user): run_isolated.sh INPUTS CODE OUTPUT RUNTIME AUTHORS|- COMMAND [ARGS...]' >&2
  exit 2
fi
inp=$(realpath "$1")
code=$(realpath "$2")
out=$(mkdir -p "$3" && realpath "$3")
runtime=$(realpath "$4")
source=$5
shift 5
test -d "$inp" && test -d "$code" && test -d "$runtime"
bindings=()
if [ "$source" != "-" ]; then
  source=$(realpath "$source")
  test -d "$source"
  bindings+=(--ro-bind "$source" /authors)
fi
command=( "$@" )
/usr/bin/bwrap --die-with-parent --new-session --unshare-user --unshare-net --unshare-pid --unshare-ipc --unshare-uts \
  --cap-drop ALL --ro-bind /usr /usr --symlink usr/bin /bin --symlink usr/lib /lib --symlink usr/lib64 /lib64 \
  --dev /dev --proc /proc --tmpfs /tmp --clearenv --setenv HOME /tmp --setenv LC_ALL C.UTF-8 \
  --setenv PATH "$runtime/bin:/usr/bin" --setenv OPENBLAS_NUM_THREADS 1 --setenv OMP_NUM_THREADS 1 \
  --setenv MKL_NUM_THREADS 1 --setenv PYTHONHASHSEED 0 --setenv PYTHONDONTWRITEBYTECODE 1 \
  --ro-bind "$runtime" "$runtime" --ro-bind "$code" /code \
  --ro-bind "$inp" /inputs --symlink /inputs /input --bind "$out" /out --chdir /inputs \
  "${bindings[@]}" "${command[@]}"
