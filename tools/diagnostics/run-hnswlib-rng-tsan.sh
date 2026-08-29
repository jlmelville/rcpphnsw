#!/usr/bin/env bash

set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: run-hnswlib-rng-tsan.sh [rounds [threads [points-per-thread [M [seed]]]]]

Build and run the bounded hnswlib insertion diagnostic under ThreadSanitizer.
Defaults: 4 rounds, 8 threads, 128 points per thread, M = 8, seed = 100.
EOF
}

if [[ ${1:-} == "--help" ]]; then
  show_usage
  exit 0
fi
if [[ $# -gt 5 ]]; then
  show_usage >&2
  exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)
source_file="${script_dir}/hnswlib-rng-tsan.cpp"
compiler=${CXX:-clang++}

if ! command -v "${compiler}" >/dev/null 2>&1; then
  printf 'required C++ compiler not found: %s\n' "${compiler}" >&2
  exit 127
fi

build_dir=$(mktemp -d "${TMPDIR:-/tmp}/rcpphnsw-tsan.XXXXXX")
binary_file="${build_dir}/hnswlib-rng-tsan"
cleanup() {
  rm -f -- "${binary_file}"
  rmdir -- "${build_dir}"
}
trap cleanup EXIT

"${compiler}" \
  -std=c++17 \
  -O1 \
  -fno-inline \
  -g \
  -fno-omit-frame-pointer \
  -fsanitize=thread \
  -pthread \
  -DNO_MANUAL_VECTORIZATION \
  "-I${repository_root}/inst/include" \
  "${source_file}" \
  -o "${binary_file}"

TSAN_OPTIONS=${TSAN_OPTIONS:-halt_on_error=1:exitcode=66:detect_deadlocks=0} \
  "${binary_file}" "$@"
