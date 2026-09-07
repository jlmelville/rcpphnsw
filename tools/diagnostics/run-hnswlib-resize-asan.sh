#!/usr/bin/env bash

set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: run-hnswlib-resize-asan.sh

Build and run the bounded hnswlib resize lifecycle diagnostic under
AddressSanitizer. The runner requires a Linux linker supporting --wrap.
EOF
}

if [[ ${1:-} == "--help" ]]; then
  show_usage
  exit 0
fi
if [[ $# -ne 0 ]]; then
  show_usage >&2
  exit 2
fi
if [[ $(uname -s) != Linux ]]; then
  printf '%s\n' 'run-hnswlib-resize-asan.sh requires Linux' >&2
  exit 77
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)
source_file="${script_dir}/hnswlib-resize-asan.cpp"
compiler=${CXX:-g++}

if ! command -v "${compiler}" >/dev/null 2>&1; then
  printf 'required C++ compiler not found: %s\n' "${compiler}" >&2
  exit 127
fi

build_dir=$(mktemp -d "${TMPDIR:-/tmp}/rcpphnsw-resize-asan.XXXXXX")
binary_file="${build_dir}/hnswlib-resize-asan"
cleanup() {
  rm -f -- "${binary_file}"
  rmdir -- "${build_dir}"
}
trap cleanup EXIT

"${compiler}" \
  -std=c++17 \
  -O1 \
  -g \
  -fno-omit-frame-pointer \
  -fsanitize=address \
  -pthread \
  -Wl,--wrap=realloc \
  -DNO_MANUAL_VECTORIZATION \
  "-I${repository_root}/inst/include" \
  "${source_file}" \
  -o "${binary_file}"

asan_options=${ASAN_OPTIONS:-detect_leaks=0:halt_on_error=1}
ASAN_OPTIONS=${asan_options} "${binary_file}" zero

failure_cases=0
for resize_mode in growth shrink; do
  count_output=$(ASAN_OPTIONS=${asan_options} "${binary_file}" "${resize_mode}")
  allocation_count=${count_output##*=}
  if [[ ! ${allocation_count} =~ ^[1-9][0-9]*$ ]]; then
    printf 'invalid allocation count from %s: %s\n' \
      "${resize_mode}" "${count_output}" >&2
    exit 2
  fi
  failure_index=1
  while [[ ${failure_index} -le ${allocation_count} ]]; do
    ASAN_OPTIONS=${asan_options} \
      "${binary_file}" "${resize_mode}" "${failure_index}"
    failure_cases=$((failure_cases + 1))
    failure_index=$((failure_index + 1))
  done
done

printf 'resize lifecycle diagnostic passed %d injected failures\n' \
  "${failure_cases}"
