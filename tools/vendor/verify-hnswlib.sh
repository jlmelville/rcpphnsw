#!/usr/bin/env bash

set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: verify-hnswlib.sh [--archive FILE] [--materialized-dir DIR]

Reconstruct the pinned, patched hnswlib headers and compare them byte-for-byte
with inst/include or a supplied disposable materialized directory.
EOF
}

fail() {
  printf 'verify-hnswlib: %s\n' "$1" >&2
  exit 2
}

archive_file=
materialized_dir=
while [[ $# -gt 0 ]]; do
  case $1 in
    --archive)
      [[ $# -ge 2 ]] || fail "--archive requires a file"
      archive_file=$2
      shift 2
      ;;
    --help)
      show_usage
      exit 0
      ;;
    --materialized-dir)
      [[ $# -ge 2 ]] || fail "--materialized-dir requires a directory"
      materialized_dir=$2
      shift 2
      ;;
    *)
      fail "unknown argument: $1"
      ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)
checksums_file="${script_dir}/upstream-files.sha256"
if [[ -z ${materialized_dir} ]]; then
  materialized_dir="${repository_root}/inst/include"
fi
[[ -d ${materialized_dir} ]] || fail "materialized directory not found: ${materialized_dir}"

temporary_root=$(mktemp -d "${TMPDIR:-/tmp}/rcpphnsw-vendor-verify.XXXXXX")
cleanup() {
  rm -r -- "${temporary_root}"
}
trap cleanup EXIT

reconstructed_dir="${temporary_root}/hnswlib"
reconstruct_args=(--output "${reconstructed_dir}")
if [[ -n ${archive_file} ]]; then
  reconstruct_args+=(--archive "${archive_file}")
fi
"${script_dir}/reconstruct-hnswlib.sh" "${reconstruct_args[@]}"

is_expected_header() {
  local candidate_name=$1
  local relative_path
  while read -r _ relative_path; do
    if [[ ${candidate_name} == "${relative_path#hnswlib/}" ]]; then
      return 0
    fi
  done < "${checksums_file}"
  return 1
}

shopt -s dotglob nullglob
for materialized_path in "${materialized_dir}"/*; do
  materialized_name=${materialized_path##*/}
  if [[ -L ${materialized_path} ]]; then
    fail "unexplained top-level vendored symlink: ${materialized_name}"
  elif [[ -f ${materialized_path} ]]; then
    is_expected_header "${materialized_name}" || \
      fail "unexplained top-level vendored file: ${materialized_name}"
  elif [[ -d ${materialized_path} ]]; then
    [[ ${materialized_name} == pforr ]] || \
      fail "unexplained top-level vendored directory: ${materialized_name}"
  else
    fail "unexplained top-level vendored entry: ${materialized_name}"
  fi
done
shopt -u dotglob nullglob

while read -r _ relative_path; do
  header_name=${relative_path#hnswlib/}
  reconstructed_header="${reconstructed_dir}/${header_name}"
  materialized_header="${materialized_dir}/${header_name}"
  [[ -f ${materialized_header} ]] || fail "materialized header is missing: ${header_name}"
  if ! cmp -s "${reconstructed_header}" "${materialized_header}"; then
    diff -u "${reconstructed_header}" "${materialized_header}" >&2 || true
    fail "materialized header differs from the pinned patch series: ${header_name}"
  fi
done < "${checksums_file}"

printf 'verified hnswlib %s headers against the pinned patch series\n' \
  "$(wc -l < "${checksums_file}" | tr -d ' ')"
