#!/usr/bin/env bash

set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: refresh-hnswlib.sh [--archive FILE]

Materialize the pinned hnswlib source and ordered patch series in inst/include.
The separate inst/include/pforr directory is not modified.
EOF
}

fail() {
  printf 'refresh-hnswlib: %s\n' "$1" >&2
  exit 2
}

archive_file=
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
    *)
      fail "unknown argument: $1"
      ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)
checksums_file="${script_dir}/upstream-files.sha256"
materialized_dir="${repository_root}/inst/include"
[[ -d ${materialized_dir}/pforr ]] || \
  fail "expected separate pforr vendor directory is missing"

temporary_root=$(mktemp -d "${TMPDIR:-/tmp}/rcpphnsw-vendor-refresh.XXXXXX")
cleanup() {
  rm -r -- "${temporary_root}"
}
trap cleanup EXIT

reconstructed_dir="${temporary_root}/hnswlib"
reconstruct_args=(--output "${reconstructed_dir}")
verify_args=()
if [[ -n ${archive_file} ]]; then
  reconstruct_args+=(--archive "${archive_file}")
  verify_args+=(--archive "${archive_file}")
fi
"${script_dir}/reconstruct-hnswlib.sh" "${reconstruct_args[@]}"

while read -r _ relative_path; do
  header_name=${relative_path#hnswlib/}
  cp -- "${reconstructed_dir}/${header_name}" \
    "${materialized_dir}/${header_name}"
done < "${checksums_file}"

"${script_dir}/verify-hnswlib.sh" "${verify_args[@]}"
printf 'refreshed materialized hnswlib headers in %s\n' "${materialized_dir}"
