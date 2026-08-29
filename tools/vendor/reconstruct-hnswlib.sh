#!/usr/bin/env bash

set -euo pipefail

show_usage() {
  cat <<'EOF'
Usage: reconstruct-hnswlib.sh --output DIR [--archive FILE]

Reconstruct the patched hnswlib header set in a new or empty directory.
Without --archive, download the pinned archive into temporary storage.
EOF
}

fail() {
  printf 'reconstruct-hnswlib: %s\n' "$1" >&2
  exit 2
}

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    fail "required command not found: $1"
  fi
}

manifest_value() {
  local key=$1
  local value
  value=$(awk -F ': ' -v key="${key}" '$1 == key { print substr($0, length($1) + 3) }' "${manifest_file}")
  if [[ -z ${value} || ${value} == *$'\n'* ]]; then
    fail "manifest key must occur exactly once with a value: ${key}"
  fi
  printf '%s\n' "${value}"
}

sha256_file() {
  local input_file=$1
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "${input_file}" | awk '{ print $1 }'
  elif command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "${input_file}" | awk '{ print $1 }'
  else
    fail "required SHA-256 tool not found: sha256sum or shasum"
  fi
}

output_dir=
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
    --output)
      [[ $# -ge 2 ]] || fail "--output requires a directory"
      output_dir=$2
      shift 2
      ;;
    *)
      fail "unknown argument: $1"
      ;;
  esac
done

[[ -n ${output_dir} ]] || fail "--output is required"
if [[ -e ${output_dir} ]]; then
  [[ -d ${output_dir} ]] || fail "output exists and is not a directory: ${output_dir}"
  [[ -z $(ls -A "${output_dir}") ]] || fail "output directory is not empty: ${output_dir}"
else
  mkdir -p -- "${output_dir}"
fi
if [[ -n ${archive_file} && ! -f ${archive_file} ]]; then
  fail "archive is not a regular file: ${archive_file}"
fi

require_command awk
require_command git
require_command mktemp
require_command tar

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
manifest_file="${script_dir}/upstream.yml"
checksums_file="${script_dir}/upstream-files.sha256"
patch_directory="${script_dir}/patches"
series_file="${patch_directory}/series"

[[ -f ${manifest_file} ]] || fail "missing manifest: ${manifest_file}"
[[ -f ${checksums_file} ]] || fail "missing header checksums: ${checksums_file}"
[[ -f ${series_file} ]] || fail "missing patch series: ${series_file}"

archive_url=$(manifest_value archive_url)
archive_root=$(manifest_value archive_root)
archive_sha256=$(manifest_value archive_sha256)
source_subdirectory=$(manifest_value source_subdirectory)

temporary_root=$(mktemp -d "${TMPDIR:-/tmp}/rcpphnsw-vendor.XXXXXX")
cleanup() {
  rm -r -- "${temporary_root}"
}
trap cleanup EXIT

if [[ -z ${archive_file} ]]; then
  require_command curl
  archive_file="${temporary_root}/upstream.tar.gz"
  curl --location --fail --silent --show-error \
    --output "${archive_file}" "${archive_url}"
fi

observed_archive_sha256=$(sha256_file "${archive_file}")
if [[ ${observed_archive_sha256} != "${archive_sha256}" ]]; then
  fail "archive SHA-256 mismatch: expected ${archive_sha256}, got ${observed_archive_sha256}"
fi

extract_directory="${temporary_root}/extract"
mkdir -p -- "${extract_directory}"
tar -xzf "${archive_file}" -C "${extract_directory}"
upstream_directory="${extract_directory}/${archive_root}/${source_subdirectory}"
[[ -d ${upstream_directory} ]] || fail "archive is missing ${archive_root}/${source_subdirectory}"

header_count=0
while read -r expected_sha256 relative_path extra; do
  [[ -n ${expected_sha256} ]] || continue
  [[ -z ${extra:-} ]] || fail "invalid checksum row: ${relative_path} ${extra}"
  case ${relative_path} in
    "${source_subdirectory}/"*.h) ;;
    *) fail "invalid vendored header path: ${relative_path}" ;;
  esac
  header_name=${relative_path#"${source_subdirectory}/"}
  [[ ${header_name} != */* ]] || fail "nested hnswlib header is not supported: ${relative_path}"
  upstream_header="${extract_directory}/${archive_root}/${relative_path}"
  [[ -f ${upstream_header} ]] || fail "archive is missing header: ${relative_path}"
  observed_sha256=$(sha256_file "${upstream_header}")
  if [[ ${observed_sha256} != "${expected_sha256}" ]]; then
    fail "header SHA-256 mismatch for ${relative_path}"
  fi
  cp -- "${upstream_header}" "${output_dir}/${header_name}"
  header_count=$((header_count + 1))
done < "${checksums_file}"
[[ ${header_count} -gt 0 ]] || fail "header checksum list is empty"

patch_count=0
while IFS= read -r patch_name || [[ -n ${patch_name} ]]; do
  [[ -n ${patch_name} ]] || fail "blank patch-series entry"
  case ${patch_name} in
    */* | *[!A-Za-z0-9._-]*) fail "invalid patch-series entry: ${patch_name}" ;;
  esac
  patch_file="${patch_directory}/${patch_name}"
  [[ -f ${patch_file} ]] || fail "patch-series entry is missing: ${patch_name}"
  git -C "${output_dir}" apply --no-index --check --whitespace=error-all "${patch_file}"
  git -C "${output_dir}" apply --no-index --whitespace=error-all "${patch_file}"
  patch_count=$((patch_count + 1))
done < "${series_file}"

for patch_file in "${patch_directory}"/*.patch; do
  patch_name=${patch_file##*/}
  listed=false
  while IFS= read -r series_name || [[ -n ${series_name} ]]; do
    if [[ ${patch_name} == "${series_name}" ]]; then
      listed=true
      break
    fi
  done < "${series_file}"
  [[ ${listed} == true ]] || fail "patch is not listed in series: ${patch_name}"
done
[[ ${patch_count} -gt 0 ]] || fail "patch series is empty"
