#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
api_url=${ZENODO_API_URL:-https://zenodo.org/api}
access_token=${ZENODO_ACCESS_TOKEN:-}
deposit_id=${1:-}

if [[ -z "$access_token" ]]; then
  echo "Set ZENODO_ACCESS_TOKEN to a Zenodo token with deposit:write scope." >&2
  exit 2
fi

archive="$repo_dir/zenodo/ehi_data_release_2_mags.tar.gz"
archive_parts=("$repo_dir"/zenodo/ehi_data_release_2_mags.tar.gz.part-*)
readme="$repo_dir/zenodo/readme"
metadata="$repo_dir/zenodo/metadata.json"
draft_state="$repo_dir/zenodo/zenodo_draft.json"

if [[ -e "${archive_parts[0]}" ]]; then
  payloads=("$readme" "${archive_parts[@]}")
else
  payloads=("$readme" "$archive")
fi

for path in "$metadata" "${payloads[@]}"; do
  if [[ ! -f "$path" ]]; then
    echo "Missing required payload: $path" >&2
    exit 2
  fi
done

temp_response=$(mktemp -t ehr02-zenodo.XXXXXX)
trap 'rm -f "$temp_response"' EXIT
auth_header="Authorization: Bearer $access_token"

if [[ -z "$deposit_id" ]]; then
  curl --fail --silent --show-error \
    --request POST \
    --header "$auth_header" \
    --header 'Content-Type: application/json' \
    --data '{}' \
    "$api_url/deposit/depositions" > "$temp_response"
  deposit_id=$(jq -er '.id' "$temp_response")
  cp "$temp_response" "$draft_state"
  echo "Created Zenodo draft $deposit_id" >&2
  echo "If this script stops, rerun it with the draft ID: $deposit_id" >&2
fi

curl --fail --silent --show-error \
  --request PUT \
  --header "$auth_header" \
  --header 'Content-Type: application/json' \
  --data-binary "@$metadata" \
  "$api_url/deposit/depositions/$deposit_id" > "$temp_response"
cp "$temp_response" "$draft_state"

bucket_url=$(jq -er '.links.bucket' "$temp_response")
reserved_doi=$(jq -er '.metadata.prereserve_doi.doi' "$temp_response")
draft_url=$(jq -er '.links.html' "$temp_response")

file_md5() {
  if command -v md5 >/dev/null 2>&1; then
    md5 -q "$1"
  else
    md5sum "$1" | awk '{print $1}'
  fi
}

for path in "${payloads[@]}"; do
  filename=$(basename -- "$path")
  local_size=$(wc -c < "$path" | tr -d ' ')
  remote_size=$(jq -r --arg filename "$filename" \
    '.files[]? | select(.filename == $filename) | .filesize' "$draft_state")
  remote_md5=$(jq -r --arg filename "$filename" \
    '.files[]? | select(.filename == $filename) | .checksum' "$draft_state")

  if [[ "$remote_size" == "$local_size" && -n "$remote_md5" ]]; then
    local_md5=$(file_md5 "$path")
    if [[ "$remote_md5" == "$local_md5" ]]; then
      echo "Already uploaded and verified: $filename" >&2
      continue
    fi
  fi

  echo "Uploading $filename to draft $deposit_id" >&2
  curl --fail --silent --show-error \
    --retry 5 \
    --retry-all-errors \
    --request PUT \
    --header "$auth_header" \
    --upload-file "$path" \
    "$bucket_url/$filename" > "$temp_response"
done

curl --fail --silent --show-error \
  --header "$auth_header" \
  "$api_url/deposit/depositions/$deposit_id" > "$draft_state"

jq -n \
  --arg deposit_id "$deposit_id" \
  --arg doi "$reserved_doi" \
  --arg draft_url "$draft_url" \
  '{deposit_id: $deposit_id, reserved_doi: $doi, draft_url: $draft_url, published: false}'

echo "Draft only: review the Zenodo page before publishing." >&2
