#!/usr/bin/env bash
set -euo pipefail

# Pushes branches listed in good_branches.txt to target remote "vgk2"
# while skipping any branch listed in bad_branches.txt.
#
# This script does NOT rely on local branch history. For each branch, it fetches
# directly from the source remote into a temporary bare repo, then pushes to target.
#
# Optional environment variables:
# - GOOD_FILE: path to good branches list (default: good_branches.txt)
# - BAD_FILE: path to bad branches list (default: bad_branches.txt)
# - SOURCE_REMOTE: source remote name in this repo (default: origin)
# - SOURCE_REMOTE_URL: explicit source remote URL (overrides SOURCE_REMOTE)
# - TARGET_REMOTE: target remote name in this repo (default: vgk2)
# - TARGET_NAMESPACE: namespace used to derive target URL from origin (default: vgk2)
# - TARGET_REMOTE_URL: explicit target URL (default: git@gitlab.opengeosys.org:VinayGK/ogs-2.git)
# - DRY_RUN: if set to "1", prints commands without executing fetch/push

GOOD_FILE="${GOOD_FILE:-good_branches.txt}"
BAD_FILE="${BAD_FILE:-bad_branches.txt}"
SOURCE_REMOTE="${SOURCE_REMOTE:-origin}"
SOURCE_REMOTE_URL="${SOURCE_REMOTE_URL:-}"
TARGET_REMOTE="${TARGET_REMOTE:-vgk2}"
TARGET_NAMESPACE="${TARGET_NAMESPACE:-vgk2}"
TARGET_REMOTE_URL="${TARGET_REMOTE_URL:-git@gitlab.opengeosys.org:VinayGK/ogs-2.git}"
DRY_RUN="${DRY_RUN:-0}"

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "Error: run this script inside a git repository." >&2
  exit 1
fi

if [[ ! -f "${GOOD_FILE}" ]]; then
  echo "Error: ${GOOD_FILE} not found." >&2
  exit 1
fi

derive_remote_url() {
  local origin_url repo host path
  origin_url="$(git remote get-url origin)"

  if [[ "${origin_url}" =~ ^https?://([^/]+)/(.+)$ ]]; then
    host="${BASH_REMATCH[1]}"
    path="${BASH_REMATCH[2]}"
    repo="${path##*/}"
    echo "https://${host}/${TARGET_NAMESPACE}/${repo}"
    return
  fi

  if [[ "${origin_url}" =~ ^git@([^:]+):(.+)$ ]]; then
    host="${BASH_REMATCH[1]}"
    path="${BASH_REMATCH[2]}"
    repo="${path##*/}"
    echo "git@${host}:${TARGET_NAMESPACE}/${repo}"
    return
  fi

  echo "Error: unsupported origin URL format: ${origin_url}" >&2
  exit 1
}

ensure_target_remote() {
  local remote_url existing_url

  if git remote get-url "${TARGET_REMOTE}" >/dev/null 2>&1; then
    existing_url="$(git remote get-url "${TARGET_REMOTE}")"
    if [[ -n "${TARGET_REMOTE_URL}" && "${TARGET_REMOTE_URL}" != "${existing_url}" ]]; then
      echo "Updating remote '${TARGET_REMOTE}' URL:"
      echo "  from: ${existing_url}"
      echo "    to: ${TARGET_REMOTE_URL}"
      git remote set-url "${TARGET_REMOTE}" "${TARGET_REMOTE_URL}"
      existing_url="${TARGET_REMOTE_URL}"
    fi
    echo "Using existing remote '${TARGET_REMOTE}': ${existing_url}"
    return
  fi

  if [[ -n "${TARGET_REMOTE_URL}" ]]; then
    remote_url="${TARGET_REMOTE_URL}"
  else
    remote_url="$(derive_remote_url)"
  fi

  echo "Adding remote '${TARGET_REMOTE}' -> ${remote_url}"
  git remote add "${TARGET_REMOTE}" "${remote_url}"
}

resolve_source_url() {
  if [[ -n "${SOURCE_REMOTE_URL}" ]]; then
    SOURCE_REMOTE_URL_RESOLVED="${SOURCE_REMOTE_URL}"
    return
  fi

  if git remote get-url "${SOURCE_REMOTE}" >/dev/null 2>&1; then
    SOURCE_REMOTE_URL_RESOLVED="$(git remote get-url "${SOURCE_REMOTE}")"
    return
  fi

  echo "Error: source remote '${SOURCE_REMOTE}' not found and SOURCE_REMOTE_URL is not set." >&2
  exit 1
}

check_remote_access() {
  if [[ "${DRY_RUN}" == "1" ]]; then
    return
  fi

  if ! git ls-remote --heads "${SOURCE_REMOTE_URL_RESOLVED}" >/dev/null 2>&1; then
    echo "Error: cannot access source remote URL: ${SOURCE_REMOTE_URL_RESOLVED}" >&2
    exit 1
  fi

  if ! git ls-remote --heads "${TARGET_REMOTE_URL_RESOLVED}" >/dev/null 2>&1; then
    echo "Error: cannot access target remote URL: ${TARGET_REMOTE_URL_RESOLVED}" >&2
    echo "Make sure the fork exists and you have push permission." >&2
    exit 1
  fi
}

is_bad_branch() {
  local branch="$1"
  if [[ ! -f "${BAD_FILE}" ]]; then
    return 1
  fi
  grep -Fxq "${branch}" "${BAD_FILE}"
}

list_good_branches() {
  sed -e 's/[[:space:]]*$//' -e 's/^[[:space:]]*//' "${GOOD_FILE}" | awk 'NF'
}

source_branch_exists() {
  local branch="$1"
  git ls-remote --exit-code "${SOURCE_REMOTE_URL_RESOLVED}" "refs/heads/${branch}" >/dev/null 2>&1
}

push_branch_via_temp_repo() {
  local branch="$1"
  local fetch_cmd=(git -C "${TMP_REPO_DIR}" fetch --force --no-tags "${SOURCE_REMOTE_URL_RESOLVED}" "refs/heads/${branch}:refs/heads/${branch}")
  local push_cmd=(git -C "${TMP_REPO_DIR}" push "${TARGET_REMOTE_URL_RESOLVED}" "refs/heads/${branch}:refs/heads/${branch}")

  if [[ "${DRY_RUN}" == "1" ]]; then
    printf 'DRY_RUN: %q ' "${fetch_cmd[@]}"
    printf '\n'
    printf 'DRY_RUN: %q ' "${push_cmd[@]}"
    printf '\n'
    return
  fi

  "${fetch_cmd[@]}"
  "${push_cmd[@]}"
}

ensure_target_remote
resolve_source_url
TARGET_REMOTE_URL_RESOLVED="$(git remote get-url "${TARGET_REMOTE}")"
check_remote_access

TMP_REPO_DIR="$(mktemp -d "${TMPDIR:-/tmp}/push-good-branches.XXXXXX")"
trap 'rm -rf "${TMP_REPO_DIR}"' EXIT

git init --bare "${TMP_REPO_DIR}" >/dev/null

pushed_count=0
skipped_count=0
missing_count=0
failed_count=0

while IFS= read -r branch; do
  [[ -z "${branch}" ]] && continue

  if is_bad_branch "${branch}"; then
    echo "Skipping '${branch}' (listed in ${BAD_FILE})"
    skipped_count=$((skipped_count + 1))
    continue
  fi

  if ! source_branch_exists "${branch}"; then
    echo "Skipping '${branch}' (not found on source remote)"
    missing_count=$((missing_count + 1))
    continue
  fi

  echo "Syncing '${branch}' from source to target..."
  if push_branch_via_temp_repo "${branch}"; then
    pushed_count=$((pushed_count + 1))
  else
    echo "Failed to sync '${branch}'"
    failed_count=$((failed_count + 1))
  fi
done < <(list_good_branches)

echo
echo "Done."
echo "Pushed: ${pushed_count}"
echo "Skipped (bad list): ${skipped_count}"
echo "Skipped (missing on source): ${missing_count}"
echo "Failed: ${failed_count}"

if [[ "${failed_count}" -gt 0 ]]; then
  exit 1
fi
