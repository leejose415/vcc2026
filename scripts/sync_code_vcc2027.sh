#!/usr/bin/env bash
# Ship committed code from this VM (vcc2026) to vcc2027 through the bucket: git bundle -> gs://vcc2026/vm/repo/ ->
# fast-forward the checkout on vcc2027. Only committed work travels; vcc2027's working tree must be clean of edits to
# tracked files. Usage: bash scripts/sync_code_vcc2027.sh
set -euo pipefail
cd "$(dirname "$0")/.."
branch=$(git branch --show-current)
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
git bundle create "$tmp/vcc.bundle" --all 2>/dev/null
gcloud storage cp "$tmp/vcc.bundle" gs://vcc2026/vm/repo/vcc.bundle --quiet
gcloud compute ssh vcc2027 --zone us-central1-c --quiet --command \
  "set -e; cd ~/vcc; gcloud storage cp gs://vcc2026/vm/repo/vcc.bundle /tmp/vcc.bundle --quiet; \
   git fetch -q /tmp/vcc.bundle $branch:refs/remotes/bundle/$branch; git merge -q --ff-only bundle/$branch; git log --oneline -1"
