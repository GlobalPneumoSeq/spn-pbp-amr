#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 1 ]]; then
    printf 'usage: %s SAMPLE_DIRECTORY\n' "$0" >&2
    exit 2
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
bash "${script_dir}/PBP_AA_sampledir_to_MIC_20180710.sh" "$1" "$script_dir"
