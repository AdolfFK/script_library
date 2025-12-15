#!/usr/bin/env bash
#==============================================================================
# 脚本名: dusub_submission.sh
# 功能:  Batch-submit every *.sh script under <scripts_dir> via nohup.
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-15  adolf1  增加了head的描述
# 依赖:  bash≥4.0  curl≥7.0
# 用法:  ./dusub_submission.sh <scripts_dir>
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================


set -euo pipefail
shopt -s nullglob          # 无匹配时返回空数组，避免 *.sh 被当成字面量

readonly VERSION="0.0.1"
readonly SCRIPT_NAME=$(basename "$0")

#--------------------------------------------------
# Usage
#--------------------------------------------------
usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} <scripts_dir>

Batch-submit every *.sh script under <scripts_dir> via nohup.
Logs are written into <scripts_dir>/logs/.

Options:
  -h, --help    Show this help and exit
EOF
}

#--------------------------------------------------
# Argument parsing
#--------------------------------------------------
if [[ $# -ne 1 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    usage
    exit 0
fi

readonly SCRIPTS_DIR=$1
readonly LOG_DIR="${SCRIPTS_DIR}/logs"

#--------------------------------------------------
# Sanity checks
#--------------------------------------------------
if [[ ! -d ${SCRIPTS_DIR} ]]; then
    echo "[ERROR] Directory does not exist: ${SCRIPTS_DIR}" >&2
    exit 1
fi

mkdir -p "${LOG_DIR}"

mapfile -t SH_FILES < <(find "${SCRIPTS_DIR}" -maxdepth 1 -type f -name "*.sh" | sort)
if [[ ${#SH_FILES[@]} -eq 0 ]]; then
    echo "[WARN] No *.sh scripts found under ${SCRIPTS_DIR}" >&2
    exit 0
fi

#--------------------------------------------------
# Submit
#--------------------------------------------------
for script in "${SH_FILES[@]}"; do
    script_name=$(basename "${script}" .sh)
    log_file="${LOG_DIR}/${script_name}.log"

    nohup bash "${script}" > "${log_file}" 2>&1 &
    echo "[INFO] Submitted ${script}  -->  ${log_file}"
done

echo "[INFO] All scripts submitted. Logs under ${LOG_DIR}/"
