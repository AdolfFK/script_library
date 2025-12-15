#!/usr/bin/env bash
#==============================================================================
# 脚本名: create_bio_user_change_passwd.sh
# 功能:  创建生信用户并初始化工作目录
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-15  adolf1  增加head说明
# 依赖:  bash≥4.0  curl≥7.0
# 用法:  ############支持明文 第一次登录后强制修改############
#       Usage:  sudo ./create_bio_user_change_passwd.sh  USERNAME  PASSWORD

# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================

set -euo pipefail

# ---------- 参数检查 ----------
if [[ $# -ne 2 ]]; then
    echo "Usage: sudo $0 <USERNAME> <PASSWORD>"
    exit 1
fi

USER=$1
PASS=$2

# ---------- 权限检查 ----------
if [[ $EUID -ne 0 ]]; then
    echo "Error: This script must be run as root!"
    exit 1
fi

# ---------- 用户是否存在 ----------
if id "$USER" &>/dev/null; then
    echo "Error: User '$USER' already exists."
    exit 2
fi

# ---------- 1. 建用户并加入附加组 ----------
useradd -m -G bioinformaticsG "$USER"
echo "${USER}:${PASS}" | chpasswd

# ---------- 2. 创建工作目录 ----------
WORK_PARENT="/data/workdata"
WORK_DIR="${WORK_PARENT}/${USER}"

mkdir -p "$WORK_DIR"
chown -R "${USER}:bioinformaticsG" "$WORK_DIR"

# ---------- 3. 强制首次登录改密 ----------
passwd -e "$USER" >/dev/null

echo "----------------------------------------"
echo "User '$USER' created successfully."
echo "Home : /home/${USER}"
echo "Work : ${WORK_DIR}"
echo "----------------------------------------"
