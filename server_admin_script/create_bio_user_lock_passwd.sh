
#!/usr/bin/env bash
#==============================================================================
# 脚本名: create_bio_user_lock_passwd.sh
# 功能:  创建生信用户，允许正常登录，但禁止用户自行修改密码
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2026-11-04  修复密码锁定问题，允许 SSH 正常登录
#
# 依赖:  bash >= 4.0
#
# 用法:
#       sudo ./create_bio_user_lock_passwd.sh USERNAME PASSWORD
#
# 特点:
#       1. 允许 SSH 密码登录
#       2. 禁止普通用户自行修改密码
#       3. root 可以随时修改用户密码
#       4. 不锁定用户账号
#       5. 不设置密码过期时间
#
# 仓库: https://github.com/AdolfFK/script_library
#==============================================================================

set -euo pipefail

# ---------- 参数检查 ----------
if [[ $# -ne 2 ]]; then
    echo "Usage: sudo $0 <USERNAME> <PASSWORD>"
    exit 1
fi

USER="$1"
PASS="$2"

# ---------- 权限检查 ----------
if [[ $EUID -ne 0 ]]; then
    echo "Error: This script must be run as root!"
    exit 1
fi

# ---------- 用户名检查 ----------
if [[ ! "$USER" =~ ^[a-zA-Z_][a-zA-Z0-9_-]*$ ]]; then
    echo "Error: Invalid username: $USER"
    exit 1
fi

# ---------- 密码检查 ----------
if [[ -z "$PASS" ]]; then
    echo "Error: Password cannot be empty!"
    exit 1
fi

# ---------- 用户是否存在 ----------
if id "$USER" &>/dev/null; then
    echo "Error: User '$USER' already exists."
    exit 2
fi

# ---------- 检查用户组 ----------
if ! getent group bioinformaticsG >/dev/null; then
    echo "Error: Group 'bioinformaticsG' does not exist."
    exit 1
fi

# ---------- 1. 创建用户并加入附加组 ----------
useradd -m -G bioinformaticsG "$USER"

# 设置初始密码
echo "${USER}:${PASS}" | chpasswd

# ---------- 2. 设置密码策略 ----------

# 不锁定密码，允许正常 SSH 密码认证
# 禁止普通用户在 99999 天内修改密码
# 同时取消密码过期限制
chage -d "$(date +%F)" \
      -m 99999 \
      -M 99999 \
      -W 0 \
      -I -1 \
      -E -1 "$USER" >/dev/null

# ---------- 3. 创建工作目录 ----------
WORK_PARENT="/data/workdata"
WORK_DIR="${WORK_PARENT}/${USER}"

mkdir -p "$WORK_DIR"
chown -R "${USER}:bioinformaticsG" "$WORK_DIR"

# ---------- 4. 检查用户状态 ----------
echo
echo "----------------------------------------"
echo "User '$USER' created successfully."
echo "----------------------------------------"
echo "Home : /home/${USER}"
echo "Work : ${WORK_DIR}"
echo "Group: bioinformaticsG"
echo
echo "Password status:"
passwd -S "$USER"

echo
echo "Password aging:"
chage -l "$USER"

echo "----------------------------------------"
echo "SSH login       : ENABLED (subject to SSH/PAM policy)"
echo "Password change : RESTRICTED"
echo "Account lock    : NO"
echo "Root reset      : ENABLED"
echo "----------------------------------------"

# 清理密码变量
unset PASS
