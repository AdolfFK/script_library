#!/usr/bin/env bash
#==============================================================================
# 脚本名: sync2gh.sh
# 功能:  提交并推送本地 Git 仓库目录 ($LOCAL_DIR) 中的改动；
#        新增“版本管理”：通过 -v 参数指定版本号，若与当前 VERSION 文件不同则创建新版本。
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-11  adolf1  Version 2：增加版本控制逻辑
# 依赖:  bash≥4.0  git≥2.0
# 用法:
#        # 相同版本内日常更新（原用法）
#        bash sync2gh.sh "日常更新"
#        bash sync2gh.sh script.py "修复 typo"
#
#        # 发布新版本
#        bash sync2gh.sh -v 1.1.0 "发布 v1.1.0：新增某某功能"
#        bash sync2gh.sh -v 1.1.0 script.py "发布 v1.1.0：更新脚本"
#
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================

set -euo pipefail

# --- 配置区 ---
REPO_URL="https://github.com/AdolfFK/script_library"
LOCAL_DIR="/Users/adolf1/Documents/Work/Script_Library"
BRANCH_NAME="main"
VERSION_FILE="VERSION"          # 版本号文件，位于仓库根

# --- 解析新增 -v 参数 ---
VERSION_ARG=""
while getopts "v:" opt; do
  case $opt in
    v) VERSION_ARG="$OPTARG" ;;
    *) echo "❌ 未知选项: -$OPTARG" >&2; exit 1 ;;
  esac
done
shift $((OPTIND-1))             # 移除已解析的选项，保留原位置参数

# --- 原有模式参数检查 ---
MODE=""
TARGET_FULL_PATH=""
COMMIT_MSG=""

if [[ $# -eq 1 ]]; then
    MODE="BULK"; COMMIT_MSG="$1"
    echo "💡 模式: 整体同步${VERSION_ARG:+ （版本: $VERSION_ARG）}。"
elif [[ $# -eq 2 ]]; then
    MODE="SINGLE_FILE"
    TARGET_PATH_ARG="$1"; COMMIT_MSG="$2"
    [[ -e "$TARGET_PATH_ARG" ]] || { echo "❌ 指定路径不存在: $TARGET_PATH_ARG"; exit 1; }
    TARGET_FULL_PATH=$(realpath "$TARGET_PATH_ARG")
    echo "💡 模式: 单文件/目录提交 -> $TARGET_FULL_PATH${VERSION_ARG:+ （版本: $VERSION_ARG）}。"
else
    echo "❌ 用法错误！"
    echo "整体同步: $0 [-v <版本>] '<提交信息>'"
    echo "单路径提交: $0 [-v <版本>] <文件/目录> '<提交信息>'"
    exit 1
fi

# --- 确保仓库存在并更新 ---
mkdir -p "$LOCAL_DIR"
cd "$LOCAL_DIR"

if [[ ! -d .git ]]; then
    echo "📥 首次运行，正在克隆仓库..."
    git clone --depth 1 "$REPO_URL" .
    git checkout "$BRANCH_NAME"
fi

echo "🔄 拉取远程最新内容..."
git pull origin "$BRANCH_NAME"

# --- 版本处理逻辑 ---
if [[ -n "$VERSION_ARG" ]]; then
    # 读取旧版本（若文件不存在则视为空）
    OLD_VER=$(test -f "$VERSION_FILE" && cat "$VERSION_FILE" || echo "")
    if [[ "$VERSION_ARG" != "$OLD_VER" ]]; then
        echo "🆕 检测到新版本号：$OLD_VER -> $VERSION_ARG"
        echo "$VERSION_ARG" > "$VERSION_FILE"
        # 确保 VERSION 文件被加入暂存
        git add "$VERSION_FILE"
    else
        echo "ℹ️  版本号未变化，执行日常更新。"
    fi
fi

# --- 暂存改动 ---
echo "📦 正在暂存改动..."
if [[ "$MODE" == "BULK" ]]; then
    git add -A .
else
    git add -A "$TARGET_FULL_PATH"
fi

# --- 若无改动则退出 ---
if git diff --staged --quiet; then
    echo "ℹ️  无改动，无需提交。"
    exit 0
fi

# --- 提交并推送 ---
echo "🚀 正在提交：$COMMIT_MSG"
git commit -m "$COMMIT_MSG"
git push origin "$BRANCH_NAME"
echo "✅ 成功推送改动到 $REPO_URL 的 $BRANCH_NAME 分支。"
