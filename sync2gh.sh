#!/usr/bin/env bash
#==============================================================================
# 脚本名: sync2gh.sh
# 功能: 将本地已存在的工作目录同步到 GitHub 的 main 分支
#      - 本地是权威源，绝不从远程拉取覆盖
#      - 若本地尚未初始化为 Git 仓库，则自动初始化并关联远程
# 作者: adolf1
# 邮箱: wangjyafk@126.com
# 创建: 2025-06-01
# 更新: 2025-12-13  按用户需求重写：本地为主，不 clone、不 pull
#==============================================================================

set -euo pipefail

# ---------------- 配置区 ----------------
REPO_URL="https://github.com/AdolfFK/script_library"
LOCAL_DIR="${LOCAL_DIR:-/Users/adolf1/Documents/Work/Script_Library}"
BRANCH_NAME="main"

# ---------------- 工具检测 ----------------
command -v git >/dev/null 2>&1 || {
  echo "❌ 请先安装 Git 并确保其在 PATH 中。"
  exit 1
}

# ---------------- 参数检查 ----------------
MODE=""
TARGET_FULL_PATH=""
COMMIT_MSG=""

if [[ $# -eq 1 ]]; then
  MODE="BULK"; COMMIT_MSG="$1"
  echo "💡 模式: 整体同步本地目录。"
elif [[ $# -eq 2 ]]; then
  MODE="SINGLE_FILE"
  TARGET_PATH_ARG="$1"; COMMIT_MSG="$2"
  [[ -e "$TARGET_PATH_ARG" ]] || { echo "❌ 指定路径不存在: $TARGET_PATH_ARG"; exit 1; }
  if command -v realpath >/dev/null 2>&1; then
    TARGET_FULL_PATH=$(realpath "$TARGET_PATH_ARG")
  elif command -v python3 >/dev/null 2>&1; then
    TARGET_FULL_PATH=$(python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$TARGET_PATH_ARG")
  else
    TARGET_FULL_PATH="$TARGET_PATH_ARG"
  fi
  echo "💡 模式: 单文件/目录提交 -> $TARGET_FULL_PATH。"
else
  echo "❌ 用法错误！"
  echo "整体同步: $0 '<提交信息>'"
  echo "单路径提交: $0 <文件/目录> '<提交信息>'"
  exit 1
fi

# ---------------- 进入本地目录 ----------------
mkdir -p "$LOCAL_DIR"  # 确保目录存在（但不干扰已有内容）
cd "$LOCAL_DIR" || exit 1

# ---------------- 初始化为 Git 仓库（如需要）----------------
if [[ ! -d .git ]]; then
  echo "🔧 当前目录不是 Git 仓库，正在初始化..."
  git init -q
  git remote add origin "$REPO_URL"
  # 设置默认分支名为 main（兼容新旧 Git）
  git symbolic-ref HEAD refs/heads/"$BRANCH_NAME"
fi

# ---------------- 暂存 & 提交 ----------------
echo "📦 正在暂存改动..."
if [[ "$MODE" == "BULK" ]]; then
  git add -A .
else
  # 确保目标路径在当前目录内（安全检查）
  case "$TARGET_FULL_PATH" in
    "$(pwd)"/*|"$PWD") ;;
    *) echo "❌ 错误：路径不在当前仓库内：$TARGET_FULL_PATH"; exit 1 ;;
  esac
  git add -A "$TARGET_FULL_PATH"
fi

if git diff --staged --quiet; then
  echo "ℹ️  无改动，无需提交。"
  exit 0
fi

echo "🚀 正在提交：$COMMIT_MSG"
git commit -m "$COMMIT_MSG"

# ---------------- 推送 ----------------
echo "📤 正在推送到远程仓库..."
# 首次推送需设置上游
git push --force-with-lease origin "$BRANCH_NAME" || \
git push -u origin "$BRANCH_NAME"

echo "✅ 成功推送至 $REPO_URL 的 $BRANCH_NAME 分支。"
