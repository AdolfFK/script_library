#!/usr/bin/env bash
#==============================================================================
# 脚本名: sync2gh.sh
# 功能:  极简版——永远向 main 分支提交并推送；版本管理完全由用户自行处理。
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-13  修复变量未定义与克隆目录问题
# 依赖:  bash≥3.2  git≥2.0
# 用法:
#        ./sync2gh.sh "日常提交说明"
#        ./sync2gh.sh script.py "修复 typo"
# 仓库:  https://github.com/AdolfFK/script_library
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
TARGET_FULL_PATH=""   # ←←← 关键：显式初始化，避免 set -u 报错
COMMIT_MSG=""

if [[ $# -eq 1 ]]; then
  MODE="BULK"
  COMMIT_MSG="$1"
  echo "💡 模式: 整体同步。"
elif [[ $# -eq 2 ]]; then
  MODE="SINGLE_FILE"
  TARGET_PATH_ARG="$1"
  COMMIT_MSG="$2"
  [[ -e "$TARGET_PATH_ARG" ]] || { echo "❌ 指定路径不存在: $TARGET_PATH_ARG"; exit 1; }

  # 跨平台获取绝对路径
  if command -v realpath >/dev/null 2>&1; then
    TARGET_FULL_PATH=$(realpath "$TARGET_PATH_ARG")
  elif command -v python3 >/dev/null 2>&1; then
    TARGET_FULL_PATH=$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$TARGET_PATH_ARG")
  else
    # 最后手段：使用原始路径（可能是相对路径）
    TARGET_FULL_PATH="$TARGET_PATH_ARG"
  fi
  echo "💡 模式: 单文件/目录提交 -> $TARGET_FULL_PATH。"
else
  echo "❌ 用法错误！"
  echo "整体同步: $0 '<提交信息>'"
  echo "单路径提交: $0 <文件/目录> '<提交信息>'"
  exit 1
fi

# ---------------- 仓库初始化 ----------------
mkdir -p "$LOCAL_DIR"
cd "$LOCAL_DIR" || exit 1

if [[ ! -d .git ]]; then
  echo "📥 首次运行，正在克隆仓库..."
  # 检查当前目录是否为空（避免 git clone . 失败）
  if [[ -n "$(ls -A . 2>/dev/null)" ]]; then
    echo "❌ 错误：\$LOCAL_DIR 非空且不是 Git 仓库。"
    echo "   请清空 $LOCAL_DIR 或删除后重试。"
    exit 1
  fi
  git clone --depth 1 "$REPO_URL" .
fi

echo "🔄 拉取远程最新内容..."
if ! git pull origin "$BRANCH_NAME"; then
  echo "⚠️ 警告：拉取远程更新失败（可能有冲突），但仍继续提交。"
fi

# ---------------- 暂存 & 提交 ----------------
echo "📦 正在暂存改动..."
if [[ "$MODE" == "BULK" ]]; then
  git add -A .
else
  # 确保目标路径在当前工作目录内（防止添加外部文件）
  # 将路径标准化为绝对路径进行比较
  REPO_ROOT=$(pwd)
  # 使用 dirname + basename 无法可靠判断，改用字符串前缀匹配（需确保 TARGET_FULL_PATH 是绝对路径）
  if [[ "$TARGET_FULL_PATH" != "$REPO_ROOT"* ]]; then
    echo "❌ 错误：指定路径不在仓库目录内：$TARGET_FULL_PATH"
    echo "   仓库根目录：$REPO_ROOT"
    exit 1
  fi
  git add -A "$TARGET_FULL_PATH"
fi

# 检查是否有实际改动
if git diff --staged --quiet; then
  echo "ℹ️  无改动，无需提交。"
  exit 0
fi

echo "🚀 正在提交：$COMMIT_MSG"
git commit -m "$COMMIT_MSG"
git push origin "$BRANCH_NAME"
echo "✅ 成功推送改动到 $REPO_URL 的 $BRANCH_NAME 分支。"

