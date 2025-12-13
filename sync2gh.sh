#!/usr/bin/env bash
#==============================================================================
# 脚本名: sync2gh.sh
# 功能:  极简版——永远向 main 分支提交并推送；版本管理完全由用户自行处理。
#        本地目录为权威源，脚本仅负责向上同步，绝不从远程拉取任何内容。
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-13  adolf1  重写逻辑：本地为主，不 clone、不 pull，彻底修复变量未定义问题
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
if [[ $# -eq 1 ]]; then
  MODE="BULK"
  COMMIT_MSG="$1"
  echo "💡 模式: 整体同步。"
elif [[ $# -eq 2 ]]; then
  MODE="SINGLE_FILE"
  TARGET_PATH_ARG="$1"
  COMMIT_MSG="$2"
  # 检查用户指定的路径是否存在
  if [[ ! -e "$TARGET_PATH_ARG" ]]; then
    echo "❌ 指定路径不存在: $TARGET_PATH_ARG"
    exit 1
  fi
  echo "💡 模式: 单文件/目录提交。"
else
  echo "❌ 用法错误！"
  echo "整体同步: $0 '<提交信息>'"
  echo "单路径提交: $0 <文件/目录> '<提交信息>'"
  exit 1
fi

# ---------------- 进入本地工作目录 ----------------
# 确保目录存在（但不干扰已有内容）
mkdir -p "$LOCAL_DIR"
cd "$LOCAL_DIR" || { echo "❌ 无法进入目录: $LOCAL_DIR"; exit 1; }

# ---------------- 初始化为 Git 仓库（如尚未初始化）----------------
if [[ ! -d .git ]]; then
  echo "🔧 当前目录不是 Git 仓库，正在初始化..."
  git init -q
  git remote add origin "$REPO_URL"
  # 设置默认分支名为 main（兼容 Git 2.28+ 及旧版本）
  git symbolic-ref HEAD refs/heads/"$BRANCH_NAME"
  echo "✅ 已初始化本地 Git 仓库并关联远程: $REPO_URL"
fi

# ---------------- 暂存改动 ----------------
echo "📦 正在暂存改动..."
if [[ "$MODE" == "BULK" ]]; then
  # 整体同步：添加所有变更
  git add -A .
else
  # 单文件/目录模式：解析绝对路径
  if command -v realpath >/dev/null 2>&1; then
    TARGET_FULL_PATH=$(realpath "$TARGET_PATH_ARG")
  elif command -v python3 >/dev/null 2>&1; then
    TARGET_FULL_PATH=$(python3 -c "import os, sys; print(os.path.realpath(sys.argv[1]))" "$TARGET_PATH_ARG")
  else
    # 最后手段：使用原始路径（可能是相对路径）
    TARGET_FULL_PATH="$TARGET_PATH_ARG"
  fi

  # 安全检查：确保目标路径位于当前仓库目录内
  # 使用字符串前缀匹配判断是否在当前目录下
  case "$TARGET_FULL_PATH" in
    "$(pwd)"/*|"$PWD") 
      # 路径合法
      ;;
    "$PWD") 
      # 路径正好是仓库根目录
      ;;
    *)
      echo "❌ 错误：指定路径不在当前仓库目录内！"
      echo "   仓库根目录: $PWD"
      echo "   指定路径:   $TARGET_FULL_PATH"
      exit 1
      ;;
  esac

  # 添加指定文件或目录
  git add -A "$TARGET_FULL_PATH"
fi

# ---------------- 检查是否有实际改动 ----------------
if git diff --staged --quiet; then
  echo "ℹ️  无改动，无需提交。"
  exit 0
fi

# ---------------- 提交 ----------------
echo "🚀 正在提交：$COMMIT_MSG"
git commit -m "$COMMIT_MSG"

# ---------------- 推送至远程仓库 ----------------
echo "📤 正在推送到远程仓库 ($REPO_URL) 的 $BRANCH_NAME 分支..."
# 尝试带 lease 的强制推送（安全），若失败则尝试首次推送（设置 upstream）
if ! git push --force-with-lease origin "$BRANCH_NAME" 2>/dev/null; then
  # 可能是首次推送，没有 upstream
  git push -u origin "$BRANCH_NAME"
fi

echo "✅ 成功推送改动到 $REPO_URL 的 $BRANCH_NAME 分支。"

