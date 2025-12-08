#!/usr/bin/env bash
# 提交并推送本地 Git 仓库目录 ($LOCAL_DIR) 中的改动。
# 支持两种模式：整体同步 (1个参数) 和 单文件提交 (2个参数)。

# 严格模式
set -euo pipefail

# --- 配置区 ---
REPO_URL="https://github.com/AdolfFK/script_library" # 仅用于首次克隆
LOCAL_DIR="/Users/adolf1/Documents/Work/Script_Library"
BRANCH_NAME="main"

# --- 模式和输入检查 ---

MODE=""
TARGET_PATH=""
COMMIT_MSG=""

if [[ $# -eq 1 ]]; then
    # 模式 1: 整体同步 (1 个参数)
    MODE="BULK"
    COMMIT_MSG="$1"
    echo "💡 当前模式: 整体同步所有文件。"

elif [[ $# -eq 2 ]]; then
    # 模式 2: 单文件提交 (2 个参数)
    MODE="SINGLE_FILE"
    TARGET_PATH="$1"
    COMMIT_MSG="$2"
    
    # 检查目标路径是否存在
    if [[ ! -e "$TARGET_PATH" ]]; then
        echo "❌ 错误: 指定的文件或目录不存在: $TARGET_PATH"
        exit 1
    fi
    echo "💡 当前模式: 单文件/目录提交 -> $TARGET_PATH"

else
    # 错误用法
    echo "❌ 用法错误: 必须提供 1 个参数 (整体同步) 或 2 个参数 (单文件提交)。"
    echo ""
    echo "用法 1 (整体同步): $0 \"<整体提交信息>\""
    echo "用法 2 (单文件提交): $0 <文件/目录路径> \"<文件提交信息>\""
    exit 1
fi

# 1. 确保本地仓库目录存在，并切换到该目录
mkdir -p "$LOCAL_DIR"
cd "$LOCAL_DIR"

# 2. 若本地仓库不存在，则克隆
if [[ ! -d ".git" ]]; then
    echo "📥 首次运行，正在克隆仓库到 $LOCAL_DIR ..."
    git clone --depth 1 "$REPO_URL" .
    
    if [[ $? -ne 0 ]]; then
        echo "❌ 克隆仓库失败！请检查配置或网络连接。"
        exit 1
    fi
    
    git checkout "$BRANCH_NAME" || { echo "❌ 切换到 $BRANCH_NAME 分支失败"; exit 1; }
fi

# 3. 拉取远程最新内容（避免冲突）
echo "🔄 正在拉取远程最新内容..."
git pull origin "$BRANCH_NAME"

# 4. 暂存改动
echo "📦 正在暂存改动..."
if [[ "$MODE" == "BULK" ]]; then
    # 暂存所有（新增、修改、删除）
    git add -A .
elif [[ "$MODE" == "SINGLE_FILE" ]]; then
    # 仅暂存指定的文件/目录
    # 注意: TARGET_PATH 必须是绝对路径，或相对于脚本运行目录的路径
    git add -A "$TARGET_PATH"
fi

# 5. 检查是否有实际的暂存改动
if git diff --staged --quiet; then
    echo "ℹ️  无改动，无需提交。"
    exit 0
fi

# 6. 提交并推送
echo "🚀 正在提交改动，提交信息: \"$COMMIT_MSG\""
git commit -m "$COMMIT_MSG"
git push origin "$BRANCH_NAME"

echo "✅ 成功推送所有改动到 $REPO_URL 的 $BRANCH_NAME 分支"
