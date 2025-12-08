#!/usr/bin/env bash
# 提交并推送本地 Git 仓库目录 ($LOCAL_DIR) 中的改动。
# 支持两种模式：整体同步 (1个参数) 和 单文件/目录提交 (2个参数)。

# 严格模式
set -euo pipefail

# --- 配置区 ---
REPO_URL="https://github.com/AdolfFK/script_library" # 仅用于首次克隆
LOCAL_DIR="/Users/adolf1/Documents/Work/Script_Library"
BRANCH_NAME="main"

# --- 模式和输入检查 ---

MODE=""
TARGET_FULL_PATH="" # 用于存储目标文件/目录的绝对路径
COMMIT_MSG=""

if [[ $# -eq 1 ]]; then
    # 模式 1: 整体同步 (1 个参数)
    MODE="BULK"
    COMMIT_MSG="$1"
    echo "💡 当前模式: 整体同步所有文件。"

elif [[ $# -eq 2 ]]; then
    # 模式 2: 单文件/目录提交 (2 个参数)
    MODE="SINGLE_FILE"
    TARGET_PATH_ARG="$1" # 用户输入的路径
    COMMIT_MSG="$2"
    
    # 1. 检查目标路径是否存在 (相对于脚本执行目录)
    if [[ ! -e "$TARGET_PATH_ARG" ]]; then
        echo "❌ 错误: 指定的文件或目录不存在: $TARGET_PATH_ARG"
        exit 1
    fi
    
    # 2. 将目标路径解析为绝对路径。
    # 这样在脚本切换到 $LOCAL_DIR 后，我们仍能准确引用它。
    # 注意: realpath 命令在大多数 Unix/Linux/Mac 系统上都可用。
    TARGET_FULL_PATH=$(realpath "$TARGET_PATH_ARG")
    
    echo "💡 当前模式: 单文件/目录提交 -> $TARGET_FULL_PATH"

else
    # 错误用法
    echo "❌ 用法错误: 必须提供 1 个参数 (整体同步) 或 2 个参数 (单文件/目录提交)。"
    echo ""
    echo "用法 1 (整体同步): $0 \"<整体提交信息>\""
    echo "用法 2 (单文件/目录提交): $0 <文件/目录路径> \"<文件提交信息>\""
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
    # 整体同步模式：暂存所有
    git add -A .
elif [[ "$MODE" == "SINGLE_FILE" ]]; then
    # 单文件/目录模式：
    # git add -A <目录路径> 会递归地暂存目录及其所有内容（包括新增、修改和删除）。
    git add -A "$TARGET_FULL_PATH"
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
