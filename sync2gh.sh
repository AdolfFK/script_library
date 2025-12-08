#!/usr/bin/env bash
# 提交并推送本地 Git 仓库目录 ($LOCAL_DIR) 中所有已有的改动。
# 不涉及从外部复制文件。
#
# 用法:
#   ./sync_local.sh "我的脚本库日常更新"

# 严格模式
set -euo pipefail

# --- 配置区 ---
REPO_URL="https://github.com/AdolfFK/script_library" # 仅用于首次克隆
LOCAL_DIR="/Users/adolf1/Documents/Work/Script_Library"
BRANCH_NAME="main" 

# --- 输入检查 ---
if [[ $# -ne 1 ]]; then
    echo "❌ 用法错误: 必须提供一个 提交信息。"
    echo "用法: $0 \"<提交信息>\""
    exit 1
fi

COMMIT_MSG="$1"

# 1. 确保本地仓库目录存在，并切换到该目录
mkdir -p "$LOCAL_DIR"
cd "$LOCAL_DIR"

# 2. 若本地仓库不存在，则克隆
if [[ ! -d ".git" ]]; then
    echo "📥 首次运行，正在克隆仓库到 $LOCAL_DIR ..."
    # 使用 --depth 1 优化克隆速度
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

# 4. 提交并推送所有改动
# 暂存所有（新增、修改、删除）
echo "📦 正在暂存所有改动..."
git add -A .

# 检查是否有实际的暂存改动
if git diff --staged --quiet; then
    echo "ℹ️  无改动，无需提交。"
    exit 0
fi

echo "🚀 正在提交和推送改动..."
git commit -m "$COMMIT_MSG"
git push origin "$BRANCH_NAME"

echo "✅ 成功推送所有改动到 $REPO_URL 的 $BRANCH_NAME 分支"