#!/usr/bin/env bash
#==============================================================================
# 脚本名： map_generator.sh
# 功能：  递归扫描子目录，生成 SCRIPT_MAP.md（含目录列）
# 作者：  adolf1
# 邮箱:   wangjyafk@126.com
# 创建：  2025-06-01
# 更新：  2025-12-14  adolf1  新增“目录”列
# 依赖：  bash≥4.0  curl≥7.0
# 用法：  
#        bash map_generator.sh
#
# 仓库：  https://github.com/AdolfFK/script_library
#==============================================================================

set -euo pipefail

# 1. 写表头（新增“目录”列）
cat > SCRIPT_MAP.md <<'EOF'
| 目录 | 脚本 | 功能 | 更新日期 |
|------|------|------|----------|
EOF

# 2. find 递归找文件
while IFS= read -r -d '' f; do
  # 去 BOM
  tail -c +4 "$f" 2>/dev/null > /tmp/$$.tmp || cp "$f" /tmp/$$.tmp
  
  # 提取脚本名（支持多行注释）
  name=$(sed -n -E '1,20s/^[[:space:]]*#?[[:space:]]*脚本名[:：][[:space:]]*//p' /tmp/$$.tmp | head -n1 | tr -d '\n\r')
  
  # 提取功能（处理多行情况）
  func=$(sed -n -E '1,20s/^[[:space:]]*#?[[:space:]]*功能[:：][[:space:]]*//p' /tmp/$$.tmp | head -n1 | tr -d '\n\r')
  
  # 如果功能描述为空或太短，尝试从脚本注释中提取更多内容
  if [[ -z "$func" || ${#func} -lt 5 ]]; then
    # 尝试从Python多行注释中提取功能描述
    if [[ "$f" == *.py ]]; then
      func=$(grep -A 1 "功能" /tmp/$$.tmp | tail -n1 | sed 's/^[[:space:]]*//' | tr -d '\n\r')
    fi
  fi
  
  # 提取更新信息
  upd=$(sed -n -E '1,20s/^[[:space:]]*#?[[:space:]]*更新[:：][[:space:]]*//p' /tmp/$$.tmp | head -n1 | tr -d '\n\r')
  date=${upd%% *}

  # 获取脚本所在目录，并标准化 . → ./
  dir=$(dirname "$f")
  [[ "$dir" == "." ]] && dir="./"

  [[ -z "$name" ]] && name="$(basename "$f")"
  [[ -z "$func" ]] && func="—"
  [[ -z "$date" ]] && date="—"
  
  # 清理功能描述中的额外空格和换行
  func=$(echo "$func" | sed 's/[[:space:]]+/ /g' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

  # 输出：目录 | 脚本名 | 功能 | 日期
  printf "| %s | %s | %s | %s |\n" "$dir" "$name" "$func" "$date" >> SCRIPT_MAP.md
done < <(find . -type f \( -name "*.py" -o -name "*.pl" -o -name "*.R" -o -name "*.sh" \) -print0)

rm -f /tmp/$$.tmp
echo "✅ SCRIPT_MAP.md 已生成（含子目录和目录列）"
