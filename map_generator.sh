#!/usr/bin/env bash
# 自动生成 SCRIPT_MAP.md  —— macOS / Linux 通用
set -euo pipefail

echo "| 脚本 | 功能 | 更新日期 |"
echo "|------|------|----------|"

for f in *.py *.pl *.R *.sh; do
  [[ -f "$f" ]] || continue

  # 用 grep + sed 提取，避开了 awk 高级语法
  name=$(grep -m1 -E "^脚本名：" "$f" | sed 's/^脚本名： *//')
  func=$(grep -m1 -E "^功能："  "$f" | sed 's/^功能： *//')
  upd=$(grep  -m1 -E "^更新："  "$f" | sed 's/^更新： *//')
  date=${upd%% *}                  # 取第一个空格前的日期

  # 如果某行缺失，给默认值
  [[ -z "$name" ]] && name="$f"
  [[ -z "$func" ]] && func="—"
  [[ -z "$date" ]] && date="—"

  printf "| [%s](%s) | %s | %s |\n" "$name" "$f" "$func" "$date"
done
