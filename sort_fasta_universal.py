#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名：sort_fasta_universal.py
功能: 对 FASTA 文件多维度排序（ID 字母、ID 数字、序列长度、正/反序）
作者: adolf
邮箱: wangjyafk@126.com
创建: 2025-12-10
更新: 2025-12-10  adolf  初版：整合字母/数字/长度+正反向
依赖: Python ≥ 3.6（仅标准库）
用法:
    # 按 ID 字母正向（默认）
    python sort_fasta_universal.py -i in.fasta -o out.fasta

    # 按 ID 字母反向
    python sort_fasta_universal.py -i in.fasta -o out.fasta --mode id_alpha --reverse

    # 按 ID 中的数字大小正向
    python sort_fasta_universal.py -i in.fasta -o out.fasta --mode id_number

    # 按序列长度反向（最长在前）
    python sort_fasta_universal.py -i in.fasta -o out.fasta --mode seqlen --reverse

    # 不换行输出
    python sort_fasta_universal.py -i in.fasta -o out.fasta --no-wrap

仓库: https://github.com/AdolfFK/script_library
==============================================================================
"""

import argparse
import os
import sys
import re
from typing import Dict, List, Tuple


# ---------- 通用子函数 ----------
def read_fasta(fasta_file: str) -> Dict[str, str]:
    """读取 FASTA 文件，返回 {seq_id: sequence} 字典"""
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"Input file not found: {fasta_file}")

    fasta_dict: Dict[str, str] = {}
    current_id: str | None = None

    with open(fasta_file, "r", encoding="utf-8") as inh:
        for line in inh:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                if current_id in fasta_dict:
                    print(f"Warning: Duplicate ID '{current_id}' found. Overwriting.",
                          file=sys.stderr)
                fasta_dict[current_id] = ""
            elif current_id:
                fasta_dict[current_id] += line
            else:
                raise ValueError("FASTA format error: sequence before header.")
    return fasta_dict


def write_fasta_record(outfile, seq_id: str, sequence: str, wrap: int | None = 60):
    """写入单条 FASTA 记录，支持换行控制"""
    outfile.write(f">{seq_id}\n")
    if wrap and wrap > 0:
        for i in range(0, len(sequence), wrap):
            outfile.write(sequence[i:i + wrap] + "\n")
    else:
        outfile.write(sequence + "\n")


# ---------- 排序关键函数 ----------
def _id_number_key(seq_id: str) -> int:
    """提取 ID 中第一组连续数字作为排序键，无数字则返回 0"""
    nums = re.findall(r"\d+", seq_id)
    return int(nums[0]) if nums else 0


def sort_items(items: List[Tuple[str, str]], mode: str, reverse: bool = False):
    """
    对 [(seq_id, sequence), ...] 列表进行指定规则排序
    mode: id_alpha / id_number / seqlen
    """
    if mode == "id_alpha":
        return sorted(items, key=lambda x: x[0], reverse=reverse)
    elif mode == "id_number":
        return sorted(items, key=lambda x: _id_number_key(x[0]), reverse=reverse)
    elif mode == "seqlen":
        return sorted(items, key=lambda x: len(x[1]), reverse=reverse)
    else:
        raise ValueError(f"Unknown sort mode: {mode}")


# ---------- 主入口 ----------
def main():
    parser = argparse.ArgumentParser(
        description="Sort FASTA entries by ID (alphabetical/numerical) or sequence length."
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA file path.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output sorted FASTA file path.")
    parser.add_argument("-m", "--mode", default="id_alpha",
                        choices=["id_alpha", "id_number", "seqlen"],
                        help="Sort mode: id_alpha=ID字母(default), id_number=ID数字, seqlen=序列长度")
    parser.add_argument("-r", "--reverse", action="store_true",
                        help="Reverse sort order.")
    parser.add_argument("--no-wrap", action="store_true",
                        help="Do not wrap sequences (write entire sequence on one line).")
    parser.add_argument("--wrap-length", type=int, default=60,
                        help="Number of characters per line for sequence (default: 60).")

    args = parser.parse_args()

    try:
        fasta_dict = read_fasta(args.input)
        items = list(fasta_dict.items())          # [(id, seq), ...]
        sorted_items = sort_items(items, args.mode, args.reverse)

        wrap_len = None if args.no_wrap else args.wrap_length
        with open(args.output, "w", encoding="utf-8") as outh:
            for seq_id, seq in sorted_items:
                write_fasta_record(outh, seq_id, seq, wrap=wrap_len)

        print(f"Successfully sorted {len(sorted_items)} sequences "
              f"({args.mode}, reverse={args.reverse}) to '{args.output}'.",
              file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
