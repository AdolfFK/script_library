#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名：sort_fasta_by_id.py
功能：  按序列ID对FASTA文件进行排序，并按标准格式输出（每行60字符）
作者：  adolf <wangjyafk@126.com>
创建：  2025-06-01
更新：  2025-12-09  adolf  改用 argparse 长参数，优化输出格式与健壮性
依赖：  Python ≥ 3.6（仅标准库）
用法：
    python sort_fasta_by_id.py --input input.fasta --output sorted.fasta
    python sort_fasta_by_id.py -i input.fasta -o sorted.fasta
    
    # 基本用法（默认每行60字符）
    python sort_fasta_by_id.py --input input.fasta --output sorted.fasta

    # 不换行输出序列
    python sort_fasta_by_id.py -i input.fasta -o sorted.fasta --no-wrap

    # 自定义每行80字符
    python sort_fasta_by_id.py -i input.fasta -o sorted.fasta --wrap-length 80
仓库：  https://github.com/AdolfFK/script_library
==============================================================================
"""

import argparse
import os
import sys


def read_fasta(fasta_file):
    """
    读取FASTA文件，返回字典：{sequence_id: sequence}
    """
    if not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"Input file not found: {fasta_file}")

    fasta_dict = {}
    current_id = None

    with open(fasta_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:].split()[0]  # 只取ID部分（空格前）
                if current_id in fasta_dict:
                    print(f"Warning: Duplicate ID '{current_id}' found. Overwriting.", file=sys.stderr)
                fasta_dict[current_id] = ""
            elif current_id is not None:
                fasta_dict[current_id] += line
            else:
                # 序列出现在任何ID之前，视为格式错误
                raise ValueError("FASTA format error: sequence before header.")

    return fasta_dict


def write_fasta_record(outfile, seq_id, sequence, wrap=60):
    """
    写入单条FASTA记录，序列按 wrap 字符每行换行（wrap=None 则不换行）
    """
    outfile.write(f">{seq_id}\n")
    if wrap and wrap > 0:
        for i in range(0, len(sequence), wrap):
            outfile.write(sequence[i:i + wrap] + "\n")
    else:
        outfile.write(sequence + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Sort sequences in a FASTA file by their IDs (lexicographically)."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file path."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output sorted FASTA file path."
    )
    parser.add_argument(
        "--no-wrap",
        action="store_true",
        help="Do not wrap sequences (write entire sequence on one line)."
    )
    parser.add_argument(
        "--wrap-length",
        type=int,
        default=60,
        help="Number of characters per line for sequence (default: 60)."
    )

    args = parser.parse_args()

    try:
        fasta_dict = read_fasta(args.input)
        wrap_len = None if args.no_wrap else args.wrap_length

        with open(args.output, "w") as outfile:
            for seq_id in sorted(fasta_dict.keys()):
                write_fasta_record(outfile, seq_id, fasta_dict[seq_id], wrap=wrap_len)

        print(f"Successfully sorted {len(fasta_dict)} sequences to '{args.output}'.", file=sys.stderr)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
