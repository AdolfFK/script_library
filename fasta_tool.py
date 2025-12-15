#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名：fasta_tool.py
功能: 该脚本是一个功能全面的FASTA序列处理工具，支持按ID提取/剔除、长度过滤、去重、格式标准化、序列变换（反向互补、大小写、子序列）、碱基统计及合法性验证等
作者: adolf
邮箱: wangjyafk@126.com
创建: 2025-12-10
更新: 2025-12-14  adolf  初版：整合字母/数字/长度+正反向
依赖: Python ≥ 3.6（仅标准库）
用法:
    # 提取长度在 100~500 之间的序列，并去重
    python fasta_tool.py -i input.fasta --min-len 100 --max-len 500 --dedup-seq -o filtered.fasta

    # 反向互补所有序列并转为大写
    python fasta_tool.py -i dna.fasta --reverse-complement --uppercase -o rc_upper.fasta

    # 检查FASTA是否全是合法DNA
    python fasta_tool.py -i test.fasta --validate-dna

    # 提取每个序列的第 10 到 50 位
    python fasta_tool.py -i genes.fasta --subseq 10 50 -o subseqs.fasta

    # 统计碱基组成
    python fasta_tool.py -i genome.fasta --base-composition

仓库: https://github.com/AdolfFK/script_library
==============================================================================
"""


import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict

def read_id_file(filepath):
    ids = set()
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                ids.add(line)
    return ids

def write_fasta(records, output_handle):
    SeqIO.write(records, output_handle, 'fasta')

def validate_dna_seq(seq_str):
    valid_bases = set('ATCGatcgNn')
    return set(seq_str) <= valid_bases

def main():
    parser = argparse.ArgumentParser(
        description="增强版 FASTA 序列处理工具",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('--input', '-i', required=True,
                        help='输入的FASTA文件路径')
    parser.add_argument('--output', '-o',
                        help='输出FASTA文件路径（若未指定则输出到stdout）')
    parser.add_argument('--remove-ids',
                        help='剔除该文件中列出的ID（每行一个）')
    parser.add_argument('--extract-ids',
                        help='仅提取该文件中列出的ID对应的序列')
    parser.add_argument('--list-ids', action='store_true',
                        help='列出所有序列ID')
    parser.add_argument('--stats', action='store_true',
                        help='输出每个序列的ID和长度')
    parser.add_argument('--rename-prefix', help='为ID添加前缀')
    parser.add_argument('--rename-suffix', help='为ID添加后缀')

    parser.add_argument('--min-len', type=int,
                        help='保留长度 >= 此值的序列')
    parser.add_argument('--max-len', type=int,
                        help='保留长度 <= 此值的序列')
    parser.add_argument('--dedup-seq', action='store_true',
                        help='基于序列内容去重（保留首次出现）')
    parser.add_argument('--dedup-id', action='store_true',
                        help='基于ID去重（保留首次出现）')
    parser.add_argument('--to-single-line', action='store_true',
                        help='将多行FASTA序列转为单行（标准格式）')
    parser.add_argument('--reverse-complement', action='store_true',
                        help='对DNA序列进行反向互补（仅核苷酸序列）')
    parser.add_argument('--uppercase', action='store_true',
                        help='将序列转为大写')
    parser.add_argument('--lowercase', action='store_true',
                        help='将序列转为小写')
    parser.add_argument('--subseq', nargs=2, metavar=('START', 'END'), type=int,
                        help='提取子序列（1-based，闭区间，如 10 20）')
    parser.add_argument('--base-composition', action='store_true',
                        help='统计每个序列的碱基组成（A/T/C/G/N等）')
    parser.add_argument('--validate-dna', action='store_true',
                        help='检查序列是否只包含合法DNA字符（ATCGN）')

    args = parser.parse_args()

    # 读取输入
    try:
        records = list(SeqIO.parse(args.input, 'fasta'))
    except Exception as e:
        sys.exit(f"错误：无法读取FASTA文件 {args.input} - {e}")

    if not records:
        sys.exit(f"错误：{args.input} 中无有效FASTA序列")

    # 功能：列出ID
    if args.list_ids:
        for rec in records:
            print(rec.id)
        return

    # 功能：碱基组成统计
    if args.base_composition:
        for rec in records:
            seq_str = str(rec.seq).upper()
            counts = {base: seq_str.count(base) for base in "ATCGN"}
            total = len(seq_str)
            others = total - sum(counts.values())
            if others > 0:
                counts['Other'] = others
            comp_str = ' '.join(f"{k}:{v}" for k, v in counts.items() if v > 0)
            print(f"{rec.id}\tLen:{total}\t{comp_str}")
        return

    # 功能：验证DNA合法性
    if args.validate_dna:
        invalid = []
        for rec in records:
            if not validate_dna_seq(str(rec.seq)):
                invalid.append(rec.id)
        if invalid:
            print("以下序列包含非法字符:", file=sys.stderr)
            for rid in invalid:
                print(rid, file=sys.stderr)
            sys.exit(1)
        else:
            print("所有序列均为合法DNA字符", file=sys.stderr)
            return

    # Step 1: ID 过滤（extract / remove）
    keep_set = read_id_file(args.extract_ids) if args.extract_ids else None
    remove_set = read_id_file(args.remove_ids) if args.remove_ids else None

    filtered = []
    for rec in records:
        if keep_set is not None and rec.id not in keep_set:
            continue
        if remove_set is not None and rec.id in remove_set:
            continue
        filtered.append(rec)

    # Step 2: 长度过滤
    if args.min_len is not None or args.max_len is not None:
        new_filtered = []
        for rec in filtered:
            l = len(rec.seq)
            if args.min_len is not None and l < args.min_len:
                continue
            if args.max_len is not None and l > args.max_len:
                continue
            new_filtered.append(rec)
        filtered = new_filtered

    # Step 3: 去重
    if args.dedup_id:
        seen = set()
        uniq = []
        for rec in filtered:
            if rec.id not in seen:
                seen.add(rec.id)
                uniq.append(rec)
        filtered = uniq

    if args.dedup_seq:
        seen_seq = set()
        uniq = []
        for rec in filtered:
            s = str(rec.seq)
            if s not in seen_seq:
                seen_seq.add(s)
                uniq.append(rec)
        filtered = uniq

    # Step 4: 序列变换（注意顺序）
    transformed = []
    for rec in filtered:
        seq_str = str(rec.seq)

        # 子序列提取（1-based inclusive）
        if args.subseq:
            start, end = args.subseq
            if start < 1:
                start = 1
            seq_str = seq_str[start-1:end]  # Python 是 0-based, end exclusive

        # 反向互补
        if args.reverse_complement:
            seq_obj = Seq(seq_str)
            seq_str = str(seq_obj.reverse_complement())

        # 大小写转换
        if args.uppercase:
            seq_str = seq_str.upper()
        elif args.lowercase:
            seq_str = seq_str.lower()

        # 更新序列
        rec.seq = Seq(seq_str)

        # 重命名ID
        if args.rename_prefix or args.rename_suffix:
            prefix = args.rename_prefix or ""
            suffix = args.rename_suffix or ""
            rec.id = prefix + rec.id + suffix
            rec.description = ""

        transformed.append(rec)

    # Step 5: 输出
    output_handle = open(args.output, 'w') if args.output else sys.stdout

    if args.to_single_line:
        # 手动写单行FASTA（避免Biopython默认折行）
        for rec in transformed:
            output_handle.write(f">{rec.id}\n{str(rec.seq)}\n")
    else:
        write_fasta(transformed, output_handle)

    if args.output:
        output_handle.close()
        print(f"处理完成，结果已保存至: {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
