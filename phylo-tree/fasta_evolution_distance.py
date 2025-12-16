#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: fasta_evolution_distance.py
功能:  蛋白 FASTA 多序列比对 + 双指标相似性（序列一致率 & BLOSUM62 得分）批量计算，并自动找出进化距离最远的一对序列
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-12-15  adolf1  增加head说明
依赖:  Python≥3.8 pandas≥1.5
用法:  python fasta_evolution_distance.py \
       --input demo.fa \
       --output-result report.txt \
       --output-alignment demo.aln.fasta
       
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""

import argparse
import os
import sys
import subprocess
from Bio import AlignIO
from Bio.Align import substitution_matrices


def run_mafft(input_fasta, output_aln):
    cmd = ["mafft", "--auto", "--quiet", input_fasta]
    try:
        with open(output_aln, "w") as out:
            subprocess.run(cmd, stdout=out, check=True)
    except subprocess.CalledProcessError as e:
        print(f"MAFFT failed: {e}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("Error: MAFFT is not installed or not in PATH.", file=sys.stderr)
        sys.exit(1)


def calculate_combined_metrics(alignment_file, output_result):
    alignment = AlignIO.read(alignment_file, "fasta")
    seq_count = len(alignment)
    if seq_count < 2:
        raise ValueError("FASTA 文件中至少需要两条序列。")

    try:
        blosum62 = substitution_matrices.load("BLOSUM62")
    except Exception as e:
        print(f"无法加载 BLOSUM62 矩阵: {e}", file=sys.stderr)
        sys.exit(1)

    aligned_seqs = {record.id: str(record.seq) for record in alignment}
    ids = list(aligned_seqs.keys())
    pair_count = 0

    # 存储每对的结果
    results = []

    total_identity_percent = 0.0
    total_blosum_score = 0.0

    min_identity = 100.0          # 最低相似性百分比 → 最大进化距离
    max_p_distance_pair = ("", "")  # 基于 p-distance 的最远对

    min_blosum_avg = float('inf')
    min_blosum_pair = ("", "")

    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            id1, id2 = ids[i], ids[j]
            seq1 = aligned_seqs[id1]
            seq2 = aligned_seqs[id2]

            identical = 0
            valid_positions = 0
            blosum_scores = []

            for a, b in zip(seq1, seq2):
                if a == '-' and b == '-':
                    continue
                valid_positions += 1

                # 计算 identity
                if a == b:
                    identical += 1

                # 计算 BLOSUM62
                try:
                    score = blosum62[a, b]
                except KeyError:
                    score = 0
                blosum_scores.append(score)

            # 相似性百分比
            identity_percent = (identical / valid_positions * 100) if valid_positions > 0 else 0.0

            # BLOSUM62 平均得分
            avg_blosum = sum(blosum_scores) / len(blosum_scores) if blosum_scores else 0.0

            results.append({
                'id1': id1,
                'id2': id2,
                'identity_percent': identity_percent,
                'avg_blosum': avg_blosum
            })

            total_identity_percent += identity_percent
            total_blosum_score += avg_blosum
            pair_count += 1

            # 找 p-distance 最大（即 identity 最小）的一对 → 进化距离最远
            if identity_percent < min_identity:
                min_identity = identity_percent
                max_p_distance_pair = (id1, id2)

            # 找 BLOSUM62 最低的一对（备选）
            if avg_blosum < min_blosum_avg:
                min_blosum_avg = avg_blosum
                min_blosum_pair = (id1, id2)

    avg_identity_over_pairs = total_identity_percent / pair_count if pair_count > 0 else 0.0
    avg_blosum_over_pairs = total_blosum_score / pair_count if pair_count > 0 else 0.0

    # === 写入结果文件 ===
    with open(output_result, 'w') as f:
        f.write("=== 蛋白序列多维度相似性与进化距离分析报告 ===\n\n")

        f.write("【指标定义】\n")
        f.write("• 相似性百分比 = (相同氨基酸位点数 / 有效比对位点数) × 100%\n")
        f.write("• BLOSUM62 相似性 = 比对后每个有效位点的 BLOSUM62 得分平均值\n")
        f.write("• 进化距离最远的序列对 = 相似性百分比最低的一对（即 p-distance 最大）\n\n")

        f.write(f"输入比对文件: {alignment_file}\n")
        f.write(f"序列总数: {seq_count}\n")
        f.write(f"序列对总数: {pair_count}\n\n")

        # ====== 相似性百分比结果 ======
        f.write("【相似性百分比分析】\n")
        f.write(f"总相似性百分比（所有对之和）: {total_identity_percent:.6f}%\n")
        f.write(f"平均相似性百分比: {avg_identity_over_pairs:.6f}%\n")
        f.write(f"最低相似性百分比（最大进化距离）: {min_identity:.6f}%\n")
        f.write(f"对应序列（进化距离最远）: {max_p_distance_pair[0]} 与 {max_p_distance_pair[1]}\n\n")

        # ====== BLOSUM62 结果 ======
        f.write("【BLOSUM62 相似性分析】\n")
        f.write(f"总 BLOSUM62 平均得分（所有对之和）: {total_blosum_score:.6f}\n")
        f.write(f"平均 BLOSUM62 相似性: {avg_blosum_over_pairs:.6f}\n")
        f.write(f"最低 BLOSUM62 得分: {min_blosum_avg:.6f}\n")
        f.write(f"对应序列（生化差异最大）: {min_blosum_pair[0]} 与 {min_blosum_pair[1]}\n\n")

        # ====== 详细列表 ======
        f.write("【所有序列对详细结果】\n")
        f.write(f"{'Seq1':<25} {'Seq2':<25} {'Identity (%)':<15} {'Avg BLOSUM62':<15}\n")
        f.write("-" * 80 + "\n")
        for r in results:
            f.write(f"{r['id1']:<25} {r['id2']:<25} {r['identity_percent']:>14.4f}% {r['avg_blosum']:>14.4f}\n")

    # === 终端摘要 ===
    print("\n✅ 分析完成！结果已保存至:", output_result)
    print("\n=== 核心摘要 ===")
    print(f"平均相似性百分比       : {avg_identity_over_pairs:.4f}%")
    print(f"总相似性百分比         : {total_identity_percent:.4f}%")
    print(f"进化距离最远的序列对   : {max_p_distance_pair[0]} ↔ {max_p_distance_pair[1]} ({min_identity:.4f}% identity)")
    print(f"平均 BLOSUM62 相似性   : {avg_blosum_over_pairs:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="比对蛋白 FASTA，同时计算相似性百分比和 BLOSUM62 相似性，并找出进化距离最远的序列对。"
    )
    parser.add_argument("--input", required=True, help="输入的蛋白 FASTA 文件")
    parser.add_argument("--output-result", required=True, help="输出结果文件路径")
    parser.add_argument("--output-alignment", default=None, help="可选：保存比对结果")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"错误：输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)

    aln_file = args.output_alignment or (args.input + ".aligned.fasta")

    print("正在运行 MAFFT 比对...")
    run_mafft(args.input, aln_file)
    print(f"比对完成: {aln_file}")

    print("正在计算相似性百分比与 BLOSUM62 相似性...")
    calculate_combined_metrics(aln_file, args.output_result)


if __name__ == "__main__":
    main()
    