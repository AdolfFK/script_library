#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: mutation_counter_v1.py
功能:  批量扫描 TSV 突变列表，计算各突变/组合频率并筛选高频（≥30 %）主要组合，输出三张频率表
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-06-10  adolf1  增加head说明
依赖:  Python≥3.8 pandas≥1.5
用法:  python mutation_counter_v1.py --input_dir ./in_dir --output_dir ./out_dir
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""

import os
import argparse
from collections import defaultdict

def parse_mutations(line: str) -> (str, list):
    parts = line.strip().split('\t')
    if len(parts) < 2:
        return None, []
    strain_id = parts[0]
    raw_mutations = parts[1]
    if ',' in raw_mutations:
        mutations = raw_mutations.split(',')
    else:
        mutations = [m.strip() for m in raw_mutations.split('\t') if m.strip()]
    return strain_id, mutations

def process_file(filepath):
    mutation_counts = defaultdict(int)
    mutation_strains = defaultdict(set)
    combination_counts = defaultdict(int)
    combination_strains = defaultdict(set)
    strain_total = 0
    strains = []

    # 第一次读取：计算所有突变频率并收集毒株数据
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            strain_id, mutations = parse_mutations(line)
            if not mutations or not strain_id:
                continue
            strains.append((strain_id, mutations))
            strain_total += 1
            for mut in mutations:
                mutation_counts[mut] += 1
                mutation_strains[mut].add(strain_id)
            combo_key = '-'.join(sorted(mutations))
            combination_counts[combo_key] += 1
            combination_strains[combo_key].add(strain_id)

    # 如果没有毒株，直接返回
    if strain_total == 0:
        return mutation_counts, mutation_strains, combination_counts, combination_strains, {}, {}

    # 筛选主要突变（频率 > 30%）
    major_mutations = set()
    for mut, count in mutation_counts.items():
        if count / strain_total >= 0.3:
            major_mutations.add(mut)
    
    # 使用主要突变重新计算组合频率
    major_combo_counts = defaultdict(int)
    major_combo_strains = defaultdict(set)
    for strain_id, mutations in strains:
        # 只保留主要突变
        major_only = [mut for mut in mutations if mut in major_mutations]
        if not major_only:
            continue
        combo_key = '-'.join(sorted(major_only))
        major_combo_counts[combo_key] += 1
        major_combo_strains[combo_key].add(strain_id)
    
    return mutation_counts, mutation_strains, combination_counts, combination_strains, major_combo_counts, major_combo_strains

def save_results(base_name, output_dir, mutation_counts, mutation_strains, combination_counts, combination_strains, major_combo_counts, major_combo_strains):
    os.makedirs(output_dir, exist_ok=True)

    # 原始突变频率输出
    mut_outfile = os.path.join(output_dir, f"{base_name}_mutation_frequencies.tsv")
    with open(mut_outfile, "w", encoding='utf-8') as f:
        f.write("Mutation\tCount\tStrains\n")
        for mut, count in sorted(mutation_counts.items(), key=lambda x: x[1], reverse=True):
            strains = ",".join(sorted(mutation_strains[mut]))
            f.write(f"{mut}\t{count}\t{strains}\n")

    # 原始组合频率输出
    combo_outfile = os.path.join(output_dir, f"{base_name}_combination_frequencies.tsv")
    with open(combo_outfile, "w", encoding='utf-8') as f:
        f.write("Combination\tCount\tStrains\n")
        for combo, count in sorted(combination_counts.items(), key=lambda x: x[1], reverse=True):
            strains = ",".join(sorted(combination_strains[combo]))
            f.write(f"{combo}\t{count}\t{strains}\n")
    
    # 主要突变组合输出（新增）
    major_outfile = os.path.join(output_dir, f"{base_name}_major_combination_frequencies.tsv")
    with open(major_outfile, "w", encoding='utf-8') as f:
        f.write("MajorCombination\tCount\tStrains\n")
        for combo, count in sorted(major_combo_counts.items(), key=lambda x: x[1], reverse=True):
            strains = ",".join(sorted(major_combo_strains[combo]))
            f.write(f"{combo}\t{count}\t{strains}\n")

def main():
    parser = argparse.ArgumentParser(description="逐个处理TSV文件，统计突变频率和组合频率")
    parser.add_argument('--input_dir', required=True, help='包含多个.tsv文件的目录')
    parser.add_argument('--output_dir', required=True, help='输出结果的目录')
    args = parser.parse_args()

    for filename in os.listdir(args.input_dir):
        if filename.endswith(".tsv"):
            filepath = os.path.join(args.input_dir, filename)
            base_name = os.path.splitext(filename)[0]
            print(f"处理文件: {filename}")
            results = process_file(filepath)
            save_results(base_name, args.output_dir, *results)

if __name__ == "__main__":
    main()

