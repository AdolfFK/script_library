#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: run_cdr_genotype.py
功能: 一键式抗体 CDR + Germline 分析工具（AbNumber 提取 CDR, ANARCI 基因型分配,结果合并）
整合了：
    1. AbNumber 提取 CDR（支持 IMGT/Kabat）
    2. ANARCI 基因型分配（自动调用）
    3. 结果合并
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-12-13  adolf1  整合了功能
依赖:  Python≥3.8 pandas≥1.5
用法:  
        # 默认（人类抗体，最常用）
        python run_cdr_genotype.py -i input.fasta -o result.tsv -s imgt

        # 指定人类
        python run_cdr_genotype.py -i input.fasta -o result.tsv --species human

        # 指定小鼠（如果你的数据是 mouse 抗体）
        python run_cdr_genotype.py -i input.fasta -o result.tsv --species mouse

        # 不指定物种（让 ANARCI 自动在 human/mouse 中选最佳匹配）
        python run_cdr_genotype.py -i input.fasta -o result.tsv --species auto
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""

# 脚本详情
import argparse
import os
import sys
import subprocess
import tempfile
from Bio import SeqIO
from abnumber import Chain


def run_anarci(input_fasta, species="human", anarci_output_dir=None):
    """
    调用 ANARCI，确保最终返回一个以 .anarci 结尾的文件路径。
    兼容 ANARCI 不同版本（有的加后缀，有的不加）。
    """
    if anarci_output_dir is None:
        anarci_output_dir = tempfile.mkdtemp(prefix="anarci_")
    
    # 使用输入文件名作为基础（可选，也可保留 Numbered_sequences）
    base_name = os.path.splitext(os.path.basename(input_fasta))[0]
    prefix = os.path.join(anarci_output_dir, base_name)
    
    cmd = [
        "ANARCI",
        "-i", input_fasta,
        "-o", prefix,
        "--scheme", "imgt",
        "--assign_germline"
    ]
    if species.lower() != "auto":
        cmd += ["--use_species", species.lower()]
    
    print(f"Running ANARCI: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("❌ ANARCI 执行失败！")
        print(result.stderr)
        sys.exit(1)
    
    # ANARCI 可能生成 prefix.anarci 或 prefix（取决于版本）
    anarci_file_with_suffix = prefix + ".anarci"
    anarci_file_no_suffix = prefix
    
    final_anarci_file = None
    
    if os.path.exists(anarci_file_with_suffix):
        final_anarci_file = anarci_file_with_suffix
    elif os.path.exists(anarci_file_no_suffix):
        # 重命名为标准 .anarci 后缀
        final_anarci_file = anarci_file_with_suffix
        shutil.move(anarci_file_no_suffix, final_anarci_file)
        print(f"ℹ️  重命名无后缀文件为: {final_anarci_file}")
    else:
        print("❌ ANARCI 未生成任何输出文件！")
        print(f"  检查路径: {anarci_file_with_suffix} 或 {anarci_file_no_suffix}")
        sys.exit(1)
    
    return final_anarci_file, anarci_output_dir


def parse_anarci_dot_anarci(anarci_file):
    """适用于单 domain 输入（每条 FASTA 序列只对应一个 germline）"""
    germline_dict = {}
    with open(anarci_file, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("# ") and not line.startswith("##") and len(line) > 2:
            seq_id = line[2:].strip()
            # 查找接下来的 germline 行
            j = i + 1
            germline_line = None
            while j < len(lines):
                l = lines[j].strip()
                if l == "# Most sequence-identical germlines":
                    # 下下一行是数据（跳过 header）
                    if j + 2 < len(lines):
                        germline_line = lines[j + 2].strip()
                        break
                j += 1

            if germline_line and germline_line.startswith("#|"):
                parts = [p.strip() for p in germline_line[2:].split("|")]
                if len(parts) >= 5:
                    germline_dict[seq_id] = {
                        "species": parts[0],
                        "v_gene": parts[1] or "N/A",
                        "v_identity": parts[2] or "N/A",
                        "j_gene": parts[3] or "N/A",
                        "j_identity": parts[4] or "N/A",
                    }
        i += 1

    return germline_dict


def extract_cdrs_with_abnumber(input_fasta, scheme='imgt'):
    """AbNumber 提取 CDR"""
    cdr_records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_id = record.id
        sequence = str(record.seq).replace("-", "")  # 去除 gap

        try:
            chain = Chain(sequence, scheme=scheme.lower())
            chain_type = getattr(chain, 'chain_type', 'H' if len(chain.cdr3_seq or '') > 15 else 'L')

            cdr_records.append({
                "id": seq_id,
                "chain_type": chain_type,
                "CDR1": chain.cdr1_seq or "",
                "CDR2": chain.cdr2_seq or "",
                "CDR3": chain.cdr3_seq or "",
            })
        except Exception as e:
            print(f"[AbNumber] Error {seq_id}: {e}")
            cdr_records.append({
                "id": seq_id,
                "chain_type": "?",
                "CDR1": "ERROR",
                "CDR2": "ERROR",
                "CDR3": "ERROR",
            })
    return cdr_records


def main():
    parser = argparse.ArgumentParser(
        description="一键抗体分析：CDR提取 + 从 .anarci 文件提取 Germline")
    parser.add_argument("-i", "--input", required=True, help="输入 FASTA 文件")
    parser.add_argument("-o", "--output", required=True, help="输出 TSV 文件")
    parser.add_argument("-s", "--scheme", choices=["imgt", "kabat"], default="imgt",
                        help="CDR 方案（默认: imgt）")
    parser.add_argument("--species", choices=["human", "mouse", "auto"], default="human",
                        help="物种（默认: human）")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"❌ 文件不存在: {args.input}")
        sys.exit(1)

    print(f"步骤 1/3: AbNumber 提取 CDR ({args.scheme.upper()})")
    cdr_data = extract_cdrs_with_abnumber(args.input, scheme=args.scheme)

    print(f"步骤 2/3: ANARCI 生成 .anarci 文件 (物种: {args.species})")
    anarci_file, _ = run_anarci(args.input, species=args.species)

    print("步骤 3/3: 从 .anarci 文件解析 germline 并合并")
    germline_data = parse_anarci_dot_anarci(anarci_file)

    with open(args.output, "w") as out:
        header = ["id", "chain_type", "CDR1", "CDR2", "CDR3",
                  "v_gene", "j_gene", "v_identity", "j_identity", "species"]
        out.write("\t".join(header) + "\n")

        for rec in cdr_data:
            seq_id = rec["id"]
            g = germline_data.get(seq_id, {
                "v_gene": "N/A", "j_gene": "N/A",
                "v_identity": "N/A", "j_identity": "N/A", "species": "N/A"
            })
            line = [
                seq_id, rec["chain_type"],
                rec["CDR1"], rec["CDR2"], rec["CDR3"],
                g["v_gene"], g["j_gene"],
                str(g["v_identity"]), str(g["j_identity"]),
                g["species"]
            ]
            out.write("\t".join(line) + "\n")

    print(f"✅ 完成！结果保存至: {args.output}")
    print(f"   总序列: {len(cdr_data)}，Germline 成功解析: {len(germline_data)} 条")


if __name__ == "__main__":
    try:
        import abnumber
    except ImportError:
        print("❌ 请安装: pip install abnumber")
        sys.exit(1)

    import shutil
    if shutil.which("ANARCI") is None:
        print("❌ 未找到 ANARCI，请安装: conda install -c conda-forge anarci")
        sys.exit(1)

    main()



