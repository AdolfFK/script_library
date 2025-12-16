#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: 4_calculate_mutations_no_meta_transcript.py
功能:  比对和分析DNA序列突变
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-12-16  adolf1  增加head功能
依赖:  Python≥3.8 pandas≥1.5
用法:  python 4_calculate_mutations_no_meta_transcript.py -h
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import os
import shutil
from collections import defaultdict
import math
import subprocess
import tempfile


class MutationAnalyzer:
    """A class to analyze mutations in sequences using MAFFT alignment, with DNA translation support."""

    def __init__(self, log_file):
        """Initialize the analyzer with a log file."""
        self.log_file = log_file
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        with open(self.log_file, "w") as f:
            f.write("突变分析日志\n")
            f.write("=" * 20 + "\n")

    def log(self, message):
        """Write a message to the log file."""
        with open(self.log_file, "a") as f:
            f.write(f"{message}\n")

    def detect_sequence_type(self, sequences):
        """Detect if sequences are DNA or protein."""
        test_seq = str(sequences[0].seq).upper()
        dna_chars = set('ATCGN-')
        return 'dna' if all(c in dna_chars for c in test_seq) else 'protein'

    def translate_dna_sequences(self, dna_sequences):
        """Translate DNA sequences to protein sequences, removing those with premature stop codons.
        If the length of DNA sequence is not a multiple of 3, only translate the first multiple of 3 nucleotides.
        """
        protein_sequences = []
        valid_indices = []

        for i, record in enumerate(dna_sequences):
            dna_seq = str(record.seq)
            # 确保只取 3 的倍数长度的 DNA 序列进行翻译
            dna_seq = dna_seq[:len(dna_seq) // 3 * 3]
            try:
                protein_seq = str(Seq(dna_seq).translate(to_stop=True))
                protein_sequences.append(SeqRecord(
                    Seq(protein_seq),
                    id=record.id,
                    description=""
                ))
                valid_indices.append(i)
            except Exception as e:
                self.log(f"序列 {record.id} 翻译失败，将被跳过: {str(e)}")

        valid_dna_sequences = [dna_sequences[i] for i in valid_indices]
        return valid_dna_sequences, protein_sequences

    def prepend_reference(self, input_fasta, ref_fasta):
        """Prepend reference sequence to input FASTA if not already present."""
        input_sequences = list(SeqIO.parse(input_fasta, "fasta"))
        ref_sequence = list(SeqIO.parse(ref_fasta, "fasta"))

        if len(ref_sequence) != 1:
            raise ValueError("Reference FASTA file must contain exactly one sequence!")

        ref_seq = ref_sequence[0]
        ref_id = ref_seq.id

        if ref_id in {seq.id for seq in input_sequences}:
            self.log(f"参考序列 ID {ref_id} 已存在于输入文件中，无需添加")
            return input_fasta

        self.log(f"参考序列 ID {ref_id} 未在输入文件中找到，将添加到文件第一行")

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_fasta:
            SeqIO.write(ref_seq, temp_fasta, "fasta")
            for seq in input_sequences:
                SeqIO.write(seq, temp_fasta, "fasta")

        return temp_fasta.name

    def run_mafft(self, input_file, output_file, mafft_path):
        """Run MAFFT to align sequences."""
        try:
            # 添加 --leavegappyregion 选项
            cmd = [mafft_path, "--auto", "--leavegappyregion", input_file]
            with open(output_file, "w") as f:
                subprocess.run(cmd, stdout=f, check=True)
            self.log(f"MAFFT 对齐完成，输出文件: {output_file}")
        except FileNotFoundError:
            raise RuntimeError(f"MAFFT executable not found at {mafft_path}.")
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MAFFT alignment failed: {e}")

    def calculate_shannon_entropy(self, column, exclude_gaps=True, valid_chars=None):
        """Calculate Shannon entropy for a column of residues."""
        counts = defaultdict(int)
        total = 0

        valid_chars = valid_chars or (set("ACDEFGHIKLMNPQRSTVWY") if valid_chars is None else valid_chars)

        for aa in column:
            if exclude_gaps and aa == "-":
                continue
            if aa in valid_chars or (not exclude_gaps and aa == "-"):
                counts[aa] += 1
                total += 1

        if total == 0:
            return 0.0

        entropy = 0
        for count in counts.values():
            p = count / total
            entropy -= p * math.log2(p) if p > 0 else 0

        return entropy

    def analyze_mutations(self, msa_file, ref_id, seq_type='protein'):
        """Analyze mutations in aligned sequences."""
        sequences = list(SeqIO.parse(msa_file, "fasta"))
        if not sequences:
            raise ValueError("No sequences found in the aligned file!")

        ref_seq = next((str(seq.seq).upper() for seq in sequences if seq.id == ref_id), None)
        if ref_seq is None:
            raise ValueError(f"Reference sequence ID {ref_id} not found!")

        seq_length = len(ref_seq)
        total_sequences = len(sequences)

        mutations = []
        site_stats = defaultdict(lambda: defaultdict(int))
        columns = [[] for _ in range(seq_length)]
        mutation_matrix = []

        valid_chars = set("ACDEFGHIKLMNPQRSTVWY") if seq_type == 'protein' else set("ATCG")

        for seq_record in sequences:
            seq_id = seq_record.id
            seq = str(seq_record.seq).upper()

            if len(seq) != seq_length:
                raise ValueError(f"Sequence {seq_id} length mismatch!")

            matrix_row = [seq_id]
            for i, res in enumerate(seq):
                columns[i].append(res)
                if seq_id == ref_id:
                    matrix_row.append(0)
                else:
                    ref_res = ref_seq[i]
                    if res == ref_res or res == "-":
                        matrix_row.append(0)
                    else:
                        matrix_row.append(1)

            mutation_matrix.append(matrix_row)

            if seq_id == ref_id:
                continue

            for i in range(seq_length):
                ref_res = ref_seq[i]
                seq_res = seq[i]
                if seq_res != ref_res and seq_res != "-":
                    mutations.append((seq_id, i + 1, ref_res, seq_res))
                    site_stats[i + 1][seq_res] += 1

        site_entropy = [self.calculate_shannon_entropy(col, valid_chars=valid_chars) for col in columns]

        results = []
        for pos in range(1, seq_length + 1):
            ref_res = ref_seq[pos - 1]
            if ref_res == "-":
                continue

            mutations_at_site = site_stats.get(pos, {})
            total_mutations = sum(mutations_at_site.values())
            mutation_rate = total_mutations / total_sequences
            max_mutation = max(mutations_at_site.items(), key=lambda x: x[1], default=(None, 0))
            max_res, max_count = max_mutation if max_mutation[0] else (None, 0)

            detail = max_res if max_res else ""
            all_detail = ":".join([f"{res}={count}" for res, count in sorted(mutations_at_site.items())])
            all_detail = f"{ref_res}={total_sequences - total_mutations}:{all_detail}:" if all_detail else f"{ref_res}={total_sequences}:"

            results.append({
                "Number": pos,
                "refRes": ref_res,
                "TotalSites": total_sequences,
                "mutations": mutation_rate,
                "max": max_count,
                "detail": detail,
                "all_detail": all_detail,
                "entropy": site_entropy[pos - 1]
            })

        self.log(f"{seq_type.upper()} 突变分析完成: 参考序列={ref_id}, 长度={seq_length}, 总序列数={total_sequences}, 总突变数={len(mutations)}")

        return ref_id, ref_seq, mutations, results, site_entropy, mutation_matrix

    def write_output_files(self, output_prefix, ref_id, mutations, results, mutation_matrix, seq_length, seq_type):
        """Write all output files for a given sequence type."""
        output_dir = os.path.dirname(output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        prefix = f"{output_prefix}_{seq_type}"

        mut_file = f"{prefix}_mutations.tsv"
        self.write_mutations(mutations, mut_file)

        stats_file = f"{prefix}_site_stats.tsv"
        self.write_site_stats(results, stats_file)

        matrix_file = f"{prefix}_matrix.tsv"
        self.write_mutation_matrix(mutation_matrix, seq_length, matrix_file)

        return mut_file, stats_file, matrix_file

    def write_mutations(self, mutations, output_file):
        """Write per-sequence mutations to a TSV file."""
        mutation_dict = defaultdict(list)
        for seq_id, pos, ref_res, mut_res in mutations:
            mutation_dict[seq_id].append(f"{ref_res}{pos}{mut_res}")

        with open(output_file, "w") as f:
            f.write("id\tmutations\n")
            for seq_id, mutation_list in sorted(mutation_dict.items()):
                f.write(f"{seq_id}\t{','.join(mutation_list)}\n")
        self.log(f"突变记录写入: {output_file}")

    def write_site_stats(self, results, output_file):
        """Write site-wise mutation statistics to a TSV file."""
        with open(output_file, "w") as f:
            f.write("Number\trefRes\tTotalSites\tmutations\tmax\tdetail\tall_detail\tentropy\n")
            for result in results:
                f.write(f"{result['Number']}\t{result['refRes']}\t{result['TotalSites']}\t"
                        f"{result['mutations']}\t{result['max']}\t{result['detail']}\t"
                        f"{result['all_detail']}\t{result['entropy']:.4f}\n")
        self.log(f"位点统计写入: {output_file}")

    def write_mutation_matrix(self, mutation_matrix, seq_length, output_file):
        """Write mutation matrix to a TSV file."""
        with open(output_file, "w") as f:
            header = ["id"] + [str(i) for i in range(1, seq_length + 1)]
            f.write("\t".join(header) + "\n")
            for row in mutation_matrix:
                f.write("\t".join(str(x) for x in row) + "\n")
        self.log(f"突变矩阵写入: {output_file}")

    def run(self, input_fasta, ref_fasta=None, mafft_path="mafft",
            output_prefix="output/mutation_results",
            aligned_output="aligned_sequences.fasta"):
        """Run the full mutation analysis pipeline."""
        if not os.path.exists(input_fasta):
            raise FileNotFoundError(f"Input file {input_fasta} not found!")

        output_dir = os.path.dirname(output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        sequences = list(SeqIO.parse(input_fasta, "fasta"))
        if not sequences:
            raise ValueError("Input FASTA file is empty!")

        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_fasta:
            for seq in sequences:
                seq.description = ""
                SeqIO.write(seq, temp_fasta, "fasta")

        try:
            seq_type = self.detect_sequence_type(sequences)
            self.log(f"检测到序列类型: {seq_type.upper()}")

            if ref_fasta:
                if not os.path.exists(ref_fasta):
                    raise FileNotFoundError(f"Reference file {ref_fasta} not found!")
                input_fasta = self.prepend_reference(temp_fasta.name, ref_fasta)
                ref_seq = list(SeqIO.parse(ref_fasta, "fasta"))[0]
                ref_id = ref_seq.id
            else:
                ref_id = sequences[0].id
                input_fasta = temp_fasta.name
                self.log(f"使用第一个序列作为参考: {ref_id}")

            os.makedirs(os.path.dirname(aligned_output), exist_ok=True)

            if seq_type == 'protein':
                # 蛋白序列：直接进行比对和突变分析
                self.log("开始蛋白序列比对")
                self.run_mafft(input_fasta, aligned_output, mafft_path)

                self.log("分析蛋白突变")
                aa_results = self.analyze_mutations(aligned_output, ref_id, 'protein')
                # 修改此处的参数传递
                aa_files = self.write_output_files(output_prefix, aa_results[0], aa_results[2], aa_results[3], 
                                                   aa_results[5], len(aa_results[1]), 'aa')
                nuc_files = None
            else:  # seq_type == 'dna'
                # 核酸序列：先翻译，然后核酸和蛋白都进行比对和突变分析
                self.log("翻译核酸序列为蛋白序列")
                dna_sequences = list(SeqIO.parse(input_fasta, "fasta"))
                valid_dna, protein_sequences = self.translate_dna_sequences(dna_sequences)

                if len(protein_sequences) < 2:
                    self.log("警告: 有效蛋白质序列不足，跳过蛋白质分析")
                    aa_files = None
                else:
                    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as temp_aa:
                        SeqIO.write(protein_sequences, temp_aa, "fasta")

                    aa_aligned = f"{output_prefix}_aligned_aa.fasta"
                    self.log("开始蛋白序列比对")
                    self.run_mafft(temp_aa.name, aa_aligned, mafft_path)

                    self.log("分析蛋白突变")
                    aa_results = self.analyze_mutations(aa_aligned, ref_id, 'protein')
                    # 修改此处的参数传递
                    aa_files = self.write_output_files(output_prefix, aa_results[0], aa_results[2], aa_results[3], 
                                                       aa_results[5], len(aa_results[1]), 'aa')
                    os.remove(temp_aa.name)

                # 核酸序列比对和分析
                self.log("开始核酸序列比对")
                self.run_mafft(input_fasta, aligned_output, mafft_path)

                self.log("分析核酸突变")
                nuc_results = self.analyze_mutations(aligned_output, ref_id, 'dna')
                # 修改此处的参数传递
                nuc_files = self.write_output_files(output_prefix, nuc_results[0], nuc_results[2], nuc_results[3], 
                                                    nuc_results[5], len(nuc_results[1]), 'nuc')

            self.log(f"分析完成。结果文件前缀: {output_prefix}")
            return nuc_files, aa_files

        finally:
            os.remove(temp_fasta.name)
            if ref_fasta and input_fasta != temp_fasta.name:
                os.remove(input_fasta)



def main():
    """Command-line interface for mutation analysis."""
    parser = argparse.ArgumentParser(
        description="Align sequences and analyze mutations with DNA translation support.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="输入 FASTA 文件")
    parser.add_argument("--ref", help="参考 FASTA 文件（单一序列）")
    parser.add_argument("--mafft-path", default="mafft", help="MAFFT 可执行文件路径")
    parser.add_argument("--output-prefix", default="output/mutation_results",
                        help="输出目录和文件名前缀")
    parser.add_argument("--aligned-output", help="对齐序列输出文件")
    parser.add_argument("--log-file", help="日志文件路径")

    args = parser.parse_args()

    output_dir = os.path.dirname(args.output_prefix) or "output"
    aligned_output = args.aligned_output or os.path.join(output_dir, "aligned_sequences.fasta")
    log_file = args.log_file or os.path.join(output_dir, "mutation_analysis.log")

    analyzer = MutationAnalyzer(log_file=log_file)
    analyzer.run(
        input_fasta=args.input,
        ref_fasta=args.ref,
        mafft_path=args.mafft_path,
        output_prefix=args.output_prefix,
        aligned_output=aligned_output
    )


if __name__ == "__main__":
    main()
    
