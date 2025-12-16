#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: calculate_average_similarity_v2.py
功能:  计算比对好的msa文件中序列的相似性，输出总的相似性和平均相似性
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-06-10  adolf1  增加head说明
依赖:  Python≥3.8 pandas≥1.5
用法:  python calculate_average_similarity_v2.py -h
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""

import argparse
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import numpy as np
import os
import logging
from datetime import datetime

def setup_logger(output_dir, output_prefix):
    """Set up logging to both console and file"""
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    # Remove existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Create console handler
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # Create file handler
    log_file = os.path.join(output_dir, f"{output_prefix}.log")
    fh = logging.FileHandler(log_file)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    return logger

def read_and_preprocess_fasta(input_file, min_length=10):
    """
    Read FASTA file and preprocess sequences:
    - Merge multi-line sequences into single line
    - Filter sequences shorter than min_length
    Returns list of (header, sequence) tuples
    """
    sequences = []
    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq).upper().replace("-", "")  # Remove gaps if present
        if len(seq) >= min_length:
            sequences.append((record.id, seq))
        else:
            logging.warning(f"Sequence {record.id} is too short ({len(seq)} < {min_length}), skipped")
    return sequences

def detect_sequence_type(sequences):
    """
    Auto-detect sequence type (DNA or protein)
    """
    protein_chars = set("ACDEFGHIKLMNPQRSTVWY*")
    dna_chars = set("ACGTN")
    
    for _, seq in sequences:
        unique_chars = set(seq)
        if unique_chars - dna_chars:
            if unique_chars - protein_chars:
                raise ValueError(f"Unrecognized sequence characters: {unique_chars - protein_chars}")
            return 'protein'
    return 'dna'

def calculate_similarity_matrix(sequences, seq_type):
    """
    Calculate pairwise sequence similarity matrix
    """
    n = len(sequences)
    matrix = np.zeros((n, n))
    
    # Set up aligner based on sequence type
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # global alignment
    
    if seq_type == 'dna':
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        logging.info("Using DNA alignment parameters")
    else:  # protein
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        logging.info("Using protein alignment with BLOSUM62 matrix")
    
    # Calculate pairwise similarities
    logging.info(f"Starting similarity calculation for {n} sequences...")
    start_time = datetime.now()
    
    for i in range(n):
        _, seq1 = sequences[i]
        matrix[i][i] = 1.0  # diagonal is 100% similarity
        
        for j in range(i+1, n):
            _, seq2 = sequences[j]
            try:
                # Calculate alignment scores
                score = aligner.score(seq1, seq2)
                score1 = aligner.score(seq1, seq1)
                score2 = aligner.score(seq2, seq2)
                
                # Normalized similarity score
                norm_score = 2 * score / (score1 + score2)
                matrix[i][j] = norm_score
                matrix[j][i] = norm_score
            except Exception as e:
                logging.warning(f"Error aligning {sequences[i][0]} and {sequences[j][0]}: {str(e)}")
                matrix[i][j] = 0
                matrix[j][i] = 0
        
        # Log progress every 10 sequences
        if (i+1) % 10 == 0 or (i+1) == n:
            elapsed = (datetime.now() - start_time).total_seconds()
            logging.info(f"Processed {i+1}/{n} sequences (elapsed: {elapsed:.2f}s)")
    
    total_time = (datetime.now() - start_time).total_seconds()
    logging.info(f"Similarity matrix calculation completed in {total_time:.2f} seconds")
    
    # Calculate and log average similarity
    upper_triangle = matrix[np.triu_indices(n, k=1)]
    if len(upper_triangle) > 0:
        avg_similarity = np.mean(upper_triangle)
        logging.info(f"Average pairwise similarity: {avg_similarity:.4f}")
    else:
        avg_similarity = 0.0
        logging.info("No pairwise comparisons made (only one sequence)")
    
    return matrix, avg_similarity

def save_results(output_dir, matrix, sequences, output_prefix):
    """
    Save similarity matrix and sequence IDs
    """
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info(f"Created output directory: {output_dir}")
    
    # Save matrix as numpy binary
    np.save(os.path.join(output_dir, f"{output_prefix}.npy"), matrix)
    logging.info(f"Saved binary matrix to {os.path.join(output_dir, f'{output_prefix}.npy')}")
    
    # Save matrix as text file with headers
    with open(os.path.join(output_dir, f"{output_prefix}.txt"), 'w') as f:
        # Write header row
        headers = [seq_id for seq_id, _ in sequences]
        f.write("\t" + "\t".join(headers) + "\n")
        
        # Write matrix rows
        for i, (seq_id, _) in enumerate(sequences):
            row = [f"{val:.4f}" for val in matrix[i]]
            f.write(seq_id + "\t" + "\t".join(row) + "\n")
    logging.info(f"Saved text matrix to {os.path.join(output_dir, f'{output_prefix}.txt')}")
    
    # Save sequence IDs
    with open(os.path.join(output_dir, f"{output_prefix}_ids.txt"), 'w') as f:
        for seq_id, _ in sequences:
            f.write(f"{seq_id}\n")
    logging.info(f"Saved sequence IDs to {os.path.join(output_dir, f'{output_prefix}_ids.txt')}")

class SequenceSimilarityCalculator:
    """Sequence Similarity Calculator class"""
    
    def __init__(self, input_file, output_dir, output_prefix='similarity_matrix', 
                 seq_type='auto', min_length=10):
        self.input_file = input_file
        self.output_dir = output_dir
        self.output_prefix = output_prefix
        self.seq_type = seq_type
        self.min_length = min_length
        self.logger = self.setup_logger()
        
    def setup_logger(self):
        """Set up logging to both console and file"""
        return setup_logger(self.output_dir, self.output_prefix)
    
    def run_analysis(self):
        """Run the entire analysis workflow"""
        try:
            # Read and preprocess sequences
            self.logger.info(f"Reading sequences from {self.input_file}...")
            sequences = read_and_preprocess_fasta(self.input_file, self.min_length)
            if not sequences:
                raise ValueError(f"No sequences found longer than {self.min_length} in {self.input_file}")
            self.logger.info(f"Found {len(sequences)} sequences after filtering.")
            
            # Determine sequence type
            if self.seq_type == 'auto':
                self.seq_type = detect_sequence_type(sequences)
                self.logger.info(f"Auto-detected sequence type: {self.seq_type}")
            else:
                self.logger.info(f"Using specified sequence type: {self.seq_type}")
            
            # Calculate similarity matrix
            self.logger.info(f"Calculating similarity matrix for {len(sequences)} {self.seq_type} sequences...")
            similarity_matrix, avg_similarity = calculate_similarity_matrix(sequences, self.seq_type)
            
            # Save results
            save_results(self.output_dir, similarity_matrix, sequences, self.output_prefix)
            self.logger.info(f"Analysis completed successfully")
            
            return similarity_matrix, avg_similarity
        
        except Exception as e:
            self.logger.error(f"Error encountered: {str(e)}", exc_info=True)
            raise

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Calculate sequence similarity matrix from FASTA file.')
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output-dir', required=True, help='Output directory')
    parser.add_argument('-p', '--prefix', default='similarity_matrix',
                       help='Output file prefix (default: similarity_matrix)')
    parser.add_argument('-t', '--type', choices=['auto', 'dna', 'protein'], default='auto',
                       help='Sequence type (auto detect by default)')
    parser.add_argument('--min-length', type=int, default=10,
                       help='Minimum sequence length to consider (default: 10)')
    
    args = parser.parse_args()
    
    calculator = SequenceSimilarityCalculator(
        args.input,
        args.output_dir,
        args.prefix,
        args.type,
        args.min_length
    )
    
    try:
        similarity_matrix, avg_similarity = calculator.run_analysis()
        return similarity_matrix, avg_similarity
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        return None, None

if __name__ == "__main__":
    main()
    