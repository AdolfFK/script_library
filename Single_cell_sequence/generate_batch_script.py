#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: generate_batch_script.py
功能:  这个脚本用来批量生成单细胞初步的脚本
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-06-10  adolf1  增加 xxx 功能
依赖:  Python≥3.8 pandas≥1.5
用法:  python generate_batch_script.py -参数
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""


import argparse
import pandas as pd
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate batch submission bash script from sample file.')
    parser.add_argument('--sample_file', required=True, help='Path to the sample input file (CSV format with id and group_name columns).')
    parser.add_argument('--output_dir', required=True, help='Directory to save the generated bash script.')
    parser.add_argument('--data_dir', default='/data/workdata/wangjy/xiey/singlecell/data', help='Base directory for input data.')
    parser.add_argument('--out_base_dir', default='/data/workdata/wangjy/xiey/singlecell/output/1_seruat', help='Base directory for output.')
    parser.add_argument('--rscript_path', default='/data/workdata/wangjy/xiey/singlecell/script/run_MultiSamples.V5_gpt.R', help='Path to the R script.')
    parser.add_argument('--conda_env_path', default='/data/workdata/wangjy/software/anaconda3/conda3/bin/activate', help='Path to conda activate script.')
    parser.add_argument('--conda_env_name', default='R_conda_env', help='Name of the conda environment.')
    parser.add_argument('--mtfilter', default='10', help='Mitochondrial filter value.')
    parser.add_argument('--umifiltermax', default='40000', help='Maximum UMI filter value.')
    parser.add_argument('--ngenefiltermax', default='4000', help='Maximum gene filter value.')
    parser.add_argument('--resolution', default='1.0', help='Resolution parameter for clustering.')
    return parser.parse_args()

def read_sample_file(sample_file):
    return pd.read_csv(sample_file, sep='\t')

def generate_bash_script(df, args):
    # Group by group_name
    grouped = df.groupby('group_name')['id'].apply(list).to_dict()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Generate a bash script for each group
    for group_name, ids in grouped.items():
        # Construct the expM parameter
        expM = ','.join([f"{args.data_dir}/{id}/count/sample_raw_feature_bc_matrix" for id in ids])
        
        # Construct datatype (all 2s based on example)
        datatype = ','.join(['2'] * len(ids))
        
        # Construct ugname (use group_name directly)
        ugname = ','.join([group_name] * len(ids))
        
        # Construct spname (list of IDs)
        spname = ','.join(map(str, ids))
        
        # Construct output directory for this group, including group_name
        outdir = os.path.join(args.out_base_dir, group_name)
        
        # Create the bash script content
        script_content = f"""#!/bin/bash

source {args.conda_env_path} {args.conda_env_name}

Rscript {args.rscript_path} \\
 --expM {expM} \\
 --datatype {datatype} \\
 --ugname {ugname} \\
 --spname {spname} \\
 --mtfilter {args.mtfilter} \\
 --umifiltermax {args.umifiltermax} \\
 --ngenefiltermax {args.ngenefiltermax} \\
 --prefix {group_name} \\
 --resolution {args.resolution} \\
 --outdir {outdir}
"""
        # Write the script to a file
        output_file = os.path.join(args.output_dir, f"run_{group_name}.sh")
        with open(output_file, 'w') as f:
            f.write(script_content)
        
        # Make the script executable
        os.chmod(output_file, 0o755)

def main():
    args = parse_arguments()
    df = read_sample_file(args.sample_file)
    generate_bash_script(df, args)
    print(f"Bash scripts generated in {args.output_dir}")

if __name__ == "__main__":
    main()

