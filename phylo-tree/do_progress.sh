#!/usr/bin/env bash
#==============================================================================
# 脚本名: do_progress.sh
# 功能:  计算clade之间相似性的流程脚本
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-06-10  adolf1  增加 xxx 功能
# 依赖:  bash≥4.0  curl≥7.0
# 用法:  ./do_progress.sh
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================
set -euo pipefail

# 提取使用的序列数据
python extract_fasta_use_ids.py pt.list.id.txt pt.1.pre.fasta pt.1.pre.use.fasta

python extract_fasta_use_ids.py pt.list.id.txt pt.2.pre.fasta pt.2.pre.use.fasta

python extract_fasta_use_ids.py pt.list.id.txt pt.3.pre.fasta pt.3.pre.use.fasta

python extract_fasta_use_ids.py pt.list.id.txt pt.4.pre.fasta pt.4.pre.use.fasta

python extract_fasta_use_ids.py pt.list.id.txt pt.5.pre.fasta pt.5.pre.use.fasta


python calculate_average_similarity.py -i pt.1.pre.use.fasta -o pt.1.pre.use.matrix.txt -t dna
#DNA similarity: 0.9995

python calculate_average_similarity.py -i pt.2.pre.use.fasta -o pt.2.pre.use.matrix.txt -t dna
#DNA similarity: 0.9998

python calculate_average_similarity.py -i pt.3.pre.use.fasta -o pt.3.pre.use.matrix.txt -t dna
#DNA similarity: 0.9992

python calculate_average_similarity.py -i pt.4.pre.use.fasta -o pt.4.pre.use.matrix.txt -t dna
#DNA similarity: 0.9998

python calculate_average_similarity.py -i pt.5.pre.use.fasta -o pt.5.pre.use.matrix.txt -t dna
#DNA similarity: 0.9998

python 4_calculate_mutations_no_meta_transcript.py --input pt.1.pre.use.fasta --mafft-path /opt/homebrew/bin/mafft --output-prefix output_p1/pt.1.mutation

python 4_calculate_mutations_no_meta_transcript.py --input pt.2.pre.use.fasta --mafft-path /opt/homebrew/bin/mafft --output-prefix output_p2/pt.2.mutation

python 4_calculate_mutations_no_meta_transcript.py --input pt.3.pre.use.fasta --mafft-path /opt/homebrew/bin/mafft --output-prefix output_p3/pt.3.mutation

python 4_calculate_mutations_no_meta_transcript.py --input pt.4.pre.use.fasta --mafft-path /opt/homebrew/bin/mafft --output-prefix output_p4/pt.4.mutation

python 4_calculate_mutations_no_meta_transcript.py --input pt.5.pre.use.fasta --mafft-path /opt/homebrew/bin/mafft --output-prefix output_p5/pt.5.mutation


#生成画热图的数据
#
python generate_similarity_matrix.py -i xsd.txt -o xsd.matirx.txt

