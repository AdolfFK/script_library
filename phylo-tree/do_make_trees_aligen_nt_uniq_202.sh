#!/usr/bin/env bash
#==============================================================================
# 脚本名: do_make_trees_aligen_nt_uniq_202.sh
# 功能:  建立已经比对好的核酸序列的进化树（增强版 + 严格过滤 + 可选跳过）
# 功能细节：
#       1. 数据质量控制与过滤
#       2. 去除冗余序列
#       3. 使用 MAFFT 重新比对
#       4. 使用 trimAl 严格过滤（可通过参数跳过）
#       5. 使用 FastTree 和 IQ-TREE 构建最大似然树
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-06-10  adolf1  增加 xxx 功能
# 依赖:  bash≥4.0  curl≥7.0
# 用法:  ./do_make_trees_aligen_nt_uniq_202.sh <alignment_file>
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================
set -euo pipefail

# ======================
# 参数检查
# ======================
if [ -z "$1" ]; then
    echo "Usage: $0 <fasta_file> [--skip-trimal]"
    exit 1
fi

filename="$1"
filehead="${filename%.*}"
skip_trimal=false

if [ "$2" == "--skip-trimal" ]; then
    skip_trimal=true
fi

# ======================
# 软件路径定义
# ======================
fasttree="/data/workdata/wangjy/tools/Fasttree/FastTreeMP"
mafft="/data/workdata/wangjy/software/mafft_7525/bin/mafft"
trimal="/data/workdata/wangjy/software/trimal-1.5.0/source/trimal"
conda_activate="/data/workdata/wangjy/software/anaconda3/conda3/bin/activate"
env_name="biosoft_env"

# ======================
# 去重与统计
# ======================
echo "[Step 1] Remove duplicate sequences..."
seqkit stat "${filename}" > "${filehead}.fasta.stat"
seqkit rmdup -s "${filename}" > "${filehead}.uniq.fasta"
seqkit stat "${filehead}.uniq.fasta" > "${filehead}.uniq.fasta.stat"

# ======================
# 质量过滤
# ======================
echo "[Step 2] Filter low-quality sequences..."
avg_len=$(seqkit fx2tab "${filehead}.uniq.fasta" | awk '{sum+=length($2)} END {print sum/NR}')
min_len=$(awk -v avg="$avg_len" 'BEGIN {print int(avg*0.8)}')

seqkit fx2tab "${filehead}.uniq.fasta" | \
awk -v min_len="$min_len" '{
    seq = toupper($2);
    len = length(seq);
    n_count = gsub(/N/, "", seq);
    if (len >= min_len && (n_count/len) <= 0.1) {
        print ">"$1"\n"$2;
    } else {
        print $1"\t"len"\t"n_count >> "'${filehead}'.filtered_out.txt";
    }
}' > "${filehead}.filtered.fasta"

seqkit stat "${filehead}.filtered.fasta" > "${filehead}.filtered.fasta.stat"

count=$(grep -c "^>" "${filehead}.filtered.fasta")
if [ "$count" -lt 4 ]; then
    echo "[Error] Too few valid sequences (<4). Cannot build tree."
    exit 1
fi
echo "Retained ${count} high-quality sequences."

# ======================
# MAFFT 比对
# ======================
echo "[Step 3] MAFFT alignment..."
start_time=$(date +%s)
"${mafft}" --auto --thread 16 "${filehead}.filtered.fasta" > "${filehead}.filtered.aln.fasta"
end_time=$(date +%s)
run_time_mafft=$(( end_time - start_time ))
echo "MAFFT alignment completed in ${run_time_mafft} seconds."

# ======================
# trimAl 严格过滤
# ======================
if [ "$skip_trimal" = false ]; then
    echo "[Step 4] Running trimAl (strict filtering)..."
    start_time=$(date +%s)
    "${trimal}" -in "${filehead}.filtered.aln.fasta" \
                -out "${filehead}.filtered.trimmed.fasta" \
                -gt 0.8 -cons 60 -w 3
    end_time=$(date +%s)
    run_time_trimal=$(( end_time - start_time ))
    echo "trimAl filtering completed in ${run_time_trimal} seconds."
    aln_file="${filehead}.filtered.trimmed.fasta"
else
    echo "[Step 4] Skipping trimAl filtering."
    aln_file="${filehead}.filtered.aln.fasta"
fi

# ======================
# FastTree 建树
# ======================
echo "[Step 5] Building ML tree with FastTree..."
start_time=$(date +%s)
"${fasttree}" -gtr -gamma "${aln_file}" > "${filehead}.fasttree.nwk"
end_time=$(date +%s)
run_time_fasttree=$(( end_time - start_time ))
echo "FastTree completed in ${run_time_fasttree} seconds."

# ======================
# IQ-TREE 建树
# ======================
echo "[Step 6] Building ML tree with IQ-TREE..."
start_time=$(date +%s)
source "${conda_activate}" "${env_name}"

iqtree2 -s "${aln_file}" \
        -m GTR+G4+I \
        -bb 1000 -alrt 1000 \
        -nt 16 \
        -pre "${filehead}.iqtree"

end_time=$(date +%s)
run_time_iqtree=$(( end_time - start_time ))
echo "IQ-TREE completed in ${run_time_iqtree} seconds."

# ======================
# 总结
# ======================
echo "----------------------------------------"
echo "Input file: ${filename}"
echo "Alignment file: ${aln_file}"
echo "FastTree tree: ${filehead}.fasttree.nwk"
echo "IQ-TREE tree: ${filehead}.iqtree.treefile"
echo "----------------------------------------"
echo "Runtime summary:"
echo "MAFFT: ${run_time_mafft} s"
if [ "$skip_trimal" = false ]; then
    echo "trimAl: ${run_time_trimal} s"
fi
echo "FastTree: ${run_time_fasttree} s"
echo "IQ-TREE: ${run_time_iqtree} s"
echo "----------------------------------------"
