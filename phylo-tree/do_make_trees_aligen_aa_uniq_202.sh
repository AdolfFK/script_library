#!/usr/bin/env bash
#==============================================================================
# 脚本名: do_make_trees_aligen_aa_uniq_202.sh
# 功能:  建立已经比对好的氨基酸序列的进化树（增强版 + 严格过滤 + 可选跳过）
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
# 用法:  ./do_make_trees_aligen_aa_uniq_202.sh <alignment_file>
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================
set -euo pipefail

# ======================
# 输入参数检查
# ======================
if [ -z "$1" ]; then
    echo "请输入氨基酸序列文件（FASTA格式）!"
    echo "Usage: $0 <alignment_file> [--no-filter]"
    exit 1
fi

filename="$1"
filehead="${filename%.*}"

# 可选参数：跳过过滤
skip_filter=false
if [ "$2" == "--no-filter" ]; then
    skip_filter=true
fi

# ======================
# 软件路径定义
# ======================
fasttree="/data/workdata/wangjy/tools/Fasttree/FastTreeMP"
mafft="/data/workdata/wangjy/software/mafft_7525/bin/mafft"
conda_env="/data/workdata/wangjy/software/anaconda3/conda3/bin/activate"
iqtree_env="biosoft_env"

# ======================
# 数据统计与预处理
# ======================
echo "统计输入序列..."
seqkit stat "${filename}" > "${filehead}.fasta.stat"

echo "去除完全重复序列..."
seqkit rmdup -s "${filename}" > "${filehead}.uniq.fasta"
seqkit stat "${filehead}.uniq.fasta" > "${filehead}.uniq.fasta.stat"

# ======================
# 数据质量过滤
# ======================
echo "开始过滤低质量序列..."
seqkit fx2tab "${filehead}.uniq.fasta" \
| awk -F'\t' '{
    id=$1; seq=toupper($2);
    len=length(seq);
    bad=gsub(/[^ACDEFGHIKLMNPQRSTVWY]/,"",seq);
    if (len>=50 && bad/len<=0.1) {
        print ">"id"\n"toupper($2);
    } else {
        print id"\t"len"\t"bad >> "'${filehead}'.filtered_out.txt";
    }
}' > "${filehead}.filtered.fasta"

seqkit stat "${filehead}.filtered.fasta" > "${filehead}.filtered.stat"

count=$(grep -c "^>" "${filehead}.filtered.fasta")
if [ "$count" -lt 4 ]; then
    echo "有效序列数量太少 (<4)，无法建树！"
    exit 1
fi
echo "过滤后保留 ${count} 条高质量序列"
echo "被剔除的序列列表已保存至: ${filehead}.filtered_out.txt"

# ======================
# MAFFT 重新比对
# ======================
echo "使用 MAFFT 重新进行氨基酸序列比对..."
start_time=$(date +%s)

"${mafft}" --auto --thread -16 "${filehead}.filtered.fasta" > "${filehead}.filtered.aligned.fasta"

end_time=$(date +%s)
run_time_mafft=$(( end_time - start_time ))
echo "MAFFT 比对完成（耗时 ${run_time_mafft} 秒）"
echo "输出文件: ${filehead}.filtered.aligned.fasta"

# ======================
# trimAl 严格过滤
# ======================
if [ "$skip_filter" = false ]; then
    echo "使用 trimAl 进行严格比对过滤..."
    source "${conda_env}" "${iqtree_env}"

    # 更严格的过滤参数：
    # -gt 0.9   保留 gap 比例 ≤ 10% 的列
    # -resoverlap 0.8   要求每列至少 80% 序列有数据
    # -seqoverlap 80    要求序列与其他序列的重叠区域 ≥80%
    trimal -in "${filehead}.filtered.aligned.fasta" \
           -out "${filehead}.filtered.trimmed.fasta" \
           -gt 0.9 -resoverlap 0.8 -seqoverlap 80

    seqkit stat "${filehead}.filtered.trimmed.fasta" > "${filehead}.filtered.trimmed.stat"
    echo "trimAl 严格过滤完成，输出文件: ${filehead}.filtered.trimmed.fasta"
    aln_for_tree="${filehead}.filtered.trimmed.fasta"
else
    echo "跳过 trimAl 过滤，直接使用 MAFFT 比对结果建树..."
    aln_for_tree="${filehead}.filtered.aligned.fasta"
fi

# ======================
# FastTree 建树
# ======================
echo "使用 FastTree 构建最大似然树..."
start_time=$(date +%s)

"${fasttree}" -wag -gamma "${aln_for_tree}" > "${filehead}.fasttree.nwk"

end_time=$(date +%s)
run_time_fasttree=$(( end_time - start_time ))
echo "FastTree 建树完毕（耗时 ${run_time_fasttree} 秒）"
echo "输出文件: ${filehead}.fasttree.nwk"

# ======================
# IQ-TREE 建树
# ======================
echo "使用 IQ-TREE 构建最大似然树（固定模型 LG+F+G4）..."
start_time=$(date +%s)

source "${conda_env}" "${iqtree_env}"

iqtree2 -s "${aln_for_tree}" -m LG+F+G4 -bb 1000 -alrt 1000 -nt AUTO -pre "${filehead}.iqtree"

end_time=$(date +%s)
run_time_iqtree=$(( end_time - start_time ))
echo "IQ-TREE 建树完毕（耗时 ${run_time_iqtree} 秒）"
echo "输出文件: ${filehead}.iqtree.treefile"

# ======================
# 结束提示
# ======================
echo "全部分析完成!"
echo "FastTree 树: ${filehead}.fasttree.nwk"
echo "IQ-TREE 树: ${filehead}.iqtree.treefile"
echo "MAFFT 比对结果: ${filehead}.filtered.aligned.fasta"
if [ "$skip_filter" = false ]; then
    echo "严格过滤后比对结果: ${filehead}.filtered.trimmed.fasta"
else
    echo "未执行 trimAl 过滤"
fi
echo "初步过滤序列: ${filehead}.filtered.fasta"
echo "被剔除序列列表: ${filehead}.filtered_out.txt"

