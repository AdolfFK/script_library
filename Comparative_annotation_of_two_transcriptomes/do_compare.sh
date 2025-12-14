#!/bin/bash

#!/usr/bin/env Rscript
# 
#==============================================================================
# 脚本名: do_compare.R
# 功能:  compare_enrichment.R脚本的运行样例
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-14  adolf1
# 依赖:  R≥4.0 data.table ggplot2
# 用法:  bash do_compare.R
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================

# 使用rsv作为背景，获得相同注释通路的结果

# 2天组和4天组
# 分析GO结果
## BP
Rscript compare_enrichment.R -t GO -a tunshu_pathway_out_csv/m2dvswt_GO_BP.csv -b rsv_pathway_out_csv/RSV075-1VSRSV075-4_up_GO_BP.csv -o output_dir_2v4_BP -g BP
## CC
Rscript compare_enrichment.R -t GO -a tunshu_pathway_out_csv/m2dvswt_GO_CC.csv -b rsv_pathway_out_csv/RSV075-1VSRSV075-4_up_GO_CC.csv -o output_dir_2v4_CC -g CC
## MF
Rscript compare_enrichment.R -t GO -a tunshu_pathway_out_csv/m2dvswt_GO_MF.csv -b rsv_pathway_out_csv/RSV075-1VSRSV075-4_up_GO_MF.csv -o output_dir_2v4_MF -g CC

# 分析KEGG结果
Rscript compare_enrichment.R -t KEGG -a tunshu_pathway_out_csv/m2dvswt_KEGG.csv -b rsv_pathway_out_csv/RSV075-1VSRSV075-4_up_KEGG.csv -o output_dir_2v4_KEGG -k mmu

# 6天组和7天组
# 分析GO结果
## BP
Rscript compare_enrichment.R -t GO -a tunshu_pathway_out_csv/m6dvswt_GO_BP.csv -b rsv_pathway_out_csv/RSV075-3VSRSV075-4_up_GO_BP.csv -o output_dir_6v7_BP -g BP
## CC
Rscript compare_enrichment.R -t GO -a tunshu_pathway_out_csv/m6dvswt_GO_CC.csv -b rsv_pathway_out_csv/RSV075-3VSRSV075-4_up_GO_CC.csv -o output_dir_6v7_CC -g CC
## MF
Rscript compare_enrichment.R -t GO -a tunshu_pathway_out_csv/m6dvswt_GO_MF.csv -b rsv_pathway_out_csv/RSV075-3VSRSV075-4_up_GO_MF.csv -o output_dir_6v7_MF -g CC

# 分析KEGG结果
Rscript compare_enrichment.R -t KEGG -a tunshu_pathway_out_csv/m6dvswt_KEGG.csv -b rsv_pathway_out_csv/RSV075-3VSRSV075-4_up_KEGG.csv -o output_dir_6v7_KEGG -k mmu
