#!/usr/bin/env Rscript
# 
#==============================================================================
# 脚本名: compare_enrichment.R
# 功能:  脚本比较两次转录组分析的结果文件，输出其中的差异
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-14  adolf1
# 依赖:  R≥4.0 data.table ggplot2
# 用法:  Rscript compare_enrichment.R -a [分析1结果] -b [分析2结果] -o [输出目录] -t [GO/KEGG]
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(ggplot2)
  library(ggvenn)
  library(GOSemSim)
  library(pheatmap)
  library(org.Mm.eg.db)
  library(readxl)
  library(tools)
  library(pathview)
  library(KEGGREST)
  library(stringr)
  library(dplyr)
  library(tidyr)
})

# 设置中文编码
Sys.setlocale("LC_ALL", "UTF-8")

# 解析命令行参数
option_list <- list(
  make_option(c("-a", "--analysis1"), type="character", default=NULL,
              help="第一次分析结果文件路径(支持CSV/XLS/XLSX)", metavar="FILE"),
  make_option(c("-b", "--analysis2"), type="character", default=NULL,
              help="第二次分析结果文件路径(支持CSV/XLS/XLSX)", metavar="FILE"),
  make_option(c("-o", "--output"), type="character", default="./enrichment_comparison",
              help="输出目录路径 [默认: %default]", metavar="DIR"),
  make_option(c("-s", "--species"), type="character", default="org.Mm.eg.db",
              help="物种注释数据库 [默认: %default]", metavar="DB"),
  make_option(c("-p", "--pvalue"), type="double", default=0.05,
              help="显著性阈值 [默认: %default]", metavar="FLOAT"),
  make_option(c("-t", "--type"), type="character", default="GO",
              help="富集类型(GO或KEGG) [默认: %default]", metavar="TYPE"),
  make_option(c("-f", "--format"), type="character", default="auto",
              help="文件格式(可选: auto, csv, excel) [默认: auto检测]", metavar="FORMAT"),
  make_option(c("-k", "--kegg_organism"), type="character", default="mmu",
              help="KEGG物种代码(如mmu, hsa) [默认: %default]", metavar="KEGG_ORG"),
  make_option(c("-g", "--go_ontology"), type="character", default="BP",
              help="GO类别(BP, MF或CC) [默认: %default]", metavar="ONTOLOGY")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$analysis1) || is.null(opt$analysis2)) {
  print_help(opt_parser)
  stop("必须提供两个分析结果文件路径", call.=FALSE)
}

# 创建输出目录
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# 加载物种数据库
load_orgdb <- function(db_name) {
  if (!require(db_name, character.only = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(db_name)
    library(db_name, character.only = TRUE)
  }
  return(get(db_name))
}

tryCatch({
  org_db <- load_orgdb(opt$species)
}, error = function(e) {
  stop(paste("无法加载物种数据库:", opt$species, "\n请确保数据库名称正确且已安装"), call.=FALSE)
})

# 日志函数
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(paste("[", timestamp, "]", msg))
  write(paste(timestamp, msg), file = file.path(opt$output, "analysis.log"), append = TRUE)
}

log_message(paste("开始富集结果比较分析"))
log_message(paste("富集类型:", opt$type))
log_message(paste("分析1:", opt$analysis1))
log_message(paste("分析2:", opt$analysis2))
log_message(paste("输出目录:", opt$output))
log_message(paste("物种数据库:", opt$species))
log_message(paste("显著性阈值:", opt$pvalue))

if (opt$type == "GO") {
  log_message(paste("GO类别:", opt$go_ontology))
} else {
  log_message(paste("KEGG物种代码:", opt$kegg_organism))
}

# ---------------------------
# 1. 通用数据读取函数
# ---------------------------
read_enrich_results <- function(file_path, analysis_name) {
  if (!file.exists(file_path)) {
    stop(paste("文件不存在:", file_path), call.=FALSE)
  }
  
  # 自动检测文件类型
  if (opt$format == "auto") {
    file_ext <- tolower(file_ext(file_path))
    if (file_ext %in% c("xls", "xlsx")) {
      file_type <- "excel"
    } else if (file_ext == "csv") {
      file_type <- "csv"
    } else {
      stop(paste("不支持的文件类型:", file_ext), call.=FALSE)
    }
  } else {
    file_type <- tolower(opt$format)
  }
  
  # 读取文件
  if (file_type %in% c("excel", "xls", "xlsx")) {
    df <- read_excel(file_path)
    log_message(paste("读取Excel文件:", file_path, "| 行数:", nrow(df)))
  } else if (file_type == "csv") {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    log_message(paste("读取CSV文件:", file_path, "| 行数:", nrow(df)))
  } else {
    stop(paste("不支持的文件格式:", file_type), call.=FALSE)
  }
  
  # 检查必要列
  required_cols <- c("ID", "Description", "p.adjust", "Count")
  if (!all(required_cols %in% colnames(df))) {
    missing <- setdiff(required_cols, colnames(df))
    stop(paste("文件缺少必要列:", paste(missing, collapse=", ")), call.=FALSE)
  }
  
  # 确保数值列正确
  df$p.adjust <- as.numeric(df$p.adjust)
  df$Count <- as.numeric(df$Count)
  
  # KEGG ID标准化
  if (opt$type == "KEGG") {
    # 确保KEGG ID格式为5位数字
    df$ID <- ifelse(grepl("^\\d+$", df$ID), 
                   paste0(opt$kegg_organism, str_pad(df$ID, 5, pad = "0")),
                   df$ID)
    
    # 添加KEGG通路链接
    df$Link <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?", df$ID)
  }
  
  # 添加分析来源标记
  df$Analysis <- analysis_name
  df$logP <- -log10(df$p.adjust)
  
  return(df)
}

# ---------------------------
# 2. 数据准备和预处理
# ---------------------------
res1 <- read_enrich_results(opt$analysis1, "Analysis1")
res2 <- read_enrich_results(opt$analysis2, "Analysis2")

# 筛选显著通路
sig_res1 <- subset(res1, p.adjust < opt$pvalue)
sig_res2 <- subset(res2, p.adjust < opt$pvalue)

log_message(paste("分析1显著通路:", nrow(sig_res1)))
log_message(paste("分析2显著通路:", nrow(sig_res2)))

# ---------------------------
# 3. 相同通路分析（交集）
# ---------------------------
common_pathways <- intersect(sig_res1$ID, sig_res2$ID)
log_message(paste("共同通路数量:", length(common_pathways)))

if (length(common_pathways) > 0) {
  merged_common <- merge(
    sig_res1[sig_res1$ID %in% common_pathways, ],
    sig_res2[sig_res2$ID %in% common_pathways, ],
    by = "ID",
    suffixes = c("_1", "_2")
  )
  
  # 添加Description列，选择Description_1（假设Description_1和Description_2相同）
  merged_common$Description <- merged_common$Description_1
  
  # 添加富集方向变化列
  merged_common$logP_diff <- -log10(merged_common$p.adjust_1) - (-log10(merged_common$p.adjust_2))
  merged_common$Enrichment_Change <- ifelse(merged_common$logP_diff > 0, "Increased", "Decreased")
} else {
  merged_common <- data.frame()
  log_message("警告: 未发现共同通路")
}

# ---------------------------
# 4. 差异通路分析
# ---------------------------
unique_res1 <- sig_res1[!sig_res1$ID %in% sig_res2$ID, ]
unique_res2 <- sig_res2[!sig_res2$ID %in% sig_res1$ID, ]

log_message(paste("分析1特有通路:", nrow(unique_res1)))
log_message(paste("分析2特有通路:", nrow(unique_res2)))

# ---------------------------
# 5. 创建综合比较结果
# ---------------------------

# 5.1 创建通路比较总表
create_comparison_table <- function() {
  # 合并所有通路
  all_pathways <- bind_rows(
    sig_res1 %>% mutate(Group = ifelse(ID %in% common_pathways, "Common", "Unique_Analysis1")),
    sig_res2 %>% mutate(Group = ifelse(ID %in% common_pathways, "Common", "Unique_Analysis2"))
  ) %>%
    distinct(ID, .keep_all = TRUE) %>%
    dplyr::select(ID, Description, Group)
  
  # 添加分析1数据
  comp_table <- merge(
    all_pathways,
    sig_res1[, c("ID", "p.adjust", "Count", "logP")],
    by = "ID", all.x = TRUE
  )
  colnames(comp_table)[colnames(comp_table) %in% c("p.adjust", "Count", "logP")] <- 
    c("p.adjust_1", "Count_1", "logP_1")
  
  # 添加分析2数据
  comp_table <- merge(
    comp_table,
    sig_res2[, c("ID", "p.adjust", "Count", "logP")],
    by = "ID", all.x = TRUE
  )
  colnames(comp_table)[colnames(comp_table) %in% c("p.adjust", "Count", "logP")] <- 
    c("p.adjust_2", "Count_2", "logP_2")
  
  # 计算富集变化
  comp_table$logP_diff <- ifelse(
    !is.na(comp_table$logP_1) & !is.na(comp_table$logP_2),
    comp_table$logP_1 - comp_table$logP_2,
    NA
  )
  
  comp_table$Enrichment_Change <- case_when(
    comp_table$logP_diff > 1 ~ "Strongly Increased",
    comp_table$logP_diff > 0.5 ~ "Moderately Increased",
    comp_table$logP_diff > 0 ~ "Slightly Increased",
    comp_table$logP_diff < -1 ~ "Strongly Decreased",
    comp_table$logP_diff < -0.5 ~ "Moderately Decreased",
    comp_table$logP_diff < 0 ~ "Slightly Decreased",
    is.na(comp_table$logP_diff) & !is.na(comp_table$logP_1) ~ "Analysis1 Only",
    is.na(comp_table$logP_diff) & !is.na(comp_table$logP_2) ~ "Analysis2 Only",
    TRUE ~ "No Change"
  )
  
  comp_table$Significance_Level <- ifelse(
    !is.na(comp_table$p.adjust_1) & !is.na(comp_table$p.adjust_2), "Common",
    ifelse(
      !is.na(comp_table$p.adjust_1), "Unique_Analysis1",
      ifelse(
        !is.na(comp_table$p.adjust_2), "Unique_Analysis2", "Unknown"
      )
    )
  )
  
  # 添加KEGG链接
  if (opt$type == "KEGG") {
    comp_table$Link <- paste0("https://www.kegg.jp/kegg-bin/show_pathway?", comp_table$ID)
  }
  
  # 重新排序列
  col_order <- c("ID", "Description", "Group", "Significance_Level", "Enrichment_Change",
                 "p.adjust_1", "logP_1", "Count_1",
                 "p.adjust_2", "logP_2", "Count_2",
                 "logP_diff")
  
  if (opt$type == "KEGG") {
    col_order <- c(col_order, "Link")
  }
  
  other_cols <- setdiff(colnames(comp_table), col_order)
  comp_table <- comp_table[, c(col_order, other_cols)]
  
  return(comp_table)
}

# 5.2 创建富集变化统计表
create_change_stats <- function(comp_table) {
  if (nrow(comp_table) == 0) {
    return(data.frame())
  }
  
  # 仅处理共同通路
  common_df <- comp_table[comp_table$Group == "Common", ]
  
  if (nrow(common_df) == 0) {
    return(data.frame())
  }
  
  # 使用dplyr::count计算Enrichment_Change的频次
  change_stats <- common_df %>%
    dplyr::count(Enrichment_Change, name = "Count")
  
  # 添加统计值
  change_stats$Mean_logP_diff <- aggregate(
    logP_diff ~ Enrichment_Change,
    data = common_df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )$logP_diff
  
  change_stats$Max_logP_diff <- aggregate(
    logP_diff ~ Enrichment_Change,
    data = common_df,
    FUN = function(x) max(x, na.rm = TRUE)
  )$logP_diff
  
  change_stats$Min_logP_diff <- aggregate(
    logP_diff ~ Enrichment_Change,
    data = common_df,
    FUN = function(x) min(x, na.rm = TRUE)
  )$logP_diff
  
  # 按平均logP差异排序
  change_stats <- change_stats[order(change_stats$Mean_logP_diff, decreasing = TRUE), ]
  
  return(change_stats)
}

# 5.3 创建通路分类统计
create_category_stats <- function(comp_table) {
  if (nrow(comp_table) == 0) {
    return(data.frame())
  }
  
  # 使用dplyr计算Significance_Level的频次和统计值
  category_stats <- comp_table %>%
    dplyr::group_by(Significance_Level) %>%
    dplyr::summarise(
      Count = n(),
      Mean_Count_1 = mean(Count_1, na.rm = TRUE),
      Mean_Count_2 = mean(Count_2, na.rm = TRUE),
      Mean_logP_1 = mean(logP_1, na.rm = TRUE),
      Mean_logP_2 = mean(logP_2, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(category_stats)
}

# 生成比较结果
comparison_table <- create_comparison_table()
change_stats <- create_change_stats(comparison_table)
category_stats <- create_category_stats(comparison_table)

# ---------------------------
# 6. 可视化分析
# ---------------------------
save_plot <- function(plot_obj, filename, width=8, height=6) {
  # 保存为PDF
  pdf_file <- file.path(opt$output, paste0(tools::file_path_sans_ext(filename), ".pdf"))
  pdf(pdf_file, width=width, height=height)
  print(plot_obj)
  dev.off()
  
  # 保存为PNG
  png_file <- file.path(opt$output, paste0(tools::file_path_sans_ext(filename), ".png"))
  png(png_file, width=width*100, height=height*100, res=150)
  print(plot_obj)
  dev.off()
  
  log_message(paste("生成图表:", basename(pdf_file), "和", basename(png_file)))
}

# 6.1 Venn图展示通路分布
venn_plot <- ggvenn(
  list(Analysis1 = sig_res1$ID, Analysis2 = sig_res2$ID),
  fill_color = c("#1f77b4", "#ff7f0e"),
  stroke_size = 0.5,
  set_name_size = 4
) + ggtitle(paste(toupper(opt$type), "通路分布比较"))

save_plot(venn_plot, "01_pathway_venn_comparison.pdf")

# 6.2 富集变化分布图
if (nrow(comparison_table) > 0) {
  # 富集变化分布
  change_dist_plot <- ggplot(comparison_table, aes(x = Enrichment_Change, fill = Enrichment_Change)) +
    geom_bar() +
    scale_fill_brewer(palette = "RdYlBu", direction = -1) +
    labs(x = "富集变化方向", y = "通路数量", 
         title = paste(toupper(opt$type), "通路富集变化分布")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  save_plot(change_dist_plot, "07_enrichment_change_distribution.pdf")
  
  # 富集变化点图
  common_df <- comparison_table[comparison_table$Group == "Common", ]
  if (nrow(common_df) > 0) {
    change_point_plot <- ggplot(common_df, 
                              aes(x = logP_1, y = logP_2, color = Enrichment_Change)) +
      geom_point(aes(size = (Count_1 + Count_2)/2), alpha = 0.7) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      scale_color_brewer(palette = "RdYlBu", direction = -1) +
      labs(x = "分析1: -log10(FDR)", y = "分析2: -log10(FDR)",
           title = "共同通路富集变化点图",
           color = "富集变化", size = "平均基因数") +
      theme_minimal()
    
    save_plot(change_point_plot, "08_enrichment_change_dotplot.pdf")
  }
}

# 6.3 富集分数对比散点图（仅当有共同通路时）
if (nrow(merged_common) > 0) {
  scatter_plot <- ggplot(merged_common, aes(x = -log10(p.adjust_1), y = -log10(p.adjust_2))) +
    geom_point(aes(size = Count_1, color = logP_diff), alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
    labs(x = "分析1: -log10(FDR)", 
         y = "分析2: -log10(FDR)",
         title = paste("共同", toupper(opt$type), "通路富集程度比较"),
         color = "富集变化", size = "基因数") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  save_plot(scatter_plot, "02_pathway_scatter_comparison.pdf")
  
  # 6.4 通路相似性热图（仅GO）
  if (opt$type == "GO" && length(common_pathways) > 1) {
    tryCatch({
      hsGO <- godata(opt$species, ont = opt$go_ontology)
      go_sim <- mgoSim(common_pathways, common_pathways, semData = hsGO, measure = "Wang", combine = NULL)
      
      heatmap_plot <- pheatmap(go_sim, 
               color = colorRampPalette(c("white", "darkred"))(100),
               main = "共同GO通路语义相似性",
               silent = TRUE)
      
      save_plot(heatmap_plot[[4]], "03_GO_semantic_similarity.pdf", width=10, height=8)
    }, error = function(e) {
      log_message(paste("无法生成GO相似性热图:", e$message))
    })
  }
}

# 6.5 特有通路富集条形图
generate_barplot <- function(data, title, color) {
  if (nrow(data) == 0) return(NULL)
  
  # 按FDR排序并取前20个
  data <- data[order(data$p.adjust), ]
  top_terms <- head(data, min(20, nrow(data)))
  
  ggplot(top_terms, 
         aes(x = reorder(Description, -log10(p.adjust)), 
             y = -log10(p.adjust))) +
    geom_bar(stat = "identity", fill = color) +
    coord_flip() +
    labs(x = "", y = "-log10(FDR)", title = title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
}

bar_plot1 <- generate_barplot(unique_res1, paste("分析1特有", toupper(opt$type), "通路"), "steelblue")
if (!is.null(bar_plot1)) save_plot(bar_plot1, "04_unique_analysis1.pdf")

bar_plot2 <- generate_barplot(unique_res2, paste("分析2特有", toupper(opt$type), "通路"), "darkorange")
if (!is.null(bar_plot2)) save_plot(bar_plot2, "05_unique_analysis2.pdf")

# 6.6 富集分数差异图
if (nrow(merged_common) > 0) {
  diff_plot <- ggplot(merged_common, aes(x = reorder(Description, logP_diff), y = logP_diff)) +
    geom_bar(aes(fill = ifelse(logP_diff > 0, "Analysis1", "Analysis2")), stat = "identity") +
    scale_fill_manual(values = c("Analysis1" = "steelblue", "Analysis2" = "darkorange"),
                      name = "富集更强") +
    coord_flip() +
    labs(x = "", y = "富集分数差异(-log10FDR)", 
         title = paste("共同", toupper(opt$type), "通路富集程度差异")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  save_plot(diff_plot, "06_enrichment_difference.pdf")
}

# 6.7 KEGG通路可视化（仅KEGG）
if (opt$type == "KEGG" && length(common_pathways) > 0) {
  # 创建KEGG可视化目录
  kegg_dir <- file.path(opt$output, "KEGG_pathways")
  if (!dir.exists(kegg_dir)) dir.create(kegg_dir)
  
  # 下载并保存KEGG通路图
  log_message("开始生成KEGG通路图...")
  
  # 尝试获取通路名称
  tryCatch({
    kegg_info <- keggGet(common_pathways)
    pathway_names <- sapply(kegg_info, function(x) x$NAME)
    pathway_names <- gsub(" - .*", "", pathway_names)  # 简化名称
    names(pathway_names) <- common_pathways
  }, error = function(e) {
    pathway_names <- common_pathways
    log_message("无法获取KEGG通路名称，将使用ID")
  })
  
  # 生成通路图
  for (path_id in common_pathways) {
    tryCatch({
      # 创建通路图
      pathview(gene.data = NULL, 
               pathway.id = gsub("^.*?0", "", path_id),  # 移除物种前缀
               species = opt$kegg_organism,
               kegg.dir = kegg_dir,
               out.suffix = "pathview")
      
      # 重命名文件
      old_file <- paste0("path", gsub("^.*?0", "", path_id), ".pathview.png")
      new_file <- paste0(path_id, "_", gsub("/", "_", pathway_names[path_id]), ".png")
      
      if (file.exists(file.path(kegg_dir, old_file))) {
        file.rename(file.path(kegg_dir, old_file), 
                   file.path(kegg_dir, new_file))
        log_message(paste("生成KEGG通路图:", new_file))
      }
    }, error = function(e) {
      log_message(paste("无法生成KEGG通路图:", path_id, "| 错误:", e$message))
    })
  }
}

# ---------------------------
# 7. 结果导出
# ---------------------------
# 7.1 保存比较结果
log_message("保存比较结果文件...")

# 综合比较表
if (nrow(comparison_table) > 0) {
  out_file <- ifelse(opt$type == "GO", "all_GO_comparison.csv", "all_KEGG_comparison.csv")
  write.csv(comparison_table, file.path(opt$output, out_file), row.names = FALSE)
  log_message(paste("保存综合比较表:", out_file))
}

# 富集变化统计
if (nrow(change_stats) > 0) {
  out_file <- ifelse(opt$type == "GO", "GO_enrichment_change_stats.csv", "KEGG_enrichment_change_stats.csv")
  write.csv(change_stats, file.path(opt$output, out_file), row.names = FALSE)
  log_message(paste("保存富集变化统计表:", out_file))
}

# 通路分类统计
if (nrow(category_stats) > 0) {
  out_file <- ifelse(opt$type == "GO", "GO_category_stats.csv", "KEGG_category_stats.csv")
  write.csv(category_stats, file.path(opt$output, out_file), row.names = FALSE)
  log_message(paste("保存通路分类统计表:", out_file))
}

# 7.2 保存原始比较数据
if (nrow(merged_common) > 0) {
  out_file <- ifelse(opt$type == "GO", "common_GO_pathways.csv", "common_KEGG_pathways.csv")
  write.csv(merged_common, file.path(opt$output, out_file), row.names = FALSE)
  log_message(paste("保存共同通路数据:", out_file))
}

if (nrow(unique_res1) > 0) {
  out_file <- ifelse(opt$type == "GO", "unique_analysis1_GO.csv", "unique_analysis1_KEGG.csv")
  write.csv(unique_res1, file.path(opt$output, out_file), row.names = FALSE)
  log_message(paste("保存分析1特有通路:", out_file))
}

if (nrow(unique_res2) > 0) {
  out_file <- ifelse(opt$type == "GO", "unique_analysis2_GO.csv", "unique_analysis2_KEGG.csv")
  write.csv(unique_res2, file.path(opt$output, out_file), row.names = FALSE)
  log_message(paste("保存分析2特有通路:", out_file))
}

# 7.3 生成HTML报告
html_report <- function(comp_table) {
  report_file <- file.path(opt$output, "comparison_report.html")
  
  # 基本报告结构
  html_content <- paste0(
    '<!DOCTYPE html>
    <html lang="en">
    <head>
      <meta charset="UTF-8">
      <title>富集分析比较报告</title>
      <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1, h2 { color: #2c3e50; }
        .summary { background-color: #f9f9f9; padding: 20px; border-radius: 8px; }
        .images { display: flex; flex-wrap: wrap; gap: 20px; }
        .image-container { flex: 1 1 300px; text-align: center; }
        img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }
        table { width: 100%; border-collapse: collapse; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .section { margin-top: 30px; }
        .downloads { margin-top: 20px; }
        .download-item { margin: 10px 0; }
      </style>
    </head>
    <body>
      <h1>富集分析比较报告</h1>
      <div class="summary">
        <h2>分析概览</h2>
        <p><strong>分析日期:</strong> ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>
        <p><strong>富集类型:</strong> ', toupper(opt$type), '</p>
        <p><strong>分析1文件:</strong> ', basename(opt$analysis1), '</p>
        <p><strong>分析2文件:</strong> ', basename(opt$analysis2), '</p>
        <p><strong>显著性阈值:</strong> ', opt$pvalue, '</p>
        <p><strong>分析1显著通路:</strong> ', nrow(sig_res1), '</p>
        <p><strong>分析2显著通路:</strong> ', nrow(sig_res2), '</p>
        <p><strong>共同通路:</strong> ', length(common_pathways), '</p>
        <p><strong>分析1特有通路:</strong> ', nrow(unique_res1), '</p>
        <p><strong>分析2特有通路:</strong> ', nrow(unique_res2), '</p>
        
        <div class="downloads">
          <h3>结果下载</h3>')
  
  # 添加下载链接
  result_files <- list.files(opt$output, pattern = "\\.csv$", full.names = FALSE)
  for (file in result_files) {
    html_content <- paste0(html_content, '
          <div class="download-item">
            <a href="', file, '">', file, '</a>
          </div>')
  }
  
  html_content <- paste0(html_content, '
        </div>
      </div>')
  
  # 添加图像部分
  html_content <- paste0(html_content, '
      <div class="section">
        <h2>分析结果可视化</h2>
        <div class="images">')
  
  # 添加所有生成的图像
  img_files <- list.files(opt$output, pattern = "\\.png$", full.names = FALSE)
  for (img in img_files) {
    if (!grepl("pathview", img)) {  # 排除KEGG通路图
      html_content <- paste0(html_content, '
          <div class="image-container">
            <h3>', tools::file_path_sans_ext(img), '</h3>
            <img src="', img, '" alt="', img, '">
          </div>')
    }
  }
  
  html_content <- paste0(html_content, '
        </div>
      </div>')
  
  # 添加富集变化统计
  if (nrow(change_stats) > 0) {
    html_content <- paste0(html_content, '
      <div class="section">
        <h2>富集变化统计</h2>
        <table>
          <tr>
            <th>富集变化方向</th>
            <th>通路数量</th>
            <th>平均logP差异</th>
            <th>最大logP差异</th>
            <th>最小logP差异</th>
          </tr>')
    
    for (i in 1:nrow(change_stats)) {
      row <- change_stats[i, ]
      html_content <- paste0(html_content, '
          <tr>
            <td>', row$Enrichment_Change, '</td>
            <td>', row$Count, '</td>
            <td>', signif(row$Mean_logP_diff, 3), '</td>
            <td>', signif(row$Max_logP_diff, 3), '</td>
            <td>', signif(row$Min_logP_diff, 3), '</td>
          </tr>')
    }
    
    html_content <- paste0(html_content, '
        </table>
      </div>')
  }
  
  # 添加通路分类统计
  if (nrow(category_stats) > 0) {
    html_content <- paste0(html_content, '
      <div class="section">
        <h2>通路分类统计</h2>
        <table>
          <tr>
            <th>通路类别</th>
            <th>数量</th>
            <th>平均基因数(分析1)</th>
            <th>平均基因数(分析2)</th>
            <th>平均logP(分析1)</th>
            <th>平均logP(分析2)</th>
          </tr>')
    
    for (i in 1:nrow(category_stats)) {
      row <- category_stats[i, ]
      html_content <- paste0(html_content, '
          <tr>
            <td>', row$Significance_Level, '</td>
            <td>', row$Count, '</td>
            <td>', signif(row$Mean_Count_1, 3), '</td>
            <td>', signif(row$Mean_Count_2, 3), '</td>
            <td>', signif(row$Mean_logP_1, 3), '</td>
            <td>', signif(row$Mean_logP_2, 3), '</td>
          </tr>')
    }
    
    html_content <- paste0(html_content, '
        </table>
      </div>')
  }
  
  # 结束HTML
  html_content <- paste0(html_content, '
    </body>
    </html>')
  
  writeLines(html_content, report_file)
  log_message(paste("生成HTML报告:", basename(report_file)))
}

# 生成HTML报告
html_report(comparison_table)

log_message("分析完成！所有结果已保存到输出目录")
