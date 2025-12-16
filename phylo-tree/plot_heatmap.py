#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
脚本名: plot_heatmap.py
功能:  画不同clade结果相似度的heatmap图
作者:  adolf1
邮箱:  wangjyafk@126.com
创建:  2025-06-01
更新:  2025-12-16  adolf1  增加head说明
依赖:  Python≥3.8
用法:  python plot_heatmap.py --input in.txt --output out.txt
仓库:  https://github.com/AdolfFK/script_library
==============================================================================
"""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def draw_heatmap(file_path, output_pdf):
    # 读取数据
    df = pd.read_csv(file_path, delimiter='\t').set_index("matrix")
    
    # 创建一个与 df 同形状的 mask 矩阵，用于遮蔽上三角
    mask = np.triu(np.ones_like(df, dtype=bool))

    # 设置画布大小
    plt.figure(figsize=(8, 6))

    # 使用 seaborn 绘制热图，并应用 mask 遮罩
    sns.heatmap(
        df,
        annot=True,            # 显示数值
        fmt=".2f",             # 数值格式
        cmap="Blues",        # 颜色映射
        mask=mask,             # 只显示下三角
        cbar_kws={"label": "Similarity (%)"}
    )

    # 设置标题和坐标轴标签
    plt.title("Lower Triangle Similarity Heatmap")
    plt.xlabel("Clades")
    plt.ylabel("Clades")

    # 保存为 PDF
    with PdfPages(output_pdf) as pdf:
        pdf.savefig()
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Draw a heatmap from a TSV file and save it as a PDF.")
    parser.add_argument('--input', dest='file_path', required=True, help='Path to the input TSV file')
    parser.add_argument('--output', dest='output_pdf', required=True, help='Path to save the output PDF')

    args = parser.parse_args()
    draw_heatmap(args.file_path, args.output_pdf)

if __name__ == "__main__":
    main()


