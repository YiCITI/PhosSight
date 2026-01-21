#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
绘制非磷酸化数据集上三大类模型的AUC对比柱状图
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 标准样式
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Liberation Sans'],
    'font.size': 7,
    'axes.titlesize': 7,
    'axes.labelsize': 7,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.titlesize': 7,
    'axes.linewidth': 1,
    'grid.linewidth': 0.8,
    'lines.linewidth': 1,
    'patch.linewidth': 1,
    'xtick.major.width': 1,
    'ytick.major.width': 1,
    'xtick.minor.width': 1,
    'ytick.minor.width': 1,
    'axes.unicode_minus': False
})

def load_data():
    """加载数据"""
    df = pd.read_csv('/data0/wangb/cd/duibi0826/0826comparison/model_performance_comparison.csv')
    
    # 处理数据集名称
    df['Dataset'] = df['Dataset'].str.replace('.csv', '')
    df['Dataset'] = df['Dataset'].str.replace('DeepDetect_', '')
    df['Dataset'] = df['Dataset'].str.capitalize()
    
    return df

def ensure_dir(path):
    import os
    os.makedirs(path, exist_ok=True)

def plot_nonphos_group(df, datasets, models, labels, colors, output_svg, fig_height=6, legend_y=1.00, legend_x=0.5):
    import os
    ensure_dir(os.path.dirname(output_svg))

    subset = df[df['Dataset'].isin(datasets)].copy()

    fig, ax = plt.subplots(figsize=(12, fig_height))

    width = 0.08
    dataset_gap = 0.45
    x = np.arange(len(datasets)) * dataset_gap

    num_models = len(models)
    offsets = (np.arange(num_models) - (num_models - 1) / 2.0) * width

    for i, (model, label, color) in enumerate(zip(models, labels, colors)):
        if model in subset.columns:
            values = subset[model].values
            ax.bar(x + offsets[i], values, width, label=label, color=color,
                   alpha=0.9, edgecolor='white', linewidth=0.6)

    # 为每个数据集组的最佳模型添加红色五角星标记
    star_color = '#d62728'
    star_size = 320
    y_margin = 0.01
    max_star_y = 0.0
    for j in range(len(datasets)):
        vals = []
        idxs = []
        for i, model in enumerate(models):
            if model in subset.columns:
                vals.append(subset[model].values[j])
                idxs.append(i)
        if not vals:
            continue
        max_k = int(np.argmax(vals))
        best_i = idxs[max_k]
        best_val = float(vals[max_k])
        star_x = x[j] + offsets[best_i]
        star_y = best_val + y_margin
        max_star_y = max(max_star_y, star_y)
        ax.scatter([star_x], [star_y], marker='*', s=star_size, c=star_color,
                   edgecolors='none', zorder=5, clip_on=False)

    ax.set_xlabel('Datasets')
    ax.set_ylabel('AUC Score')
    ax.set_xticks(x)
    ax.set_xticklabels(datasets)
    # 根据是否为finetuned图调整y轴范围，finetuned图缩小纵坐标间距
    if 'finetuned' in output_svg:
        ax.set_ylim(0.75, 1.0)  # 缩小范围，让刻度更密集
        ax.set_yticks(np.arange(0.75, 1.01, 0.05))  # 设置更密集的刻度
    else:
        ax.set_ylim(0.7, 1.0)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # 图例位置可调整，避免遮挡柱子
    ax.legend(ncol=2, frameon=False, fancybox=False, shadow=False, borderaxespad=0.2,
              loc='upper center', bbox_to_anchor=(legend_x, legend_y))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    # 确保星标不被裁剪
    ymin, ymax = ax.get_ylim()
    if max_star_y > ymax:
        ax.set_ylim(ymin, max_star_y + 0.02)

    plt.tight_layout()
    plt.savefig(output_svg, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)

def create_non_phosphorylation_plots():
    df = load_data()
    datasets = ['Ecoli', 'Human', 'Mouse', 'Yeast']

    original_models = ['DeepDetect', 'Original_DeepDetect', 'PFly', 'PhosSight']
    original_labels = ['DeepDetect', 'Original DeepDetect', 'PFly', 'PhosDetect']
    # 顶刊配色：其他三个模型使用相近的低饱和度蓝绿色系，PhosDetect使用低饱和度红色系突出
    original_colors = ['#B8D4E3', '#B0D4D4', '#A8D4C8', '#D4A5A5']

    finetuned_models = ['Fine-tuned_DeepDetect', 'Fine-tuned_Original_DeepDetect', 'Fine-tuned_PFly', 'Fine-tuned_PhosSight']
    finetuned_labels = ['Fine-tuned DeepDetect', 'Fine-tuned Original DeepDetect', 'Fine-tuned PFly', 'Fine-tuned PhosDetect']
    # Finetuned版本颜色稍深但仍保持淡色调，PhosDetect使用稍深的低饱和度红色
    finetuned_colors = ['#9CC5D8', '#98C5C5', '#90C5B0', '#C59595']

    outdir = '/data0/wangb/cd/plottu'
    plot_nonphos_group(
        df, datasets,
        original_models, original_labels, original_colors,
        os.path.join(outdir, 'non_phosphorylation_original_auc_comparison.svg')
    )
    plot_nonphos_group(
        df, datasets,
        finetuned_models, finetuned_labels, finetuned_colors,
        os.path.join(outdir, 'non_phosphorylation_finetuned_auc_comparison.svg'),
        fig_height=6,  # 压缩高度，通过缩小纵坐标间距来实现
        legend_y=1.02,  # 图例位置往上移，避免和柱子重合
        legend_x=0.65  # 图例位置往右移，避免和柱子重合
    )

def main():
    """主函数"""
    print("开始创建非磷酸化数据集AUC对比图 (原始 vs Finetuned 分开)...")
    create_non_phosphorylation_plots()
    print("\n已生成 2 个SVG:")
    print("- non_phosphorylation_original_auc_comparison.svg")
    print("- non_phosphorylation_finetuned_auc_comparison.svg")

if __name__ == "__main__":
    main()
