#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

ORIGINAL_CSV = '/data0/wangb/cd/duibi0826/0826comparison/evaluation_results/phosight_v2_summary_results.csv'
FINETUNED_CSV = '/data0/wangb/cd/duibi0826/0826comparison/evaluation_results/phosight_finetuned_summary_results.csv'
OUTDIR = '/data0/wangb/cd/plottu'


def load_merged():
    df_orig = pd.read_csv(ORIGINAL_CSV)
    df_ft = pd.read_csv(FINETUNED_CSV)

    df_orig = df_orig[df_orig['Model'] == 'PhosSight V2'].copy()
    df_ft = df_ft[df_ft['Model'] == 'PhosSight_finetuned'].copy()

    cols = ['Dataset', 'Dataset_Type', 'AUC', 'Accuracy', 'Precision', 'Recall', 'F1_Score']
    df_orig = df_orig[cols]
    df_ft = df_ft[cols]
    df = pd.merge(df_orig, df_ft, on=['Dataset', 'Dataset_Type'], suffixes=('_orig', '_ft'))
    return df


def plot_dataset_triplet(df_row: pd.Series, outpath: str, title=None):
    metrics = ['Accuracy', 'Precision', 'F1_Score']
    x = np.arange(len(metrics))
    # 减小柱宽和间距，使三张图总宽度为14
    width = 0.18
    # 减小x轴间距
    x_spacing = 0.7  # 原来默认是1.0

    orig_vals = [df_row[f'{m}_orig'] for m in metrics]
    ft_vals = [df_row[f'{m}_ft'] for m in metrics]

    # 调整图片宽度：三张图总宽度为14，每张图约4.67，高度与amino_acid_type_performance.svg一致为5
    fig, ax = plt.subplots(figsize=(4.67, 5))
    # 使用低饱和度配色（与AUC对比图一致）
    # Original用淡蓝色，Finetuned用红色系突出（与最优模型色系一致）
    ax.bar(x * x_spacing - width/2, orig_vals, width, label='Original', 
           color='#B8D4E3', edgecolor='white', linewidth=0.5, alpha=0.85)
    ax.bar(x * x_spacing + width/2, ft_vals, width, label='Finetuned', 
           color='#D4A5A5', edgecolor='white', linewidth=0.5, alpha=0.85)

    ax.set_xticks(x * x_spacing)
    ax.set_xticklabels(['Accuracy', 'Precision', 'F1-Score'])
    ax.set_ylabel('Score')
    ax.set_ylim(0.0, 1.0)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # 图例放在右上角
    ax.legend(ncol=1, frameon=False, loc='upper right')

    # 顶部标出数据集名称
    if title is None:
        ds_name = str(df_row['Dataset']).replace('DeepRescore2_', '').replace('DeepDetect_', '')
    else:
        ds_name = title
    ax.set_title(ds_name, pad=8)

    # 轴样式
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)


def plot_combined_phosphorylation(df):
    """绘制三个磷酸化数据集汇总图"""
    metrics = ['Accuracy', 'Precision', 'F1_Score']
    phos_keys = ['DeepRescore2_HCC', 'DeepRescore2_label_free', 'DeepRescore2_UCEC']
    dataset_names = ['HCC', 'Label-free', 'UCEC']
    
    x = np.arange(len(metrics))
    dataset_gap = 0.25
    width = 0.06
    
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # 为每个数据集计算位置偏移
    for i, (key, ds_name) in enumerate(zip(phos_keys, dataset_names)):
        row = df[df['Dataset'] == key].iloc[0]
        base_pos = x + i * dataset_gap
        
        orig_vals = [row[f'{m}_orig'] for m in metrics]
        ft_vals = [row[f'{m}_ft'] for m in metrics]
        
        # 蓝色配紫色对比
        ax.bar(base_pos - width/2, orig_vals, width, 
               label='Original' if i == 0 else '', 
               color='#1f77b4', edgecolor='white', linewidth=0.6, alpha=0.8)
        ax.bar(base_pos + width/2, ft_vals, width, 
               label='Finetuned' if i == 0 else '', 
               color='#9467bd', edgecolor='white', linewidth=0.6, alpha=0.8)
        
        # 在每组柱子上方添加数据集名称
        for j, metric in enumerate(metrics):
            ax.text(base_pos[j], 1.05, ds_name, ha='center', va='bottom', 
                   rotation=0)
    
    ax.set_xticks(x + dataset_gap * (len(phos_keys) - 1) / 2)
    ax.set_xticklabels(['Accuracy', 'Precision', 'F1-Score'])
    ax.set_ylabel('Score')
    ax.set_ylim(0.0, 1.15)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    # 图例在右上角
    ax.legend(ncol=1, frameon=False, loc='upper right')
    
    # 轴样式
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)
    
    plt.tight_layout()
    outsvg = os.path.join(OUTDIR, 'phossight_phos_combined_orig_vs_finetuned.svg')
    os.makedirs(os.path.dirname(outsvg), exist_ok=True)
    plt.savefig(outsvg, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f'Saved: {outsvg}')




def main():
    df = load_merged()
    # 仅三个磷酸化数据集
    phos_keys = ['DeepRescore2_HCC', 'DeepRescore2_label_free', 'DeepRescore2_UCEC']

    # 绘制三张单独的图
    for key in phos_keys:
        row = df[df['Dataset'] == key].iloc[0]
        safe_name = key.replace('DeepRescore2_', '')
        outsvg = os.path.join(OUTDIR, f'phossight_phos_{safe_name}_orig_vs_finetuned_triplet.svg')
        plot_dataset_triplet(row, outsvg)
        print(f'Saved: {outsvg}')
    
    # 绘制磷酸化数据集汇总图
    plot_combined_phosphorylation(df)
    
    # 绘制非磷酸化数据集Human的单独图
    human_key = 'DeepDetect_human'
    if human_key in df['Dataset'].values:
        row = df[df['Dataset'] == human_key].iloc[0]
        outsvg = os.path.join(OUTDIR, 'phossight_non_phos_human_orig_vs_finetuned_triplet.svg')
        plot_dataset_triplet(row, outsvg, title='Human')
        print(f'Saved: {outsvg}')


if __name__ == '__main__':
    main()


