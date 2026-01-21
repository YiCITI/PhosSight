#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy import stats
import seaborn as sns
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

def load_data(file_path):
    """加载数据"""
    print(f"📊 加载数据: {file_path}")
    
    if not os.path.exists(file_path):
        print(f"❌ 文件不存在: {file_path}")
        return None
    
    df = pd.read_csv(file_path)
    print(f"✅ 成功加载 {len(df)} 条数据")
    
    # 分离正负样本
    positive_df = df[df['label'] == 1]
    negative_df = df[df['label'] == 0]
    
    positive_scores = positive_df['detectability_score'].values
    negative_scores = negative_df['detectability_score'].values
    
    print(f"📈 数据统计:")
    print(f"   - 正样本数: {len(positive_scores)}")
    print(f"   - 负样本数: {len(negative_scores)}")
    print(f"   - 正样本平均分数: {np.mean(positive_scores):.4f}")
    print(f"   - 负样本平均分数: {np.mean(negative_scores):.4f}")
    
    return positive_scores, negative_scores

def calculate_p_value(positive_scores, negative_scores):
    """计算p值"""
    print(f"\n🔬 计算p值...")
    
    # 正态性检验
    _, p_positive = stats.shapiro(positive_scores)
    _, p_negative = stats.shapiro(negative_scores)
    
    # 方差齐性检验
    _, p_levene = stats.levene(positive_scores, negative_scores)
    
    # 选择适当的检验方法
    is_normal = p_positive > 0.05 and p_negative > 0.05
    is_equal_var = p_levene > 0.05
    
    if is_normal and is_equal_var:
        # 使用t检验
        _, p_value = stats.ttest_ind(positive_scores, negative_scores)
        test_name = "t-test"
    else:
        # 使用Mann-Whitney U检验
        _, p_value = stats.mannwhitneyu(positive_scores, negative_scores, alternative='two-sided')
        test_name = "Mann-Whitney U test"
    
    print(f"使用检验方法: {test_name}")
    print(f"p值: {p_value:.6f}")
    
    return p_value, test_name

def add_jitter(data, jitter_amount=0.1):
    """为散点添加抖动，避免重叠"""
    return data + np.random.normal(0, jitter_amount, len(data))

def create_nature_boxplot(positive_scores, negative_scores, p_value, output_file="nature_boxplot_10000.png"):
    """创建Nature期刊标准箱线图"""
    print(f"\n🎨 创建Nature标准箱线图: {output_file}")
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(6, 8))
    
    # Nature期刊标准配色方案
    colors = ['#1f77b4', '#ff7f0e']  # 专业蓝色和橙色
    edge_colors = ['#0d47a1', '#e65100']  # 深色边框
    
    # 准备数据
    data_to_plot = [positive_scores, negative_scores]
    labels = ['Positive', 'Negative']
    
    # 创建箱线图
    bp = ax.boxplot(data_to_plot, labels=labels, patch_artist=True,
                    boxprops=dict(facecolor='white', alpha=0.8, linewidth=1.2),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(color='black', linewidth=1.2),
                    capprops=dict(color='black', linewidth=1.2),
                    flierprops=dict(marker='o', markerfacecolor='gray', 
                                  markeredgecolor='black', markersize=3, alpha=0.6))
    
    # 设置箱线图颜色
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.7)
        patch.set_edgecolor(edge_colors[i])
        patch.set_linewidth(1.2)
    
    # 添加散点图 - 使用抖动避免重叠
    pos_jittered = add_jitter(np.ones(len(positive_scores)), 0.15)
    neg_jittered = add_jitter(np.ones(len(negative_scores)) * 2, 0.15)
    
    # 绘制散点
    ax.scatter(pos_jittered, positive_scores, color=colors[0], alpha=0.4, s=8, edgecolors='white', linewidth=0.5)
    ax.scatter(neg_jittered, negative_scores, color=colors[1], alpha=0.4, s=8, edgecolors='white', linewidth=0.5)
    
    # 添加均值线
    pos_mean = np.mean(positive_scores)
    neg_mean = np.mean(negative_scores)
    ax.axhline(y=pos_mean, xmin=0.7, xmax=1.3, color=colors[0], linestyle='--', alpha=0.8, linewidth=1.5)
    ax.axhline(y=neg_mean, xmin=1.7, xmax=2.3, color=colors[1], linestyle='--', alpha=0.8, linewidth=1.5)
    
    # 统计信息已移除，保持图表简洁
    
    # 添加p值标注
    if p_value < 0.001:
        p_text = "p < 0.001***"
    elif p_value < 0.01:
        p_text = f"p = {p_value:.3f}**"
    elif p_value < 0.05:
        p_text = f"p = {p_value:.3f}*"
    else:
        p_text = f"p = {p_value:.3f}"
    
    # 在图形顶部添加p值
    ax.text(0.5, 0.95, p_text, transform=ax.transAxes, ha='center', va='top', 
            color='black',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", 
                     alpha=0.9, edgecolor='gray', linewidth=0.5))
    
    # 设置标签
    ax.set_ylabel('Detectability Score')
    ax.set_xlabel('Sample Type')
    ax.set_title('Model Performance on Synthetic Dataset\n(10,000 samples)', 
                pad=15)
    
    # 设置网格
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # 美化坐标轴
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    
    # 设置刻度
    ax.tick_params(axis='both', which='major', labelsize=9, length=4, width=0.8)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存高质量图片
    plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"✅ Nature标准箱线图已保存: {output_file}")

def create_nature_violin_plot(positive_scores, negative_scores, p_value, output_file="nature_violin_10000.png"):
    """基于组合图violin plot代码，增强差异表现"""
    print(f"\n🎨 创建增强差异表现的小提琴图: {output_file}")
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # 使用组合图中的配色方案
    colors = ['#1f77b4', '#ff7f0e']  # 蓝色和橙色
    
    # 准备散点数据
    pos_jittered = add_jitter(np.ones(len(positive_scores)), 0.15)
    neg_jittered = add_jitter(np.ones(len(negative_scores)) * 2, 0.15)
    
    # 使用seaborn创建小提琴图，添加cut参数增强差异
    import seaborn as sns
    
    # 准备数据
    data = []
    labels = []
    for scores, label in [(positive_scores, 'Positive'), (negative_scores, 'Negative')]:
        data.extend(scores)
        labels.extend([label] * len(scores))
    
    df_plot = pd.DataFrame({'Score': data, 'Type': labels})
    
    # 创建小提琴图，使用cut参数增加截断效果
    sns.violinplot(data=df_plot, x='Type', y='Score', ax=ax, 
                   inner='box',  # 显示小箱线图
                   cut=2,        # 增加两边的截断，增强视觉差异
                   palette=colors,
                   linewidth=1.2,
                   width=0.8)    # 调整宽度
    
    # 设置小提琴图颜色和样式（参考组合图）
    for i, violin in enumerate(ax.collections):
        if i < 2:  # 只处理前两个小提琴图
            violin.set_alpha(0.7)
            violin.set_edgecolor('black')
            violin.set_linewidth(0.8)
    
    # 添加散点（参考组合图样式）
    ax.scatter(pos_jittered, positive_scores, color=colors[0], alpha=0.3, s=6, 
               edgecolors='white', linewidth=0.3)
    ax.scatter(neg_jittered, negative_scores, color=colors[1], alpha=0.3, s=6, 
               edgecolors='white', linewidth=0.3)
    
    # 设置标签（参考组合图）
    ax.set_ylabel('Detectability Score')
    ax.set_xlabel('Sample Type')
    ax.set_title('Model Performance Distribution\n(10,000 samples)', 
                pad=15)
    
    # 设置网格（参考组合图）
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # 美化坐标轴（参考组合图）
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    
    # 设置刻度（参考组合图）
    ax.tick_params(axis='both', which='major', labelsize=9, length=4, width=0.8)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"✅ 增强差异表现的小提琴图已保存: {output_file}")

def create_enhanced_boxplot(positive_scores, negative_scores, p_value, output_dir=None, output_basename="enhanced_boxplot_10000"):
    """创建增强差异表现的箱线图，并导出为 PNG 和 SVG。

    如果提供 output_dir，则将结果保存至该目录；否则保存到当前工作目录。
    文件名基于 output_basename 生成。
    """
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    png_path = os.path.join(output_dir, f"{output_basename}.png")
    svg_path = os.path.join(output_dir, f"{output_basename}.svg")
    print(f"\n🎨 创建增强差异表现的箱线图: {png_path} 和 {svg_path}")
    
    # 创建图形
    # 画布尺寸：宽度5，高度5
    fig, ax = plt.subplots(figsize=(5, 5))
    
    # 使用组合图中的配色方案
    colors = ['#1f77b4', '#9467bd']  # 蓝色和紫色
    edge_colors = ['#0d47a1', '#7b3f98']  # 深色边框
    
    # 先创建小提琴图
    import seaborn as sns
    
    # 准备数据
    data = []
    labels = []
    for scores, label in [(positive_scores, 'Positive'), (negative_scores, 'Negative')]:
        data.extend(scores)
        labels.extend([label] * len(scores))
    
    df_plot = pd.DataFrame({'Score': data, 'Type': labels})
    
    # 创建小提琴图，不显示内部元素
    sns.violinplot(data=df_plot, x='Type', y='Score', ax=ax, 
                   inner=None,  # 不显示内部元素，我们自己添加
                   cut=0,        # 不截断，显示完整分布
                   palette=colors,
                   linewidth=0.8,
                   width=0.55)
    
    # 设置小提琴图颜色和样式
    for i, violin in enumerate(ax.collections):
        if i < 2:  # 只处理前两个小提琴图
            violin.set_alpha(0.2)  # 更轻的填充
            violin.set_edgecolor(edge_colors[i])
            violin.set_linewidth(0.8)
    
    # 手动添加简洁的boxplot（无异常值和散点），分别控制左右宽度
    # 获取seaborn violin plot的x轴位置
    violin_positions = [0, 1]  # seaborn默认从0开始
    # 左侧（蓝色）
    bp_left = ax.boxplot([positive_scores], positions=[violin_positions[0]], patch_artist=True,
                         widths=0.18,
                         boxprops=dict(facecolor='white', alpha=0.72, linewidth=0.7),
                         medianprops=dict(color='black', linewidth=1.0),
                         whiskerprops=dict(color='black', linewidth=0.7),
                         capprops=dict(color='black', linewidth=0.7),
                         flierprops=dict(marker='', markersize=0))
    for patch in bp_left['boxes']:
        patch.set_facecolor(colors[0])
        patch.set_alpha(0.6)
        patch.set_edgecolor(edge_colors[0])
        patch.set_linewidth(1.0)
    # 右侧（紫色）进一步变窄
    bp_right = ax.boxplot([negative_scores], positions=[violin_positions[1]], patch_artist=True,
                          widths=0.12,
                          boxprops=dict(facecolor='white', alpha=0.72, linewidth=0.65),
                          medianprops=dict(color='black', linewidth=0.9),
                          whiskerprops=dict(color='black', linewidth=0.65),
                          capprops=dict(color='black', linewidth=0.65),
                          flierprops=dict(marker='', markersize=0))
    for patch in bp_right['boxes']:
        patch.set_facecolor(colors[1])
        patch.set_alpha(0.6)
        patch.set_edgecolor(edge_colors[1])
        patch.set_linewidth(0.9)
    
    # 在图形顶部添加p值
    if p_value < 0.001:
        p_text = "p < 0.001***"
    elif p_value < 0.01:
        p_text = f"p = {p_value:.3f}**"
    elif p_value < 0.05:
        p_text = f"p = {p_value:.3f}*"
    else:
        p_text = f"p = {p_value:.3f}"
    # p值标注放在坐标系内，顶部与纵坐标最上方（y=1.0）对齐
    ax.text(0.5, 1.0, p_text, transform=ax.transAxes, ha='center', va='top',
            color='black',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      alpha=0.85, edgecolor='gray', linewidth=0.5))

    # 设置标签
    ax.set_ylabel('Detectability Score')
    # 删除横坐标标签
    # ax.set_xlabel('Sample Type', fontsize=10, fontweight='bold')
    # 去掉标题，减少视觉拥挤
    # ax.set_title('Model Performance Comparison\n(10,000 samples)', fontsize=12, fontweight='bold', pad=12)
    
    # 设置网格
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # 美化坐标轴
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    
    # 设置刻度
    ax.tick_params(axis='both', which='major', labelsize=9, length=4, width=0.8)
    
    # 确保横坐标显示为Positive和Negative（但不显示x轴名称）
    ax.set_xticklabels(['Positive', 'Negative'])
    
    # 留出顶部和四周边距，防止元素被裁切
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.06 * (ymax - ymin))
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    # 保存为 PNG 和 SVG
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.savefig(svg_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"✅ 增强差异表现的箱线图已保存: {png_path}")
    print(f"✅ SVG 矢量图已保存: {svg_path}")

def create_combined_plot(positive_scores, negative_scores, p_value, output_file="nature_combined_10000.png"):
    """创建组合图（箱线图+小提琴图）"""
    print(f"\n🎨 创建组合图: {output_file}")
    
    # 创建子图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Nature期刊标准配色
    colors = ['#1f77b4', '#ff7f0e']
    
    # 左图：箱线图
    data_to_plot = [positive_scores, negative_scores]
    bp = ax1.boxplot(data_to_plot, labels=['Positive', 'Negative'], patch_artist=True,
                    boxprops=dict(facecolor='white', alpha=0.8, linewidth=1.2),
                    medianprops=dict(color='black', linewidth=1.5),
                    whiskerprops=dict(color='black', linewidth=1.2),
                    capprops=dict(color='black', linewidth=1.2),
                    flierprops=dict(marker='o', markerfacecolor='gray', 
                                  markeredgecolor='black', markersize=3, alpha=0.6))
    
    for i, patch in enumerate(bp['boxes']):
        patch.set_facecolor(colors[i])
        patch.set_alpha(0.7)
        patch.set_edgecolor(colors[i])
        patch.set_linewidth(1.2)
    
    # 添加散点
    pos_jittered = add_jitter(np.ones(len(positive_scores)), 0.15)
    neg_jittered = add_jitter(np.ones(len(negative_scores)) * 2, 0.15)
    ax1.scatter(pos_jittered, positive_scores, color=colors[0], alpha=0.4, s=8, 
                edgecolors='white', linewidth=0.5)
    ax1.scatter(neg_jittered, negative_scores, color=colors[1], alpha=0.4, s=8, 
                edgecolors='white', linewidth=0.5)
    
    ax1.set_ylabel('Detectability Score')
    ax1.set_title('Box Plot')
    ax1.grid(True, alpha=0.3)
    
    # 右图：小提琴图 - 使用seaborn增强差异表现
    import seaborn as sns
    
    # 准备数据
    data = []
    labels = []
    for scores, label in [(positive_scores, 'Positive'), (negative_scores, 'Negative')]:
        data.extend(scores)
        labels.extend([label] * len(scores))
    
    df_plot = pd.DataFrame({'Score': data, 'Type': labels})
    
    # 创建小提琴图，使用cut参数增加截断效果
    sns.violinplot(data=df_plot, x='Type', y='Score', ax=ax2, 
                   inner='box',  # 显示小箱线图
                   cut=2,        # 增加两边的截断，增强视觉差异
                   palette=colors,
                   linewidth=1.2,
                   width=0.8)    # 调整宽度
    
    # 设置小提琴图颜色和样式
    for i, violin in enumerate(ax2.collections):
        if i < 2:  # 只处理前两个小提琴图
            violin.set_alpha(0.7)
            violin.set_edgecolor('black')
            violin.set_linewidth(0.8)
    
    # 添加散点
    ax2.scatter(pos_jittered, positive_scores, color=colors[0], alpha=0.3, s=6, 
                edgecolors='white', linewidth=0.3)
    ax2.scatter(neg_jittered, negative_scores, color=colors[1], alpha=0.3, s=6, 
                edgecolors='white', linewidth=0.3)
    
    ax2.set_xticks([1, 2])
    ax2.set_xticklabels(['Positive', 'Negative'])
    ax2.set_ylabel('Detectability Score')
    ax2.set_title('Violin Plot')
    ax2.grid(True, alpha=0.3)
    
    # 添加p值到箱线图（左图）
    if p_value < 0.001:
        p_text = "p < 0.001***"
    elif p_value < 0.01:
        p_text = f"p = {p_value:.3f}**"
    elif p_value < 0.05:
        p_text = f"p = {p_value:.3f}*"
    else:
        p_text = f"p = {p_value:.3f}"
    
    # 只在箱线图（左图）添加p值标签
    ax1.text(0.5, 0.95, p_text, transform=ax1.transAxes, 
            ha='center', va='top',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9, 
                     edgecolor='gray', linewidth=0.5))
    
    plt.suptitle('Model Performance on Synthetic Dataset (10,000 samples)')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"✅ 组合图已保存: {output_file}")

def main():
    """主函数"""
    print("📊 开始创建Nature期刊标准图表...")
    
    # 数据文件路径（使用脚本所在目录的绝对路径，避免因当前工作目录不同而找不到文件）
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_file = os.path.join(script_dir, "synthetic_dataset_10000_simple_with_scores.csv")
    
    # 加载数据
    result = load_data(data_file)
    if result is None:
        return
    
    positive_scores, negative_scores = result
    
    # 计算p值
    p_value, test_name = calculate_p_value(positive_scores, negative_scores)
    
    # 创建增强差异表现的箱线图（导出 PNG 和 SVG 到指定目录）
    output_dir = "/data0/wangb/cd/plottu"
    create_enhanced_boxplot(positive_scores, negative_scores, p_value, output_dir=output_dir, output_basename="enhanced_boxplot_10000")

    print(f"\n🎉 完成!")
    print(f"📁 输出目录: {output_dir}")
    print(f"   - enhanced_boxplot_10000.png")
    print(f"   - enhanced_boxplot_10000.svg")

if __name__ == "__main__":
    main()
