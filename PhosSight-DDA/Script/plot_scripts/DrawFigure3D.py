import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print("Starting DrawFigure3D_v2.py...")

# Set style
plt.style.use('default')
sns.set_palette("husl")

# Data paths
all_features_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features.Localization.entropy.Label.txt'
psm_pga_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/PGA/psm_level/pga-peptideSummary.txt'
percolator_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2/DeepRescore2.psms.txt'
percolator_without_dlfeatures_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2/DeepRescore2WithoutDLFeatures.psms.txt'
deeprescore2_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/results_v1/DeepRescore2Results.Label.txt'

# Read data
all_features = pd.read_csv(all_features_path, sep='\t')
psm_pga_results = pd.read_csv(psm_pga_results_path, sep='\t')[['index', 'peptide', 'evalue']]
psm_pga_results.rename(columns={'index': 'Title'}, inplace=True)

data = pd.merge(psm_pga_results, all_features, on='Title', how='left')
data = data[data['PhosphoLabel'] == 1]

# Read percolator results
percolator_deep = pd.read_csv(percolator_results_path, sep='\t')
percolator_rescore = pd.read_csv(percolator_without_dlfeatures_path, sep='\t')
res_all = pd.read_csv(deeprescore2_results_path, sep='\t')

# Filter percolator results by q-value <= 0.01
percolator_deep = percolator_deep[percolator_deep['q-value'] <= 0.01]
percolator_rescore = percolator_rescore[percolator_rescore['q-value'] <= 0.01]

# Method definitions
# M1: PhosphoRS
M1 = data[data['phosphors_Sequence_Label'] == True].copy()
M1 = M1.sort_values(['PhosphoRS_IsoformProbability', 'PhosphoRS_IsoformScore'], ascending=False)

# M4: PhosphoRS + AutoRT + pDeep3
M4 = data[data['autort_pdeep_Sequence_Label'] == True].copy()
M4 = M4.sort_values(['autort_pDeep_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

# M5: PhosphoRS + PhosSight
M5 = data[data['phosSight_Sequence_Label'] == True].copy()
M5 = M5.sort_values(['PhosSightProb', 'PhosphoRS_IsoformScore'], ascending=False)

# M6: PhosphoRS + AutoRT + pDeep3 + PhosSight
M6 = data[data['autort_pdeep_phosSight_Sequence_Label'] == True].copy()
M6 = M6.sort_values(['autort_pDeep_phosSight_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

# M7: PhosphoRS + Rescore (traditional percolator)
M7 = res_all[res_all['Title'].isin(percolator_rescore['PSMId'])].copy()
M7 = M7[M7['phosphors_Sequence_Label'] == True]
M7 = M7.sort_values(['PhosphoRS_IsoformProbability', 'PhosphoRS_IsoformScore'], ascending=False)

# M8: PhosphoRS + DeepRescore (deep percolator)
M8 = res_all[res_all['Title'].isin(percolator_deep['PSMId'])].copy()
M8 = M8[M8['phosphors_Sequence_Label'] == True]
M8 = M8.sort_values(['PhosphoRS_IsoformProbability', 'PhosphoRS_IsoformScore'], ascending=False)

# M9: Method4 + DeepRescore
M9 = res_all[res_all['Title'].isin(percolator_deep['PSMId'])].copy()
M9 = M9[M9['autort_pdeep_Sequence_Label'] == True]
M9 = M9.sort_values(['autort_pDeep_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

# M10: Method6 + DeepRescore
M10 = data[data['Title'].isin(percolator_deep['PSMId'])].copy()
M10 = M10[M10['autort_pdeep_phosSight_Sequence_Label'] == True]
M10 = M10.sort_values(['autort_pDeep_phosSight_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

def compute_flr(df, prob_col, site_label_col, sequence_label_col, cutoff, use_sequence_label=True):
    """Compute FLR for a given method and probability cutoff"""
    if len(df) == 0:
        return np.nan
    
    # Filter by probability threshold first
    df2 = df[df[prob_col].notna() & (df[prob_col] >= cutoff)]
    if len(df2) == 0:
        return np.nan
    
    # Calculate FLR using the same logic as the reference script
    tl = 0
    fl = 0
    
    for i in range(len(df2)):
        if use_sequence_label and sequence_label_col is not None:
            # For methods that use both sequence and site labels
            seq_true = (df2.iloc[i][sequence_label_col] == True or 
                       df2.iloc[i][sequence_label_col] == 'TRUE')
            site_true = (df2.iloc[i][site_label_col] == 'TRUE')
            if seq_true and site_true:
                tl += 1
            else:
                fl += 1
        else:
            # Method7: only use site_label
            if df2.iloc[i][site_label_col] == 'TRUE':
                tl += 1
            else:
                fl += 1
    
    if (tl + fl) == 0:
        return np.nan
    return fl / (tl + fl)

# Thresholds
cuts = [0.25, 0.5, 0.75, 0.99]

# Method specifications (renamed and filtered)
spec = [
    {'name': 'Method1', 'df': M1, 'prob': 'PhosphoRS_IsoformProbability', 
     'site': 'phosphors_Site_Label', 'seq': 'phosphors_Sequence_Label', 
     'color': '#1f77b4', 'use_seq': True},
    {'name': 'Method2', 'df': M9, 'prob': 'autort_pDeep_Prob', 
     'site': 'autort_pdeep_Site_Label', 'seq': 'autort_pdeep_Sequence_Label', 
     'color': '#d62728', 'use_seq': True},  # Method9 -> Method2
    {'name': 'Method3', 'df': M7, 'prob': 'PhosphoRS_IsoformProbability', 
     'site': 'phosphors_Site_Label', 'seq': None, 
     'color': '#8c564b', 'use_seq': False},  # Method7 -> Method3
    {'name': 'Method4', 'df': M5, 'prob': 'PhosSightProb', 
     'site': 'phosSight_Site_Label', 'seq': 'phosSight_Sequence_Label', 
     'color': '#2ca02c', 'use_seq': True},  # Method5 -> Method4
    {'name': 'Method5', 'df': M10, 'prob': 'autort_pDeep_phosSight_Prob', 
     'site': 'autort_pdeep_phosSight_Site_Label', 'seq': 'autort_pdeep_phosSight_Sequence_Label', 
     'color': '#17becf', 'use_seq': True}  # Method10 -> Method5
]

# Compute FLR per method per cutoff
plot_data = []
for s in spec:
    flr_values = []
    for cutoff in cuts:
        flr = compute_flr(s['df'], s['prob'], s['site'], s['seq'], cutoff, s['use_seq'])
        flr_values.append(flr)
    
    for cutoff, flr in zip(cuts, flr_values):
        plot_data.append({
            'Method': s['name'],
            'cutoff': f'>{cutoff}',
            'FLR': flr * 100,  # Convert to percentage
            'color': s['color']
        })

# Create DataFrame for plotting
df_plot = pd.DataFrame(plot_data)

# Set method order (renamed methods)
method_order = ['Method1', 'Method2', 'Method3', 'Method4', 'Method5']
df_plot['Method'] = pd.Categorical(df_plot['Method'], categories=method_order, ordered=True)

# Set cutoff order
cutoff_order = [f'>{cutoff}' for cutoff in cuts]
df_plot['cutoff'] = pd.Categorical(df_plot['cutoff'], categories=cutoff_order, ordered=True)

# Create the plot
plt.figure(figsize=(11/2.54, 10/2.54))  # Convert cm to inches - same as DrawFigure3B_v2.py

# Define specific colors for each method (renamed methods)
method_colors = {
    'Method1': '#1f77b4',    # Blue
    'Method2': '#d62728',    # Purple (was Method9)
    'Method3': '#8c564b',    # Brown (was Method7)
    'Method4': '#2ca02c',    # Green (was Method5)
    'Method5': '#17becf'     # Cyan (was Method10)
}

# Plot lines for each method
for method in method_order:
    method_data = df_plot[df_plot['Method'] == method]
    if not method_data.empty:
        plt.plot(range(len(method_data)), method_data['FLR'], 
                label=method, color=method_colors[method], linewidth=1.5, marker='o', markersize=2)

# Customize plot
plt.xlabel('Site probability')
plt.ylabel('FLR(%)')
plt.ylim(0, 6)
plt.yticks(range(7), range(7))

# Set x-axis labels
plt.xticks(range(len(cutoff_order)), cutoff_order, rotation=45, ha='right')

# Add horizontal line at FLR = 1%
plt.axhline(y=1, color='grey', linestyle='--', alpha=0.7)

# Customize appearance
plt.grid(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(True)
plt.gca().spines['bottom'].set_visible(True)

# Add arrow to axes (same as DrawFigure3B_v2.py)
# ax = plt.gca()
# ax.annotate('', xy=(1, 0), xytext=(0, 0), 
#            arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
# ax.annotate('', xy=(0, 1), xytext=(0, 0), 
#            arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))

# Add legend
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()

# Save plot
output_path = 'Figure3D_v2.svg'
plt.savefig(output_path, format='svg', bbox_inches='tight', dpi=300)
print(f'Plot saved to: {output_path}')

plt.show()