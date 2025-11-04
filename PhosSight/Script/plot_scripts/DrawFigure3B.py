import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print("Starting DrawFigure3B_v2.py...")
print("Loading data files...")

# Set style
plt.style.use('default')
sns.set_palette("husl")

# Data paths
all_features_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Features.Localization.entropy.Label.txt'
psm_pga_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/PGA/psm_level/pga-peptideSummary.txt'
percolator_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2/DeepRescore2.psms.txt'
percolator_without_dlfeatures_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/Percolator/DeepRescore2/DeepRescore2WithoutDLFeatures.psms.txt'
deeprescore2_with_phossight_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/PXD000138/DeepRescore2Results.Label.txt'
deeprescore2_results_path = 'E:/MyProject/PTM-Project/datasets/ExampleData1/results_v1/DeepRescore2Results.Label.txt'

# Read data
print("Reading all_features...")
all_features = pd.read_csv(all_features_path, sep='\t')
print("Reading psm_pga_results...")
psm_pga_results = pd.read_csv(psm_pga_results_path, sep='\t')[['index', 'peptide', 'evalue']]
psm_pga_results.rename(columns={'index': 'Title'}, inplace=True)
print("Merging data...")
data = pd.merge(psm_pga_results, all_features, on='Title', how='left')
data = data[data['PhosphoLabel'] == 1]
print("Data loading completed.")

# Read additional data
deeprescore2_percolator_results = pd.read_csv(percolator_results_path, sep='\t')
rescore_percolator_results = pd.read_csv(percolator_without_dlfeatures_path, sep='\t')
deeprescore2_results = pd.read_csv(deeprescore2_results_path, sep='\t')
deeprescore2_with_phossight_results = pd.read_csv(deeprescore2_with_phossight_results_path, sep='\t')

# Method1: PhosphoRS
Method1 = data[data['phosphors_Sequence_Label'] == True].copy()
Method1 = Method1.sort_values(['PhosphoRS_IsoformProbability', 'PhosphoRS_IsoformScore'], ascending=False)

# Method4: PhosphoRS + AutoRT + pDeep3
Method4 = data[data['autort_pdeep_Sequence_Label'] == True].copy()
Method4 = Method4.sort_values(['autort_pDeep_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

# Method5: PhosphoRS + PhosSight
Method5 = data[data['phosSight_Sequence_Label'] == True].copy()
Method5 = Method5.sort_values(['PhosSightProb', 'PhosphoRS_IsoformScore'], ascending=False)

# Method6: PhosphoRS + AutoRT + pDeep3 + PhosSight
Method6 = data[data['autort_pdeep_phosSight_Sequence_Label'] == True].copy()
Method6 = Method6.sort_values(['autort_pDeep_phosSight_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

# Method7: PhosphoRS + Rescore (without deep learning features)
rescore_percolator_results = rescore_percolator_results[rescore_percolator_results['q-value'] <= 0.01]
Method7 = deeprescore2_results[deeprescore2_results['Title'].isin(rescore_percolator_results['PSMId'])].copy()
Method7 = Method7[Method7['phosphors_Sequence_Label'] == True]
Method7 = Method7.sort_values(['PhosphoRS_IsoformProbability', 'PhosphoRS_IsoformScore'], ascending=False)

# Method8: PhosphoRS + DeepRescore (with deep learning features)
deeprescore2_percolator_results = deeprescore2_percolator_results[deeprescore2_percolator_results['q-value'] <= 0.01]
Method8 = deeprescore2_results[deeprescore2_results['Title'].isin(deeprescore2_percolator_results['PSMId'])].copy()
Method8 = Method8[Method8['phosphors_Sequence_Label'] == True]
Method8 = Method8.sort_values(['PhosphoRS_IsoformProbability', 'PhosphoRS_IsoformScore'], ascending=False)

# Method9: Method4 + DeepRescore
Method9 = deeprescore2_results[deeprescore2_results['Title'].isin(deeprescore2_percolator_results['PSMId'])].copy()
Method9 = Method9[Method9['autort_pdeep_Sequence_Label'] == True]
Method9 = Method9.sort_values(['autort_pDeep_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

# Method10: Method6 + DeepRescore
Method10 = deeprescore2_with_phossight_results[deeprescore2_with_phossight_results['Title'].isin(deeprescore2_percolator_results['PSMId'])].copy()
Method10 = Method10[Method10['autort_pdeep_phosSight_Sequence_Label'] == True]
Method10 = Method10.sort_values(['autort_pDeep_phosSight_Prob', 'PhosphoRS_IsoformScore'], ascending=False)

def calculate_flr(method_data, sequence_label_col, site_label_col):
    """Calculate FLR for a method using vectorized operations"""
    print(f"  Calculating FLR for {len(method_data)} rows...")
    
    # Vectorized operations for better performance
    seq_true = (method_data[sequence_label_col] == True) | (method_data[sequence_label_col] == 'TRUE')
    site_true = (method_data[site_label_col] == 'TRUE')
    
    # Both sequence and site labels must be true for correct identification
    correct = seq_true & site_true
    
    # Calculate cumulative sums
    TL_cumsum = correct.cumsum()
    FL_cumsum = (~correct).cumsum()
    
    # Calculate FLR
    FLR = FL_cumsum / (FL_cumsum + TL_cumsum)
    
    return TL_cumsum.tolist(), FL_cumsum.tolist(), FLR.tolist()

# Calculate FLR for each method (renamed and filtered)
methods = [
    (Method1, 'phosphors_Sequence_Label', 'phosphors_Site_Label', 'Method1'),
    (Method9, 'autort_pdeep_Sequence_Label', 'autort_pdeep_Site_Label', 'Method2'),  # Method9 -> Method2
    (Method7, 'phosphors_Sequence_Label', 'phosphors_Site_Label', 'Method3'),  # Method7 -> Method3
    (Method5, 'phosSight_Sequence_Label', 'phosSight_Site_Label', 'Method4'),  # Method5 -> Method4
    (Method10, 'autort_pdeep_phosSight_Sequence_Label', 'autort_pdeep_phosSight_Site_Label', 'Method5')  # Method10 -> Method5
]

print("Starting FLR calculations...")
print(f"Data shape after filtering: {data.shape}")
print(f"Method1 rows: {len(Method1)}")
print(f"Method4 rows: {len(Method4)}")
print(f"Method5 rows: {len(Method5)}")
print(f"Method6 rows: {len(Method6)}")

# Special case for Method7 - only use site_label
def calculate_flr_method7(method_data, site_label_col):
    """Calculate FLR for Method7 (only uses site_label)"""
    tl = 0
    fl = 0
    TL_list = []
    FL_list = []
    
    for i in range(len(method_data)):
        if method_data.iloc[i][site_label_col] == 'TRUE':
            tl += 1
        else:
            fl += 1
        TL_list.append(tl)
        FL_list.append(fl)
    
    FLR = np.array(FL_list) / (np.array(FL_list) + np.array(TL_list))
    return TL_list, FL_list, FLR

# Calculate FLR for all methods
plot_data = []
print("Starting FLR calculations for all methods...")

for i, (method_data, seq_col, site_col, method_name) in enumerate(methods):
    print(f'{method_name}: {len(method_data)}')
    
    if len(method_data) == 0:
        print(f'  Skipping {method_name} - no data')
        continue
    
    print(f'  Processing {method_name}...')
    if method_name == 'Method7':
        TL, FL, FLR = calculate_flr_method7(method_data, site_col)
    else:
        TL, FL, FLR = calculate_flr(method_data, seq_col, site_col)
    
    print(f'  FLR calculation completed for {method_name}')
    
    # Find count at FLR <= 0.01
    count = max([tl for tl, flr in zip(TL, FLR) if flr <= 0.01], default=0)
    print(f'count_{method_name}: {count}')
    
    # Add to plot data
    print(f'  Adding {len(TL)} data points to plot...')
    for tl, flr in zip(TL, FLR):
        plot_data.append({
            'TL': tl,
            'FLR': flr * 100,  # Convert to percentage
            'Method': method_name
        })
    
    print(f'  Completed {method_name}')

print("All FLR calculations completed!")

# Create DataFrame for plotting
plot_df = pd.DataFrame(plot_data)

# Check if we have any data to plot
if plot_df.empty:
    print("No data to plot - all methods returned 0 rows")
    exit()

print(f"Plot data shape: {plot_df.shape}")
print(f"Methods with data: {plot_df['Method'].unique()}")

# Create the plot
plt.figure(figsize=(11/2.54, 10/2.54))  # Convert cm to inches

# Define method order and specific colors (renamed methods)
method_order = ['Method1', 'Method2', 'Method3', 'Method4', 'Method5']
# Define specific colors for each method
method_colors = {
    'Method1': '#1f77b4',    # Blue
    'Method2': '#d62728',    # Red (was Method9)
    'Method3': '#8c564b',    # Brown (was Method7)
    'Method4': '#2ca02c',    # Green (was Method5)
    'Method5': '#17becf'     # Cyan (was Method10)
}

for method in method_order:
    method_data = plot_df[plot_df['Method'] == method]
    if not method_data.empty:
        plt.plot(method_data['TL'], method_data['FLR'], 
                label=method, color=method_colors[method], linewidth=1.5)

# Customize plot
plt.xlabel('#Correct PSMs')
plt.ylabel('FLR(%)')
plt.ylim(0, 6)
plt.xlim(6000, plot_df['TL'].max())
plt.yticks(range(7), range(7))

# Add horizontal line at FLR = 1%
plt.axhline(y=1, color='grey', linestyle='--', alpha=0.7)

# Customize appearance
plt.grid(False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(True)
plt.gca().spines['bottom'].set_visible(True)

# Add arrow to axes
# ax = plt.gca()
# ax.annotate('', xy=(1, 0), xytext=(0, 0), 
#            arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))
# ax.annotate('', xy=(0, 1), xytext=(0, 0), 
#            arrowprops=dict(arrowstyle='->', lw=1.5, color='black'))

# Rotate x-axis labels
plt.xticks(rotation=45, ha='right')

# Add legend
plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()

# Save plot
output_path = 'Figure3B_v2.svg'
plt.savefig(output_path, format='svg', bbox_inches='tight', dpi=300)
print(f'Plot saved to: {output_path}')

plt.show()
