```python
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

# Load the two spreadsheets
# Replace 'file1.xlsx' and 'file2.xlsx' with the actual filenames
spreadsheet1 = pd.read_excel('MASH_human_liver_secreted.xlsx')
spreadsheet2 = pd.read_excel('WD_LS_perSteatosisliver.xlsx')

# Ensure the 'ID' column is in both datasets
spreadsheet1 = spreadsheet1.set_index('ID')
spreadsheet2 = spreadsheet2.set_index('ID')

# Select columns from C (index 2) to DXG (index 3002) in spreadsheet1
columns_to_analyze = spreadsheet1.iloc[:, 2:3002]

# Ensure 'Steatosis' is present in spreadsheet2
if 'Steatosis' not in spreadsheet2.columns:
    raise ValueError("Steatosis column not found in spreadsheet2")

# Merge the two datasets on 'ID', keeping only the relevant columns
merged_data = pd.merge(columns_to_analyze, spreadsheet2[['Steatosis']], on='ID')

# Extract the Steatosis values
Steatosis_values = merged_data['Steatosis'].values

# Initialize lists to store results
correlation_results = []
p_values = []

# Perform Spearman's correlation for Steatosis with each column in spreadsheet1
for column in columns_to_analyze.columns:
    variable_values = merged_data[column].values

    # Remove rows with NaN or Inf values in either the column or Steatosis
    mask = ~np.isnan(variable_values) & ~np.isnan(Steatosis_values) & np.isfinite(variable_values) & np.isfinite(Steatosis_values)
    valid_Steatosis = Steatosis_values[mask]
    valid_variable_values = variable_values[mask]

    # Ensure a minimum of 30 valid data points before calculating correlation
    if len(valid_Steatosis) < 30 or len(valid_variable_values) < 30:
        print(f"Skipping correlation for {column}: less than 30 valid data points")
        continue

    # Perform Spearman's correlation
    correlation, p_value = spearmanr(valid_Steatosis, valid_variable_values)

    # Store the result for Steatosis and the variable
    correlation_results.append({
        'Variable': column,
        'Correlation': correlation,
        'P-value': p_value
    })
    p_values.append(p_value)

# Apply FDR correction (Benjamini-Hochberg method) to the p-values
if len(p_values) > 0:
    rejected, q_values, _, _ = multipletests(p_values, method='fdr_bh')

    # Add q-values and significance to correlation results
    for i, result in enumerate(correlation_results):
        result['Q-value'] = q_values[i]
        result['Significant'] = rejected[i]
else:
    print("No valid correlations were found. Skipping FDR correction.")

# Convert the results to a DataFrame
results_df = pd.DataFrame(correlation_results)

# Save the results to an Excel file if there are valid correlations
if not results_df.empty:
    results_df.to_excel('steatosis_correlation_results.xlsx', index=False)
    print("Results saved to 'steatosis_correlation_results.xlsx'")
else:
    print("No results to save.")

# Display the top results for quick review
if not results_df.empty:
    print(results_df.head())

```

    Skipping correlation for TRIP11: less than 30 valid data points
    Skipping correlation for LPCAT2: less than 30 valid data points
    Skipping correlation for HECTD1: less than 30 valid data points
    Skipping correlation for IGHD: less than 30 valid data points
    Skipping correlation for SCRN1: less than 30 valid data points
    Skipping correlation for CASP14: less than 30 valid data points
    Skipping correlation for ACTB: less than 30 valid data points
    Skipping correlation for ASB15: less than 30 valid data points
    Skipping correlation for SERPINB4: less than 30 valid data points
    Skipping correlation for APOC1: less than 30 valid data points
    Skipping correlation for CAMK2G: less than 30 valid data points
    Skipping correlation for IGLV1-47: less than 30 valid data points
    Skipping correlation for NF2: less than 30 valid data points
    Skipping correlation for SLC25A1: less than 30 valid data points
    Skipping correlation for FLG2: less than 30 valid data points
    Skipping correlation for C1RL: less than 30 valid data points
    Skipping correlation for OLFM4: less than 30 valid data points
    Skipping correlation for PTRH2: less than 30 valid data points
    Skipping correlation for GOLGA4: less than 30 valid data points
    Skipping correlation for SMARCE1: less than 30 valid data points
    Skipping correlation for EFEMP1: less than 30 valid data points
    Skipping correlation for LECT2: less than 30 valid data points
    Skipping correlation for HLA-A.1: less than 30 valid data points
    Skipping correlation for CXCL2: less than 30 valid data points
    Skipping correlation for CSAD: less than 30 valid data points
    Skipping correlation for GTF2A1: less than 30 valid data points
    Skipping correlation for AUP1: less than 30 valid data points
    Skipping correlation for GTPBP4: less than 30 valid data points
    Skipping correlation for POGLUT2: less than 30 valid data points
    Skipping correlation for Unnamed: 967: less than 30 valid data points
    Skipping correlation for HNRNPD: less than 30 valid data points
    Skipping correlation for HBG2: less than 30 valid data points
    Skipping correlation for VNN2: less than 30 valid data points
    Skipping correlation for MTCH2: less than 30 valid data points
    Skipping correlation for 2024-03-01 00:00:00: less than 30 valid data points
    Skipping correlation for HPN: less than 30 valid data points
    Skipping correlation for SPTB: less than 30 valid data points
    Skipping correlation for TAGLN3: less than 30 valid data points
    Skipping correlation for ITGB2: less than 30 valid data points
    Skipping correlation for PPP1R12C: less than 30 valid data points
    Skipping correlation for NPM1.1: less than 30 valid data points
    Skipping correlation for SMARCA5: less than 30 valid data points
    Skipping correlation for PPP1R21: less than 30 valid data points
    Skipping correlation for HLA-A.2: less than 30 valid data points
    Skipping correlation for PRPF4: less than 30 valid data points
    Skipping correlation for AQR: less than 30 valid data points
    Skipping correlation for CKM: less than 30 valid data points
    Skipping correlation for CD300A: less than 30 valid data points
    Skipping correlation for C17orf49;C17orf49;RNASEK-C17orf49;BAP18;BAP18;BAP18: less than 30 valid data points
    Skipping correlation for CP.1: less than 30 valid data points
    Skipping correlation for METAP1D: less than 30 valid data points
    Skipping correlation for SORBS2: less than 30 valid data points
    Skipping correlation for BLOC1S3: less than 30 valid data points
    Skipping correlation for POLB: less than 30 valid data points
    Skipping correlation for MMP2: less than 30 valid data points
    Skipping correlation for MTHFD1.1: less than 30 valid data points
    Skipping correlation for FAU: less than 30 valid data points
    Skipping correlation for SAA1: less than 30 valid data points
    Skipping correlation for MAPK1IP1L: less than 30 valid data points
    Skipping correlation for F5: less than 30 valid data points
    Skipping correlation for APOC4: less than 30 valid data points
    Skipping correlation for GIT2: less than 30 valid data points
    Skipping correlation for ABCD3: less than 30 valid data points
    Skipping correlation for CES2.1: less than 30 valid data points
    Skipping correlation for MYL1: less than 30 valid data points
    Skipping correlation for IL6ST: less than 30 valid data points
    Skipping correlation for SLC39A5: less than 30 valid data points
    Skipping correlation for PKP2: less than 30 valid data points
    Skipping correlation for CCDC9: less than 30 valid data points
    Skipping correlation for AMY2B;AMY2B;AMY2A;AMY1C;AMY2A;AMY2A;AMY2B;AMY2B;AMY1A;AMY1B: less than 30 valid data points
    Skipping correlation for TRIOBP: less than 30 valid data points
    Skipping correlation for NFS1: less than 30 valid data points
    Skipping correlation for BPI: less than 30 valid data points
    Skipping correlation for YLPM1: less than 30 valid data points
    Skipping correlation for NDUFAF7: less than 30 valid data points
    Skipping correlation for PTPN6: less than 30 valid data points
    Skipping correlation for EFNA1: less than 30 valid data points
    Skipping correlation for PLTP: less than 30 valid data points
    Skipping correlation for YME1L1: less than 30 valid data points
    Skipping correlation for MYL3: less than 30 valid data points
    Skipping correlation for ZCCHC8: less than 30 valid data points
    Skipping correlation for RAB8A: less than 30 valid data points
    Skipping correlation for SERPINA1.1: less than 30 valid data points
    Skipping correlation for MECR: less than 30 valid data points
    Skipping correlation for DBN1: less than 30 valid data points
    Skipping correlation for CDC42BPB: less than 30 valid data points
    Skipping correlation for NEU1: less than 30 valid data points
    Skipping correlation for PPP1R11: less than 30 valid data points
    Skipping correlation for MYH6: less than 30 valid data points
    Skipping correlation for TNNI3: less than 30 valid data points
    Skipping correlation for PCNP: less than 30 valid data points
    Skipping correlation for DBH: less than 30 valid data points
    Skipping correlation for VIPR1: less than 30 valid data points
    Skipping correlation for GRN: less than 30 valid data points
    Skipping correlation for PPA2.1: less than 30 valid data points
    Skipping correlation for COL4A1: less than 30 valid data points
    Skipping correlation for IGHV4-4;IGHV4-30-2;;IGHV4-61;IGHV4-39;IGHV4-59;IGHV4-34;IGHV4-30-4;IGHV4-31;IGHV4-38-2: less than 30 valid data points
    Skipping correlation for ADPRH: less than 30 valid data points
    Skipping correlation for COQ8A: less than 30 valid data points
    Skipping correlation for HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DPA1;HLA-DQA1;HLA-DPA1;HLA-DPA1;HLA-DPA1: less than 30 valid data points
    Skipping correlation for NUP50: less than 30 valid data points
    Skipping correlation for NOP10: less than 30 valid data points
    Skipping correlation for C12orf57: less than 30 valid data points
    Skipping correlation for MCC: less than 30 valid data points
    Skipping correlation for RETSAT: less than 30 valid data points
    Skipping correlation for SART1: less than 30 valid data points
    Skipping correlation for KNG1.1: less than 30 valid data points
    Skipping correlation for NDUFA6: less than 30 valid data points
    Skipping correlation for MOCS1: less than 30 valid data points
    Skipping correlation for GAPDH: less than 30 valid data points
    Skipping correlation for SYPL1: less than 30 valid data points
    Skipping correlation for TRIM56: less than 30 valid data points
    Skipping correlation for THOC1: less than 30 valid data points
    Skipping correlation for HLA-B.1: less than 30 valid data points
    Skipping correlation for IGHV5-51: less than 30 valid data points
    Skipping correlation for CDH13: less than 30 valid data points
    Skipping correlation for FKBP7: less than 30 valid data points
    Skipping correlation for ITSN2: less than 30 valid data points
    Skipping correlation for GSTM4: less than 30 valid data points
    Skipping correlation for MROH2B: less than 30 valid data points
    Skipping correlation for MMP10: less than 30 valid data points
    Skipping correlation for CEMIP2: less than 30 valid data points
    Skipping correlation for DHRS2: less than 30 valid data points
    Skipping correlation for HLA-DRB1: less than 30 valid data points
    Skipping correlation for RBM3: less than 30 valid data points
    Skipping correlation for NDUFS6: less than 30 valid data points
    Skipping correlation for MRRF: less than 30 valid data points
    Skipping correlation for LGALS7B: less than 30 valid data points
    Skipping correlation for TNNC1: less than 30 valid data points
    Skipping correlation for WDR48: less than 30 valid data points
    Skipping correlation for CHMP4B: less than 30 valid data points
    Skipping correlation for DMTN: less than 30 valid data points
    Skipping correlation for ALDH9A1.1: less than 30 valid data points
    Skipping correlation for PCOLCE: less than 30 valid data points
    Skipping correlation for MAPRE2: less than 30 valid data points
    Skipping correlation for OXLD1: less than 30 valid data points
    Skipping correlation for TIMP2: less than 30 valid data points
    Skipping correlation for CCDC188: less than 30 valid data points
    Skipping correlation for OAF: less than 30 valid data points
    Skipping correlation for BORCS7-ASMT;BORCS7: less than 30 valid data points
    Skipping correlation for RBM12B: less than 30 valid data points
    Skipping correlation for ZNF638: less than 30 valid data points
    Skipping correlation for NDUFA12: less than 30 valid data points
    Skipping correlation for TUSC1: less than 30 valid data points
    Skipping correlation for DHFR: less than 30 valid data points
    Skipping correlation for CMA1: less than 30 valid data points
    Skipping correlation for UQCRB: less than 30 valid data points
    Skipping correlation for NECTIN3: less than 30 valid data points
    Skipping correlation for HLA-B.3: less than 30 valid data points
    Skipping correlation for MLIP: less than 30 valid data points
    Skipping correlation for ALAD: less than 30 valid data points
    Skipping correlation for HRNR: less than 30 valid data points
    Skipping correlation for HTRA4: less than 30 valid data points
    Skipping correlation for LMNA.1: less than 30 valid data points
    Skipping correlation for KRT6A: less than 30 valid data points
    Skipping correlation for BCLAF1: less than 30 valid data points
    Skipping correlation for PAF1: less than 30 valid data points
    Skipping correlation for RBM7;RBM7;RBM7;;RBM7;RBM7: less than 30 valid data points
    Skipping correlation for CSF3: less than 30 valid data points
    Skipping correlation for DNAH8: less than 30 valid data points
    Skipping correlation for DDX41: less than 30 valid data points
    Skipping correlation for KRT83: less than 30 valid data points
    Skipping correlation for SMOC1: less than 30 valid data points
    Skipping correlation for STC1: less than 30 valid data points
    Skipping correlation for ZNF326: less than 30 valid data points
    Skipping correlation for PPIG: less than 30 valid data points
    Skipping correlation for CNN1: less than 30 valid data points
    Skipping correlation for CSDE1: less than 30 valid data points
    Skipping correlation for ITGAM: less than 30 valid data points
    Skipping correlation for HLA-DRB1.2: less than 30 valid data points
    Skipping correlation for DSG1: less than 30 valid data points
    Skipping correlation for SCAF4: less than 30 valid data points
    Skipping correlation for AFDN.1: less than 30 valid data points
    Skipping correlation for MYL2: less than 30 valid data points
    Skipping correlation for PDXDC1: less than 30 valid data points
    Skipping correlation for HLA-DRB1.3: less than 30 valid data points
    Skipping correlation for CPE: less than 30 valid data points
    Skipping correlation for SCPEP1: less than 30 valid data points
    Skipping correlation for PRXL2B: less than 30 valid data points
    Skipping correlation for CA3: less than 30 valid data points
    Skipping correlation for Unnamed: 2925: less than 30 valid data points
    Skipping correlation for ENDOD1: less than 30 valid data points
    Skipping correlation for NUP133: less than 30 valid data points
    Skipping correlation for SMPDL3A: less than 30 valid data points
    Skipping correlation for PTGFRN: less than 30 valid data points
    Skipping correlation for FBLN2: less than 30 valid data points
    Skipping correlation for PELP1: less than 30 valid data points
    Skipping correlation for GAK: less than 30 valid data points
    Skipping correlation for CARNMT1: less than 30 valid data points
    Results saved to 'steatosis_correlation_results.xlsx'
      Variable  Correlation   P-value   Q-value  Significant
    0    CAPN2     0.271161  0.012067  0.118332        False
    1     CD63     0.317843  0.033364  0.183175        False
    2    IGHG4     0.504581  0.003795  0.066265        False
    3   GOLGA5     0.407559  0.001496  0.039680         True
    4   S100A8     0.362319  0.000656  0.026070         True
    


```python

```
