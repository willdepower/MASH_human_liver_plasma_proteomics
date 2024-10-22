```python
import pandas as pd
import numpy as np
import pingouin as pg
import statsmodels.stats.multitest as multi

# Step 1: Load the MASH_liver_secreted and patients (clinical) data
liver_data = pd.read_excel('')  # MASH_liver_secreted proteins data set
patients_data = pd.read_excel('')  # Clinical data

# Step 2: Set 'ID' as the index for easier merging
liver_data.set_index('ID', inplace=True)
patients_data.set_index('ID', inplace=True)

# Step 3: Impute missing values in MASH_liver_secreted data
def imputation_normal_distribution(df):
    """
    Imputes missing values using a down-shifted Gaussian distribution.
    """
    data_imputed = df.copy()
    for col in data_imputed.columns:
        if data_imputed[col].isnull().any():
            missing = data_imputed[col].isnull()
            std = data_imputed[col].std()
            mean = data_imputed[col].mean()
            sigma = std * 0.3
            mu = mean - (std * 1.8)
            np.random.seed(123)
            data_imputed.loc[missing, col] = np.random.normal(mu, sigma, size=missing.sum())
    return data_imputed

# Perform imputation on the liver dataset
imputed_liver_data = imputation_normal_distribution(liver_data)



```


```python
# Step 4: Merge the imputed liver data with clinical data (patients) on the 'ID'
merged_liver = imputed_liver_data.join(patients_data, how='inner')

merged_liver
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>VWA8</th>
      <th>TRIP11</th>
      <th>LPCAT2</th>
      <th>CAPN2</th>
      <th>CD63</th>
      <th>IGHG4</th>
      <th>GOLGA5</th>
      <th>S100A8</th>
      <th>RAI14</th>
      <th>HBG1</th>
      <th>...</th>
      <th>TFIP11</th>
      <th>ABCF2</th>
      <th>CD84</th>
      <th>CRYL1.1</th>
      <th>SIGLEC9</th>
      <th>TMED5</th>
      <th>Age</th>
      <th>BMI</th>
      <th>nas_steatosis_ordinal</th>
      <th>Gender</th>
    </tr>
    <tr>
      <th>ID</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>302</th>
      <td>10.868683</td>
      <td>15.837627</td>
      <td>18.302187</td>
      <td>27.421979</td>
      <td>14.638799</td>
      <td>19.249382</td>
      <td>14.903990</td>
      <td>31.105998</td>
      <td>16.939341</td>
      <td>22.392148</td>
      <td>...</td>
      <td>13.015443</td>
      <td>NaN</td>
      <td>13.544364</td>
      <td>NaN</td>
      <td>10.834764</td>
      <td>14.658653</td>
      <td>58</td>
      <td>39.453125</td>
      <td>2</td>
      <td>0</td>
    </tr>
    <tr>
      <th>303</th>
      <td>11.814315</td>
      <td>10.098213</td>
      <td>12.091767</td>
      <td>23.371925</td>
      <td>15.260578</td>
      <td>15.337399</td>
      <td>15.157019</td>
      <td>20.302369</td>
      <td>15.383156</td>
      <td>19.056913</td>
      <td>...</td>
      <td>14.437900</td>
      <td>NaN</td>
      <td>13.984611</td>
      <td>NaN</td>
      <td>11.410622</td>
      <td>15.261853</td>
      <td>41</td>
      <td>58.333333</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>305</th>
      <td>11.490005</td>
      <td>14.271902</td>
      <td>12.900737</td>
      <td>26.396626</td>
      <td>15.047336</td>
      <td>13.703069</td>
      <td>11.524556</td>
      <td>22.426601</td>
      <td>16.913282</td>
      <td>21.698084</td>
      <td>...</td>
      <td>13.369850</td>
      <td>NaN</td>
      <td>13.833626</td>
      <td>NaN</td>
      <td>11.213128</td>
      <td>14.134934</td>
      <td>51</td>
      <td>55.574296</td>
      <td>3</td>
      <td>0</td>
    </tr>
    <tr>
      <th>306</th>
      <td>16.453802</td>
      <td>11.203775</td>
      <td>13.331585</td>
      <td>20.710538</td>
      <td>14.513229</td>
      <td>14.560178</td>
      <td>12.341919</td>
      <td>20.416828</td>
      <td>11.687406</td>
      <td>18.957032</td>
      <td>...</td>
      <td>13.248305</td>
      <td>NaN</td>
      <td>13.455454</td>
      <td>NaN</td>
      <td>10.718467</td>
      <td>15.054983</td>
      <td>42</td>
      <td>48.069220</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>307</th>
      <td>14.352145</td>
      <td>10.824617</td>
      <td>12.623296</td>
      <td>19.788876</td>
      <td>14.790150</td>
      <td>16.207486</td>
      <td>12.061600</td>
      <td>20.426717</td>
      <td>12.605988</td>
      <td>19.651508</td>
      <td>...</td>
      <td>12.943869</td>
      <td>NaN</td>
      <td>13.651527</td>
      <td>NaN</td>
      <td>10.974937</td>
      <td>14.536835</td>
      <td>25</td>
      <td>52.892562</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>411</th>
      <td>11.574201</td>
      <td>13.113536</td>
      <td>11.998311</td>
      <td>21.841983</td>
      <td>18.006723</td>
      <td>17.508900</td>
      <td>14.442554</td>
      <td>17.589383</td>
      <td>14.168363</td>
      <td>19.437170</td>
      <td>...</td>
      <td>13.337517</td>
      <td>NaN</td>
      <td>14.073371</td>
      <td>NaN</td>
      <td>12.856697</td>
      <td>15.309594</td>
      <td>40</td>
      <td>44.444444</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>412</th>
      <td>10.984210</td>
      <td>11.096305</td>
      <td>13.063835</td>
      <td>22.241998</td>
      <td>17.764147</td>
      <td>14.589825</td>
      <td>13.596222</td>
      <td>17.923675</td>
      <td>12.536555</td>
      <td>18.395520</td>
      <td>...</td>
      <td>14.810620</td>
      <td>NaN</td>
      <td>13.944431</td>
      <td>NaN</td>
      <td>12.936167</td>
      <td>14.655302</td>
      <td>44</td>
      <td>46.708570</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>414</th>
      <td>11.889157</td>
      <td>10.841227</td>
      <td>12.826904</td>
      <td>21.041622</td>
      <td>16.755738</td>
      <td>13.831831</td>
      <td>14.903090</td>
      <td>18.478120</td>
      <td>14.730740</td>
      <td>17.598572</td>
      <td>...</td>
      <td>13.146603</td>
      <td>NaN</td>
      <td>13.783432</td>
      <td>NaN</td>
      <td>13.605878</td>
      <td>14.358217</td>
      <td>22</td>
      <td>43.983673</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>415</th>
      <td>10.863429</td>
      <td>9.970494</td>
      <td>12.531063</td>
      <td>21.181317</td>
      <td>15.982149</td>
      <td>14.476860</td>
      <td>14.920917</td>
      <td>19.225106</td>
      <td>14.628562</td>
      <td>18.089012</td>
      <td>...</td>
      <td>13.015397</td>
      <td>NaN</td>
      <td>15.642842</td>
      <td>NaN</td>
      <td>10.921274</td>
      <td>15.274126</td>
      <td>28</td>
      <td>37.920165</td>
      <td>2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>416</th>
      <td>10.397691</td>
      <td>11.426670</td>
      <td>12.422869</td>
      <td>20.746491</td>
      <td>16.029508</td>
      <td>14.279105</td>
      <td>11.895596</td>
      <td>17.798464</td>
      <td>11.536122</td>
      <td>17.383069</td>
      <td>...</td>
      <td>13.075533</td>
      <td>NaN</td>
      <td>15.312404</td>
      <td>NaN</td>
      <td>11.354644</td>
      <td>14.856227</td>
      <td>48</td>
      <td>37.537112</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>96 rows × 3337 columns</p>
</div>




```python
# Step 5: Select relevant clinical variables
clinical_vars = ['Age', 'BMI', 'Gender']  # Covariates
y_col = 'nas_steatosis_ordinal'  # Steatosis stages (ordinal)

# Ensure numeric conversion of clinical variables
merged_liver[clinical_vars + [y_col]] = merged_liver[clinical_vars + [y_col]].apply(pd.to_numeric, errors='coerce')

# Step 6: Define a function to perform partial correlation and FDR correction
def partial_corr_with_fdr(data, proteins, covariates, y_col):
    results = []
    for protein in proteins:
        # Drop rows where the protein or target variable is missing
        subset_data = data.dropna(subset=[protein, y_col])
        
        # Only run partial correlation if enough samples are available (at least 3)
        if subset_data.shape[0] > 2:
            result = pg.partial_corr(data=subset_data,
                                     x=protein,
                                     y=y_col,
                                     covar=covariates,
                                     method='spearman')  # Spearman correlation for ordinal data
            result['protein'] = protein
            results.append(result)
    
    # Convert to DataFrame
    results_df = pd.concat(results, ignore_index=True)
    
    # FDR correction for multiple hypothesis testing
    reject, pvals_corrected = multi.fdrcorrection(results_df['p-val'], alpha=0.05)
    results_df['pval_corrected'] = pvals_corrected
    results_df['significant'] = reject
    
    return results_df
```


```python

# Step 7: Get the list of liver proteins (exclude the clinical variables)
liver_proteins = merged_liver.columns.difference(clinical_vars + [y_col])

# Step 8: Run the partial correlation analysis
liver_corr_results = partial_corr_with_fdr(data=merged_liver, 
                                           proteins=liver_proteins, 
                                           covariates=clinical_vars, 
                                           y_col=y_col)

# Step 9: Filter and display significant correlations
liver_significant = liver_corr_results[liver_corr_results['significant']]

# Display results
print("Significant liver-secreted proteins correlated with steatosis stages:")
print(liver_significant)
```

    Significant liver-secreted proteins correlated with steatosis stages:
           n         r           CI95%     p-val  protein  pval_corrected  \
    9     96 -0.335475   [-0.5, -0.14]  0.001012     AASS        0.049256   
    120   96  0.397629    [0.21, 0.56]  0.000079      AFM        0.020237   
    129   96  0.409629    [0.22, 0.57]  0.000046     AGRN        0.012984   
    167   96  0.363050    [0.17, 0.53]  0.000348      ALB        0.035028   
    181   96 -0.339549  [-0.51, -0.15]  0.000869  ALDH6A1        0.048946   
    ...   ..       ...             ...       ...      ...             ...   
    3137  96 -0.356893  [-0.52, -0.17]  0.000445   UBE2D3        0.036412   
    3162  96 -0.357256  [-0.52, -0.17]  0.000439    UCHL5        0.036412   
    3164  96  0.336713    [0.14, 0.51]  0.000966     UFD1        0.048946   
    3184  96 -0.349098  [-0.52, -0.16]  0.000604     UMPS        0.042258   
    3241  96  0.356827    [0.17, 0.52]  0.000447   VPS13A        0.036412   
    
          significant  
    9            True  
    120          True  
    129          True  
    167          True  
    181          True  
    ...           ...  
    3137         True  
    3162         True  
    3164         True  
    3184         True  
    3241         True  
    
    [69 rows x 7 columns]
    


```python
# Export the significant correlations to an Excel spreadsheet
output_file = 'liver_significant_proteins_correlated_with_steatosis_partcorrel_gender_binary.xlsx'
liver_significant.to_excel(output_file, index=False)

print(f"Results successfully exported to {output_file}")

```

    Results successfully exported to liver_significant_proteins_correlated_with_steatosis_partcorrel_gender_binary.xlsx
    
