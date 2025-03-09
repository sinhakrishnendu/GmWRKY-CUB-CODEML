import pandas as pd
import numpy as np
from scipy.stats import chi2
import openpyxl
import os

# Initialize an empty list to store the results for branchsite model and branch model
branchsite_results = []
branch_model_results = []

# Critical chi-square values for 1 degree of freedom at various significance levels
significance_levels = [0.05, 0.001, 0.0003]
critical_values = {alpha: chi2.ppf(1 - alpha, df=1) for alpha in significance_levels}

# Get the list of all CSV files in the current directory
csv_files = [file for file in os.listdir() if file.endswith('.csv')]

# Loop through each CSV file
for csv_file in csv_files:
    print(f"Processing file: {csv_file}")
    
    # Load the data from the current CSV file
    df = pd.read_csv(csv_file)

    # Get the unique gene names by stripping suffixes like _BS, _BS_NULL, and _B
    gene_names = df['Folder'].str.extract(r'(Glyma\.\d+G\d+\.\d+)_')[0].unique()

    # LRT for branchsite model (_BS vs _BS_NULL)
    for gene in gene_names:
        # Extract the BS and BS_NULL data for the gene
        bs_data = df[df['Folder'].str.contains(f'{gene}_BS')]
        null_data = df[df['Folder'].str.contains(f'{gene}_BS_NULL')]

        # Check if both BS and BS_NULL entries exist
        if not bs_data.empty and not null_data.empty:
            # Extract the first matching row for BS and BS_NULL models
            bs_data = bs_data.iloc[0]
            null_data = null_data.iloc[0]

            # Extract lnL values
            lnL_bs = bs_data['lnL']
            lnL_null = null_data['lnL']

            # Calculate LRT value
            lrt_value = 2 * (lnL_null - lnL_bs)

            # Determine significance at each level
            significance_results = {alpha: "Yes" if lrt_value > crit_val else "No"
                                    for alpha, crit_val in critical_values.items()}

            # Append results to the list for branchsite model
            result = {
                'Gene': gene,
                'lnL_BS': lnL_bs,
                'lnL_NULL': lnL_null,
                'LRT': lrt_value
            }

            # Add critical values and significance results to the dictionary
            for alpha, crit_val in critical_values.items():
                result[f'Critical_Chi2_{alpha}'] = crit_val
                result[f'Significant_{alpha}'] = significance_results[alpha]

            branchsite_results.append(result)
        else:
            print(f"Missing BS or BS_NULL data for gene: {gene}")

    # LRT for branch model (_B vs M0)
    for gene in gene_names:
        # Extract the _B data for the gene and the M0 data for comparison
        branch_data = df[df['Folder'].str.contains(f'{gene}_B')]
        m0_data = df[df['Folder'].str.contains('M0')]

        # Check if both branch model (_B) and M0 entries exist
        if not branch_data.empty and not m0_data.empty:
            # Extract the first matching row for _B and M0 models
            branch_data = branch_data.iloc[0]
            m0_data = m0_data.iloc[0]

            # Extract lnL values
            lnL_branch = branch_data['lnL']
            lnL_m0 = m0_data['lnL']

            # Calculate LRT value
            lrt_value_branch = 2 * (lnL_branch - lnL_m0)

            # Determine significance at each level
            significance_results_branch = {alpha: "Yes" if lrt_value_branch > crit_val else "No"
                                           for alpha, crit_val in critical_values.items()}

            # Append results to the list for branch model
            result_branch = {
                'Gene': gene,
                'lnL_B': lnL_branch,
                'lnL_M0': lnL_m0,
                'LRT': lrt_value_branch
            }

            # Add critical values and significance results to the dictionary
            for alpha, crit_val in critical_values.items():
                result_branch[f'Critical_Chi2_{alpha}'] = crit_val
                result_branch[f'Significant_{alpha}'] = significance_results_branch[alpha]

            branch_model_results.append(result_branch)
        else:
            print(f"Missing _B or M0 data for gene: {gene}")

# Convert the results to DataFrames
branchsite_results_df = pd.DataFrame(branchsite_results)
branch_model_results_df = pd.DataFrame(branch_model_results)

# Write both sheets to the same Excel file
output_file = 'lrt_results_multiple_significance.xlsx'
with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    # Write branchsite model results
    branchsite_results_df.to_excel(writer, sheet_name='LRT for branchsite model', index=False)

    # Write branch model results
    branch_model_results_df.to_excel(writer, sheet_name='LRT for branch model', index=False)

print(f"LRT tests completed. Results saved to '{output_file}' with two sheets.")
