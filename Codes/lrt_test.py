import pandas as pd
import numpy as np
from scipy.stats import chi2, false_discovery_control
import os
import glob

# Function to perform LRT analysis on each file
def perform_lrt_analysis(csv_file):
    clade_name = os.path.splitext(os.path.basename(csv_file))[0]  # Extract clade name

    # Load the data from the current CSV file
    df = pd.read_csv(csv_file)

    # Get the unique gene names by stripping suffixes like _BS, _BS_NULL, and _B
    gene_names = df['Folder'].str.extract(r'(Glyma\.\d+G\d+\.\d+)_')[0].unique()

    # Initialize lists for storing results for branchsite model and branch model
    branchsite_results = []
    branch_model_results = []

    # LRT for branchsite model (_BS vs _BS_NULL)
    for gene in gene_names:
        bs_data = df[df['Folder'].str.contains(f'{gene}_BS')]
        null_data = df[df['Folder'].str.contains(f'{gene}_BS_NULL')]

        if not bs_data.empty and not null_data.empty:
            bs_data = bs_data.iloc[0]
            null_data = null_data.iloc[0]

            # Extract lnL and np values
            lnL_bs = bs_data['lnL']
            lnL_null = null_data['lnL']
            np_bs = bs_data['np']
            np_null = null_data['np']

            # Calculate LRT value and degrees of freedom
            lrt_value = 2 * (lnL_bs - lnL_null)
            df_branchsite = np_bs - np_null  # Degrees of freedom is the difference in np

            # Calculate the p-value for the LRT statistic
            p_value_branchsite = 1 - chi2.cdf(lrt_value, df_branchsite)

            result = {
                'Gene': gene,
                'lnL_BS': lnL_bs,
                'lnL_NULL': lnL_null,
                'LRT': lrt_value,
                'Degrees_of_Freedom': df_branchsite,
                'p_value': p_value_branchsite
            }

            branchsite_results.append(result)

    # LRT for branch model (_B vs M0)
    for gene in gene_names:
        branch_data = df[df['Folder'].str.contains(f'{gene}_B')]
        m0_data = df[df['Folder'].str.contains('M0')]

        if not branch_data.empty and not m0_data.empty:
            branch_data = branch_data.iloc[0]
            m0_data = m0_data.iloc[0]

            # Extract lnL and np values
            lnL_branch = branch_data['lnL']
            lnL_m0 = m0_data['lnL']
            np_branch = branch_data['np']
            np_m0 = m0_data['np']

            # Calculate LRT value and degrees of freedom
            lrt_value_branch = 2 * (lnL_branch - lnL_m0)
            df_branch_model = np_branch - np_m0  # Degrees of freedom is the difference in np

            # Calculate the p-value for the LRT statistic
            p_value_branch_model = 1 - chi2.cdf(lrt_value_branch, df_branch_model)

            result_branch = {
                'Gene': gene,
                'lnL_B': lnL_branch,
                'lnL_M0': lnL_m0,
                'LRT': lrt_value_branch,
                'Degrees_of_Freedom': df_branch_model,
                'p_value': p_value_branch_model
            }

            branch_model_results.append(result_branch)

    # Convert the results for the current file to DataFrames
    branchsite_results_df = pd.DataFrame(branchsite_results)
    branch_model_results_df = pd.DataFrame(branch_model_results)

    # Adjust p-values using Benjamini-Hochberg correction for both branchsite and branch model
    if not branchsite_results_df.empty:
        branchsite_results_df['adjusted_p_value'] = false_discovery_control(branchsite_results_df['p_value'].values, method='bh')

    if not branch_model_results_df.empty:
        branch_model_results_df['adjusted_p_value'] = false_discovery_control(branch_model_results_df['p_value'].values, method='bh')

    clade = clade_name.split('_')[0]

    # Write results to individual Excel file
    individual_output_file = f'{clade}_lrt_results.xlsx'
    with pd.ExcelWriter(individual_output_file, engine='openpyxl') as writer:
        branchsite_results_df.to_excel(writer, sheet_name='LRT for branchsite model', index=False)
        branch_model_results_df.to_excel(writer, sheet_name='LRT for branch model', index=False)

    print(f"Results saved to '{individual_output_file}'.")

# Get all CSV files in the current directory
csv_files = glob.glob('*.csv')

# Perform LRT analysis for each CSV file
for csv_file in csv_files:
    print(f"Processing file: {csv_file}")
    perform_lrt_analysis(csv_file)

print("LRT analysis completed for all files.")
