import pandas as pd
from scipy.stats import chi2
import os

# Get the current working directory (where the script is located)
directory = os.getcwd()  # This gets the path of the current working directory

for filename in os.listdir(directory):
    # Check if the file is a CSV file
    if filename.endswith('.csv'):
        # Construct the full file path
        file_path = os.path.join(directory, filename)
        
        # Load the CSV file into a DataFrame
        df = pd.read_csv(file_path)# Load the data from CSV

        # Initialize an empty list to store the results for branchsite model and branch model
        branchsite_results = []
        branch_model_results = []

        # Critical chi-square values for 1 degree of freedom at various significance levels
        significance_levels = [0.05, 0.001, 0.0003]
        critical_values = {alpha: chi2.ppf(1 - alpha, df=1) for alpha in significance_levels}

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
                lrt_value = 2 * (lnL_bs - lnL_null)
                
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

        # Extract the file name from the file path
        file_name = os.path.basename(file_path)

        # Split the file name by underscores and get the first part (clade name)
        clade_name = file_name.split('_')[0]

        # Write both sheets to the same Excel file
        with pd.ExcelWriter(f'{clade_name}_lrt_results_multiple_significance.xlsx', engine='openpyxl') as writer:
            # Write branchsite model results
            branchsite_results_df.to_excel(writer, sheet_name='LRT for branchsite model', index=False)
            
            # Write branch model results
            branch_model_results_df.to_excel(writer, sheet_name='LRT for branch model', index=False)

        print(f"LRT tests completed. Results saved to '{clade_name}lrt_results_multiple_significance.xlsx' with two sheets.")