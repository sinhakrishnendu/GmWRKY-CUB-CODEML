#!/bin/bash

# Loop through each clade folder (e.g., C1, C2, etc.)
for clade in C*/; do
    # Extract the clade name (removes the trailing '/')
    clade_name=$(basename "$clade")

    # Navigate to the clade folder
    cd "$clade" || exit

    # Run the lnL_np_extractor.sh script for the current clade
    bash ../lnL_np_extractor.sh

    # Move the generated CSV to the parent folder with clade name prefix
    mv lnL_np_values.csv "../${clade_name}_lrt_np_values.csv"

    # Navigate back to the parent directory
    cd ..

    # Run the Python analysis on the clade's CSV file
    #python3 lrt_test.py "${clade_name}_lrt_np_values.csv"
    
    python3 lrt_test_BHCorrection.py

done
