#this .sh file will run two-ration branch model and branch-site model from the base codeml.ctl file with necessary modification needed
#this file require one base codeml.ctl file and all the separate treefiles with the foreground branches and the MSA file

#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define the base control file
base_ctl_file="codeml.ctl"

# Check if the control file exists
if [ ! -f "$base_ctl_file" ]; then
    echo "Error: Base control file '$base_ctl_file' not found."
    exit 1
fi

# Loop through all .treefile files in the current directory
for treefile in *.treefile; do
    echo "Processing tree file: $treefile"

    # Extract base name without extension
    base_name="${treefile%.treefile}"

    # Extract seqfile from the base control file
    seqfile=$(grep -E '^ *seqfile *= *' "$base_ctl_file" | sed -E 's/.*= *([^ ]*).*/\1/')
    seqfile="${seqfile//\"/}"  # Trim possible quotes around the seqfile path

    # Check if the sequence file exists
    if [ ! -f "$seqfile" ]; then
        echo "Error: Sequence file '$seqfile' not found. Skipping $treefile."
        continue
    fi

    # ----- Branch Model -----
    # Create a new folder with prefix _B for branch model
    branch_folder="${base_name}_B"
    mkdir -p "$branch_folder"

    # Create a modified .ctl file for branch model
    branch_ctl_file="$branch_folder/$base_name.B.ctl"
    cp "$base_ctl_file" "$branch_ctl_file"
    
    # Modify the branch model parameters (model = 2, NSsites = 0)
    sed -i "s|^ *treefile *=.*|treefile = $treefile|" "$branch_ctl_file"
    sed -i "s|^ *model *=.*|model = 2|" "$branch_ctl_file"
    sed -i "s|^ *NSsites *=.*|NSsites = 0|" "$branch_ctl_file"

    # Copy necessary files into the branch model folder
    cp "$seqfile" "$treefile" "$branch_folder/"

    # Run codeml with the branch model control file
    (
        cd "$branch_folder"
        echo "  Running branch model codeml in $branch_folder/"
        codeml "$(basename "$branch_ctl_file")"
        echo "  Branch model codeml run completed for $branch_ctl_file."
    )

    # ----- Branch-Site Model -----
    # Create a new folder with prefix _BS for branch-site model
    bs_folder="${base_name}_BS"
    mkdir -p "$bs_folder"

    # Create a modified .ctl file for branch-site model
    bs_ctl_file="$bs_folder/$base_name.BS.ctl"
    cp "$base_ctl_file" "$bs_ctl_file"
    
    # Modify the branch-site model parameters (model = 2, NSsites = 2)
    sed -i "s|^ *treefile *=.*|treefile = $treefile|" "$bs_ctl_file"
    sed -i "s|^ *model *=.*|model = 2|" "$bs_ctl_file"
    sed -i "s|^ *NSsites *=.*|NSsites = 2|" "$bs_ctl_file"

    # Copy necessary files into the branch-site model folder
    cp "$seqfile" "$treefile" "$bs_folder/"

    # Run codeml with the branch-site model control file
    (
        cd "$bs_folder"
        echo "  Running branch-site model codeml in $bs_folder/"
        codeml "$(basename "$bs_ctl_file")"
        echo "  Branch-site model codeml run completed for $bs_ctl_file."
    )

    # ----- Null Model (for testing against the branch-site model) -----
    # Create a new folder with suffix _BS_NULL for null model
    null_folder="${base_name}_BS_NULL"
    mkdir -p "$null_folder"

    # Create a modified .ctl file for null model (fix_omega = 1, omega = 1)
    null_ctl_file="$null_folder/$base_name.BS_NULL.ctl"
    cp "$base_ctl_file" "$null_ctl_file"

    # Modify the null model parameters (model = 2, NSsites = 2, fix_omega = 1, omega = 1)
    sed -i "s|^ *treefile *=.*|treefile = $treefile|" "$null_ctl_file"
    sed -i "s|^ *model *=.*|model = 2|" "$null_ctl_file"
    sed -i "s|^ *NSsites *=.*|NSsites = 2|" "$null_ctl_file"
    sed -i "s|^ *fix_omega *=.*|fix_omega = 1|" "$null_ctl_file"
    sed -i "s|^ *omega *=.*|omega = 1|" "$null_ctl_file"

    # Copy necessary files into the null model folder
    cp "$seqfile" "$treefile" "$null_folder/"

    # Run codeml with the null model control file
    (
        cd "$null_folder"
        echo "  Running null model codeml in $null_folder/"
        codeml "$(basename "$null_ctl_file")"
        echo "  Null model codeml run completed for $null_ctl_file."
    )

    echo "Finished processing $treefile."
    echo "---------------------------------------"
done

echo "All processing completed."
