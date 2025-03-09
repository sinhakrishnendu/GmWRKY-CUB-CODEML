#this script scroll throuh all the folders in the director it stays and search for output.txt file. 
#then it search lnL and np values and store them in a .csv file

#!/bin/bash

# Output CSV file
output_file="lnL_np_values.csv"

# Initialize the CSV file with headers
echo "Folder,lnL,np" > "$output_file"

# Loop through all directories and subdirectories
for dir in */ ; do
    # Remove the trailing '/' from the folder name
    folder_name=$(basename "$dir")

    # Check if the output.txt file exists in the current folder
    if [[ -f "$dir/output.txt" ]]; then
        # Extract the lnL and np values from the output.txt file
        lnL_np_line=$(grep "lnL(ntime" "$dir/output.txt")

        if [[ ! -z "$lnL_np_line" ]]; then
            # Use awk to extract lnL and np values
            lnL=$(echo "$lnL_np_line" | awk '{print $5}')
            np=$(echo "$lnL_np_line" | awk -F "np:|\\)" '{print $2}' | awk '{print $1}')

            # Append the results to the CSV file
            echo "$folder_name,$lnL,$np" >> "$output_file"
        fi
    fi
done

echo "Extraction complete. Results saved to $output_file."
