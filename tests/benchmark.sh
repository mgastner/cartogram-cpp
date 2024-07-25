#!/bin/bash

# Path to the file containing the commands
test_folder=$1

# Check if the commands file exists
if [[ ! -d "./$test_tolder" ]]; then
    echo "Directoy not found: $test_folder"
    exit 1
fi

cd "$test_folder"

# Read each command from the file and execute it
while IFS= read -r cmd; do
    # Run hyperfine and export results to a temporary CSV file
    hyperfine "$cmd" --export-csv tmp.csv
    # Concatenate results while handling headers
    if [ ! -f final_results.csv ]; then
        # If the final results file does not exist, create it with the header
        cp tmp.csv final_results.csv
    else
        # Otherwise, append the new results, excluding the header
        tail -n +2 tmp.csv >> final_results.csv
    fi
done < "successful_runs.txt"

# Clean up the temporary file
rm tmp.csv

# Print full file location
echo "Results saved to: $(pwd)/final_results.csv"