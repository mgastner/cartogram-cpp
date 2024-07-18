#!/bin/bash

# Path to the file containing the commands
commands_file="commands.txt"

# Check if the commands file exists
if [[ ! -f "$commands_file" ]]; then
    echo "File not found: $commands_file"
    exit 1
fi

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
done < "$commands_file"

# Clean up the temporary file
rm tmp.csv
