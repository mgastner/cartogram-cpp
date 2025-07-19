#!/bin/bash

# Start time, and other metadata
start_date=$(date '+%Y-%m-%d_%H-%M')

# Add colors
red=1
green=2
yellow=3
blue=4
magenta=5
cyan=6
white=7
color() {
  tput setaf $1
  cat
  tput sgr0
}

# Parsing command line arguments
cli=""
if [ $# -ne 0 ]; then
  for arg in "$@"; do
    if [[ "$arg" == --results-folder=* ]]; then
      results_folder="${arg#--results-folder=}"
      printf "\nResults folder set to: $results_folder\n"
    # elif [ "$arg" == "--verbose" ]; then
    #   verbose=1
    #   printf "\VERBOSE mode turned on.\n"
    # elif [ "$arg" == "--flags" ]; then
    #   flags=1
    #   printf "\FLAGS mode turned on.\n"
    else
      cli="$cli $arg"
    fi
  done
fi

# Create results directory and change to it
if [ -n "$results_folder" ]; then
  results_dir="$results_folder"
else
  results_dir="results_${start_date}"
fi
mkdir -p "${results_dir}"
cd "${results_dir}"

# Function to benchmark map with csv
benchmark() {
  local map=$1
  local csv=$2
  local curr_cli=$cli

  if [[ "${country}" == world* ]]; then
    curr_cli="${cli} --world"
  fi

  cmd="cartogram ${map} ${csv} ${curr_cli}"
  csv_base="${csv_name%.csv}"
  # hyperfine '$cmd' --export-csv tmp.csv --command-name $csv_base
  hyperfine --parameter-scan qlcf 5 15 "${cmd} --quadtree_leaf_count_factor \$((2**{qlcf}))" --export-csv tmp.csv $(echo "--command-name ${csv_base}-"{05..15})

  if [ ! -f final_results.csv ]; then
    cp tmp.csv final_results.csv
  else
    tail -n +2 tmp.csv >> final_results.csv
  fi
}

# Iterating through folders in sample_data/
for folder in ../../sample_data/*; do
  if [[ -d "${folder}" && "${folder}" != *"sandbox"* ]]; then
    country=${folder##*/}
    echo " -------- Benchmarking ${country}" | color $magenta

    # Iterating through maps (GeoJSON(s) or JSON(s)) in country's folder
    for map in ${folder}/*.*json; do
      map_file_name=${map##*/}
      echo "Trying ${map_file_name}" | color $blue

      # Iterating through visual variable files (CSV(s)) in country's folder
      for csv in ${folder}/*.csv; do
        csv_name=${csv##*/}
        echo " - with ${csv_name}" | color $cyan
        if [[ -f "${map}" && -f "${csv}" ]]; then
          benchmark "${map}" "${csv}"
        else
          echo "Skipping invalid file: ${map} or ${csv}" | color $red
        fi
      done
    done
    echo " -------- Finished benchmarking ${country}." | color $yellow
  fi
done

# Clean up the temporary file
rm tmp.csv

# Print full file location
echo "Results saved to: $(pwd)/final_results.csv" | color $green