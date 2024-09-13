#!/usr/bin/env bash

# Start time, and other metadata
start_date=$(date '+%Y-%m-%d_%H-%M')

# Create results directory and change to it
results_dir="results_${start_date}"
mkdir -p "${results_dir}"
cd "${results_dir}"

results_file="log.txt"
successful_runs="successful_runs.txt"
failed_runs="failed_runs.txt"
tmp_file="tmp_file.txt"
SECONDS=0


# Arrays to store runtimes and map-csv pairs
runtimes=()
map_csv_pairs=()

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

# Flag to indicate if the script should skip the current test

# Function to handle SIGINT (Ctrl + C)
handle_sigint() {
  printf "\n\nSkipping current test...\n\n" | color $red
}

# Trap SIGINT (Ctrl + C)
trap handle_sigint SIGINT

# Parsing command line arguments
# From https://unix.stackexchange.com/questions/462515/
# check-command-line-arguments
if [ $# -eq 0 ] || [ -z "$1" ]; then
  cli="-QST"
else
  cli="$1"
fi

printf "\nWriting to ${results_file}\n"
printf "Tested on ${start_date}\n" >> "${results_file}"

printf "Passing" | tee -a "${results_file}"
printf " ${cli} " | tee -a "${results_file}" | color $magenta
printf "to cartogram\n\n" | tee -a "${results_file}"

# Turning on extended glob option
shopt -s extglob nocasematch

# Simple progress bar:
# Inspired from: https://github.com/pollev/bash_progress_bar

PROGRESS_BAR_WIDTH=50 # default progress bar length in characters

draw_progress_bar() {

  local __percentage=$((10#$1))

  # Rescale the bar according to the progress bar width
  local __num_bar=$(($__percentage * $PROGRESS_BAR_WIDTH / 100))

  # Draw progress bar
  for b in $(seq 1 $__num_bar); do printf "⬜"; done
  for s in $(seq 1 $(($PROGRESS_BAR_WIDTH - $__num_bar))); do printf " "; done
  printf "$__percentage%% \r"
}

# Function to test map with csv
run_map() {

  COLUMNS=$(tput cols)
  if [[ "$COLUMNS" -gt 60 ]]; then
    PROGRESS_BAR_WIDTH=50
  elif [[ "$COLUMNS" -gt 35 ]]; then
    PROGRESS_BAR_WIDTH=25
  else
    PROGRESS_BAR_WIDTH=0
  fi

  printed=0
  integration_count=0
  max_area_err=""
  start=$SECONDS
  # if world map then use -W flag
  curr_cli=$cli
  if [[ "${country}" == world* ]]; then
    curr_cli="${cli} -W"
  fi

  printf "now trying: \ncartogram ${map} ${csv} ${curr_cli}\n\n"
  echo "cartogram ${map} ${csv} ${curr_cli}" | pbcopy
  draw_progress_bar 0
  stdbuf -oL cartogram ${map} ${csv} ${curr_cli} 2>&1 | \
  # "cartogram" ${map} ${csv} ${cli} 2>&1 |
    while read line; do
      # save to temp file
      echo $line >>${tmp_file}

      # display warnings, errors etc.
      if [[ $line =~ "error" || $line =~ "warning" || $line =~ "invalid" ]]; then
        printf "\n\n%s\n\n" "$line" | tee -a "${results_file}" | color $red
      fi

      if [[ $line =~ "progress: 0." ]]; then
        draw_progress_bar ${line:12:2}
      fi

      # Check if integration is mentioned and update counter
      if [[ $line =~ "Integration number" ]]; then
        integration_count=$((integration_count + 1))
        printf "%s\n" "$line" | color $red
      fi

      if [[ $line =~ "New grid dimensions:" ]]; then
        printf "%s\n" "$line" | color $blue
      fi

      # Check if max area error is mentioned and update the variable
      if [[ $line =~ "max. area err:" ]]; then
        max_area_err=$line
        printf "%s\n" "$line" | color $yellow
      fi

      if [[ $line =~ "average area err:" ]]; then
        printf "%s\n" "$line" | color $magenta
      fi

      # Check if integration finished
      if [[ $line =~ "progress: 1" && "$printed" -eq 0 ]]; then
        draw_progress_bar 100
        printed=1
        printf "\n\n== Integration finished ==\nTotal integrations done: %d\n%s\n" "$integration_count" "$max_area_err" | tee -a "${results_file}" | color $yellow
      fi
    done
  end=$SECONDS
  runtime=$((end - start))
  printf "== Runtime ${runtime}s == \n" | tee -a "${results_file}"

  # Checking for any errors, invalid geometry or unfinished integration
  if grep -qi "invalid" ${tmp_file} || grep -qi "error" ${tmp_file} || ! grep -Fxq "Progress: 1" ${tmp_file}; then
    printf "== FAILED ==\n" | tee -a "${results_file}" | color $red
    printf "cartogram ${map} ${csv} ${cli}\n" >> ${failed_runs}

    # Printing country to failed_tmp.txt
    if [[ "${failed}" -eq 0 ]] || ! grep -qi "${country}" failed_tmp.txt; then
      printf "\n${country}\n" >>failed_tmp.txt
    fi

    # Increasing failed counter
    failed=$((failed + 1))

    # Saving data to file if FAILED
    map_wo_ext=${map_file_name%.*json}
    csv_wo_ext=${csv_name%.csv}
    err_file="DNF-${map_wo_ext}-with-${csv_wo_ext}.txt"
    printf " - with ${csv_name} in ${runtime}s\n" >> failed_tmp.txt
    printf "cartogram %s %s %s\n" "$map" "$csv" "$cli" >> failed_tmp.txt
    printf "Full output saved ro ${err_file}\n" | tee -a "${results_file}"
    mv ${tmp_file} ${err_file}

    # If no errors, and integration finished, pass
  else
    printf "== PASSED ==\n" | tee -a "${results_file}" | color $green
    printf "cartogram ${map} ${csv} ${cli}\n" >> ${successful_runs}

    # Store the runtime for the current map
    runtimes+=($runtime)
    map_csv_pairs+=("${map_file_name} with ${csv_name}")

  fi

  # Empty temporary file
  >${tmp_file}

  # Printing new line
  printf "\n" | tee -a "${results_file}"
}

# Counting number of tests
countries=0
total_tests=0
failed=0

# Iterating through folders in ..sample_data/
for folder in ../../sample_data/*; do
  if [[ -d "${folder}" && "${folder}" != *"sandbox"* ]]; then
    countries=$((countries + 1))
    country=${folder##*/}
    printf " -------- Testing ${country}\n\n" | tee -a "${results_file}" | color $magenta

    # Iterating through maps (GeoJSON(s) or JSON(s)) in country's folder
    for map in ${folder}/*.*json; do
      map_file_name=${map##*/}
      printf "Trying ${map_file_name}\n" | tee -a "${results_file}" | color $blue

      # Iterating through visual variable files (CSV(s)) in country's folder
      for csv in ${folder}/*.csv; do
        total_tests=$((total_tests + 1))
        csv_name=${csv##*/}
        printf " - with ${csv_name}\n\n" | tee -a "${results_file}" | color $cyan

        # Running cartogram, and timing it in seconds
        run_map

      done
    done
    printf " -------- Finished testing ${country}.\n\n\n\n\n\n" | tee -a "${results_file}" | color $magenta
  fi
done

# Summary report
printf "===== Finished testing all countries. =====\n\n" | tee -a "${results_file}" | color $magenta

duration="$(($SECONDS / 60))m $(($SECONDS % 60))s"
total_runtime=0
max_runtime=0
max_runtime_map=""
median_runtime=0

# Calculate total runtime, max runtime, and find the map with max runtime
for i in "${!runtimes[@]}"; do
  total_runtime=$((total_runtime + runtimes[i]))
  if (( runtimes[i] > max_runtime )); then
    max_runtime=${runtimes[i]}
    max_runtime_map=${map_csv_pairs[i]}
  fi
done

# Calculate average runtime
average_runtime=$((total_runtime / total_tests))

# Calculate median runtime
sorted_runtimes=($(printf '%s\n' "${runtimes[@]}" | sort -n))
if (( total_tests % 2 == 0 )); then
  median_runtime=$(( (sorted_runtimes[total_tests/2 - 1] + sorted_runtimes[total_tests/2]) / 2 ))
else
  median_runtime=${sorted_runtimes[total_tests/2]}
fi

# Calculate the top quartile runtime maps
# top_quartile_count=$((total_tests / 4))
# if (( top_quartile_count == 0 )); then
#   top_quartile_count=1
# fi
# sorted_indices=($(printf '%s\n' "${!runtimes[@]}" | sort -rn -k1,1 -t ' '))
# top_quartile_maps=()
# for i in $(seq 0 $((top_quartile_count - 1))); do
#   idx=${sorted_indices[$i]}
#   top_quartile_maps+=("${runtimes[$idx]}s - ${map_csv_pairs[$idx]}")
# done

printf "Finished ${total_tests} tests on ${countries} countries in ${duration}.\n"
printf "Total runtime: ${total_runtime}s\n" | tee -a "${results_file}"
printf "Average runtime: ${average_runtime}s\n" | tee -a "${results_file}"
printf "Max runtime: ${max_runtime}s (Map: ${max_runtime_map})\n" | tee -a "${results_file}"
printf "Median runtime: ${median_runtime}s\n" | tee -a "${results_file}"

# printf "Top quartile maps with longest runtime:\n" | tee -a "${results_file}"
# for map in "${top_quartile_maps[@]}"; do
#   printf " - ${map}\n" | tee -a "${results_file}"
# done

failed_per=$((100 * $failed / $total_tests))
printf "Passed [$((total_tests - failed))/${total_tests}] | $((100 - failed_per))%% \n" | tee -a "${results_file}" | color $green
printf "Failed [${failed}/${total_tests}] | ${failed_per}%% \n" | tee -a "${results_file}" | color $red

# Checking if any tests failed
if [[ "${failed}" -gt 0 ]]; then

  # Showing failed tests and removing temporary file
  printf "\nFailed tests:\n" | tee -a "${results_file}" | color $magenta

while IFS= read -r line; do
    if [[ $(echo $line | wc -w) -eq 1 ]]; then
      printf "$line\n" | tee -a "${results_file}" | color $red
    elif [[ $line == " - "* ]]; then
      printf "$line\n" | tee -a "${results_file}" | color $yellow
    else
      printf "$line\n" | tee -a "${results_file}"
    fi
  done < failed_tmp.txt
  rm failed_tmp.txt
  printf "\n" | tee -a "${results_file}"
fi

# Removing temporary files
rm ${tmp_file}

printf "Results saved to ${results_dir}\n\n"

# Prompting for file deletion
printf "Clear ALL *.geojson, *.png and *.ps files in directory? [y/N]: " | color $yellow
read -t 10 to_clear
if [[ $? -ne 0 ]]; then
  printf "\nNo response received! " | color $red
elif [ "$to_clear" == "y" ]; then
  rm *.geojson
  rm *.png
  rm *.ps
  printf "All *.geojson, *.png and *.ps files deleted.\n" | color $red
  exit ${failed}
fi

printf "Files not cleared.\n" | color $green
exit ${failed}
