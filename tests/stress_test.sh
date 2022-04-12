#!/usr/bin/env bash

# Start time, and other metadata
start_date=$(date '+%Y-%m-%d_%H-%M-%S')
results_file="results_${start_date}.txt"
tmp_file="tmp_file.txt"
SECONDS=0
printf "\nWriting to ${results_file}\n"
printf "Tested on ${start_date}\n" >> "${results_file}"

# Add colors
red=1; green=2; yellow=3; blue=4; magenta=5; cyan=6; white=7
color() { tput setaf $1; cat; tput sgr0; }

# Parsing command line arguments
# From https://unix.stackexchange.com/questions/462515/
# check-command-line-arguments
if [ $# -eq 0 ] || [ -z "$1" ]; then
  cli="-epst"
else
  cli="$1"
fi
printf "Passing" | tee -a "${results_file}"
printf " ${cli} " | tee -a "${results_file}" | color $magenta
printf "to cartogram\n\n" | tee -a "${results_file}"

# Turning on extended glob option
shopt -s extglob nocasematch

# Simple progress bar:
# Inspired from: https://github.com/pollev/bash_progress_bar

PROGRESS_BAR_WIDTH=50  # default progress bar length in characters

draw_progress_bar() {

  local __percentage=$((10#$1))

  # Rescale the bar according to the progress bar width
  local __num_bar=$(( $__percentage * $PROGRESS_BAR_WIDTH / 100 ))

  # Draw progress bar
  for b in $(seq 1 $__num_bar); do printf "⬜" ; done
  for s in $(seq 1 $(($PROGRESS_BAR_WIDTH - $__num_bar))); do printf " "; done
  printf "$__percentage%% \r"
}

# Function to test map with csv
run_map()
{

  COLUMNS=$(tput cols)
  if [[ "$COLUMNS" -gt 60 ]]; then
    PROGRESS_BAR_WIDTH=50
  elif [[ "$COLUMNS" -gt 35 ]]; then
    PROGRESS_BAR_WIDTH=25
  else
    PROGRESS_BAR_WIDTH=0
  fi

  printed=0
  draw_progress_bar 0
  start=$SECONDS
  cartogram ${map} ${csv} ${cli} 2>&1 |
    while read line
    do
      # save to temp file
      echo $line >> ${tmp_file}

      # display warnings, errors etc.
      if [[ $line =~ "error"  || $line =~ "warning" || $line =~ "invalid" ]]; then
        printf "\n\n$line\n\n" | tee -a "${results_file}" | color $red
      fi

      if [[ $line =~ "progress: 0." ]]; then
        draw_progress_bar ${line:12:2}
      fi

      # check if integration finished
      if [[ $line =~ "progress: 1" && "$printed" -eq 0 ]]; then
        draw_progress_bar 100
        printed=1
        printf "\n\n== Integration finished ==\n" | tee -a "${results_file}" | color $yellow
      fi
    done
  end=$SECONDS
  runtime=$((end-start))
  printf "== Runtime ${runtime}s == \n" | tee -a "${results_file}"

  # Checking for any errors, invalid geometry or unfinished integration
  if grep -qi "invalid" ${tmp_file} || grep -qi "error" ${tmp_file} || ! grep -Fxq "Progress: 1" ${tmp_file} ; then
    printf "== FAILED ==\n" | tee -a "${results_file}" | color $red

    # Printing country to failed_tmp.txt
    if [[ "${failed}" -eq 0 ]] || ! grep -qi "${country}" failed_tmp.txt; then
      printf "\n${country}\n" >> failed_tmp.txt
    fi

    # Increasing failed counter
    failed=$((failed+1))

    # Saving data to file if FAILED
    map_wo_ext=${map_file_name%.*json}
    csv_wo_ext=${csv_name%.csv}
    err_file="results_${start_date}-${map_wo_ext}-${csv_wo_ext}.txt"
    printf " - ${map_file_name} with ${csv_name}\n" >> failed_tmp.txt
    printf "Full output saved to ${err_file}\n" | tee -a "${results_file}"
    mv ${tmp_file} ${err_file}

    # If no errors, and integration finished, pass
  else
    printf "== PASSED ==\n" | tee -a "${results_file}" | color $green
  fi

  # Empty temporary file
  > ${tmp_file}

  # Printing new line
  printf "\n" | tee -a "${results_file}"
}

# Counting number of tests
countries=0
total_tests=0
failed=0

# Iterating through folders in ..sample_data/
for folder in ../sample_data/*; do
  if [[ -d "${folder}" && "${folder}" != *"sandbox"* ]]; then
    countries=$((countries+1))
    country=${folder##*/}
    printf " -------- Testing ${country}\n\n" | tee -a "${results_file}" | color $magenta

    # Iterating through maps (GeoJSON(s) or JSON(s)) in country's folder
    for map in ${folder}/*.*json; do
      map_file_name=${map##*/}
      printf "Trying ${map_file_name}\n" | tee -a "${results_file}" | color $blue

      # Iterating through visual variable files (CSV(s)) in country's folder
      for csv in ${folder}/*.csv; do
        total_tests=$((total_tests+1))
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
printf "Finished ${total_tests} tests on ${countries} countries in ${duration}.\n"

failed_per=$(( 100 * $failed / $total_tests))
printf "Passed [$((total_tests-failed))/${total_tests}] | $((100-failed_per))%% \n" | tee -a "${results_file}" | color $green
printf "Failed [${failed}/${total_tests}] | ${failed_per}%% \n"  | tee -a "${results_file}" | color $red

# Checking if any tests failed
if [[ "${failed}" -gt 0 ]]; then

  # Showing failed tests and removing temporary file
  printf "\nFailed tests:\n"  | tee -a "${results_file}" | color $red
  cat failed_tmp.txt | tee -a "${results_file}"
  rm failed_tmp.txt
  printf "\n" | tee -a "${results_file}"
fi

# Removing temporary files
rm ${tmp_file}

# Prompting for file deletion
printf "Clear ALL *.geojson, *.png and *.ps files in current directory? [y/N]: " | color $yellow
read -t 10 to_clear
if [[ $? -ne 0 ]]; then
  printf "\nNo response received! " | color $red
elif [ "$to_clear" == "y" ]; then
  rm *.geojson;
  rm *.png;
  rm *.ps;
  printf "All *.geojson, *.png and *.ps files deleted.\n" | color $red
  exit ${failed}
fi

printf "Files not cleared.\n" | color $green
exit ${failed}
