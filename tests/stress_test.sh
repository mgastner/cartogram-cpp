#!/usr/bin/env bash

shopt -s extglob

start_date=$(date '+%Y-%m-%d_%H:%M:%S')
results="results_${start_date}.txt"
echo -e "Writing to ${results}\n"
echo -e "Tested on ${start_date}\n" >> "${results}"

for folder in ../sample_data/*; do
  if [ -d "${folder}" ]; then
    country=${folder##*/}
    echo -e "-------- Testing ${country}...\n" | tee -a "${results}"
    for map in ${folder}/*.*json; do
      map_file_name=${map##*/}

# echo -e "ting $ all countries" | te"e"$ults}"
      echo "Trying ${map_file_name} ..." | tee -a "${results}"
      for csv in ${folder}/*.csv; do
        csv_name=${csv##*/}
        echo -e "with ${csv_name} ...\n" | tee -a "${results}"
        start=$(date +%s)
        $(cartogram ${map} -V ${csv} -est > tmp.txt 2>&1) >> tmp.txt
        end=$(date +%s)
        runtime=$((end-start))
        grep -i "error" tmp.txt | tee -a "${results}"
        grep -i "warning" tmp.txt | tee -a "${results}"
        echo -e "\nRuntime: ${runtime}s" | tee -a "${results}"
        if grep -Fxq "Progress: 1" tmp.txt; then
          echo -e "== PASS ==\n" | tee -a "${results}"
        else
          echo -e "== FAIL ==" | tee -a "${results}"
          map_wo_ext=${map_file_name%.*json}
          csv_wo_ext=${csv_name%.csv}
          err_file="results_${start_date}-${map_wo_ext}-${csv_wo_ext}.txt"
          echo -e "Full output saved to ${err_file}" | tee -a "${results}"
          $(cp tmp.txt "${err_file}")
        fi
      done
    done
    echo -e "-------- Finished testing ${country}.\n\n\n\n\n" | tee -a "${results}"
  fi
done

echo -e "===== Finished testing all countries. =====\n"  | tee -a "${results}"

# Removing temporary file
rm tmp.txt
