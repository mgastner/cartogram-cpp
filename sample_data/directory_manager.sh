#!/usr/bin/env bash

shopt -s extglob

for file in *.*json; do
  country="${file%%.*json}"
  country="${country%%_*}"
  echo "Working on ${country}..."
  mkdir -p ./${country}/
  mv ${file} ./${country}/${file}
  echo -e "# Sources\n" >> ./${country}/${country}_metadata.md
  echo -e "## ${file}\n\n\n" >> ./${country}/${country}_metadata.md
  for csv in ${country}*.csv; do
    if [ -f "${csv}" ]; then
      mv ${csv} ./${country}/${csv}
      echo -e "### ${csv}\n\n\n" >> ./${country}/${country}_metadata.md
    fi
  done
done
