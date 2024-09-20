#!/usr/bin/env bash

# Function to add color
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

# Initialize variables
printed=0
integration_count=0
max_area_err=""

# Read lines from standard input
while read line; do
  case "$line" in
  *"error"* | *"warning"* | *"invalid"*)
    printf "\n\n%s\n\n" "$line" | color $red
    ;;
  *"progress: 0."*)
    draw_progress_bar ${line:12:2}
    ;;
  *"Integration number"*)
    integration_count=$((integration_count + 1))
    printf "%s\n" "$line" | color $red
    ;;
  *"New grid dimensions:"*)
    printf "%s\n" "$line" | color $blue
    ;;
  *"Max. area err:"*)
    max_area_err=$line
    printf "%s\n" "$line" | color $yellow
    ;;
  *"Area drift:"*)
    printf "%s\n" "$line" | color $magenta
    ;;
  *"progress: 1"*)
    if [[ "$printed" -eq 0 ]]; then
      draw_progress_bar 100
      printed=1
      printf "\n\n== Integration finished ==\nTotal integrations done: %d\n%s\n" "$integration_count" "$max_area_err" | color $yellow
    fi
    ;;
  *" ms")
    printf "%s\n" "$line" | color $cyan
    ;;
  *)
    printf "%s\n" "$line"
    ;;
  esac
done
