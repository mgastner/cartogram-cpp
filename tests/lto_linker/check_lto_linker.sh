#!/usr/bin/env bash
# Checks the selection of the linker for LTO (cmake/LtoLinker.cmake) with
# clang on Linux. Each case configures, builds, and runs the small project in
# this folder, and compares the ld.lld that CMake selected with the one that
# clang ran:
#
# - By default, and when -B, COMPILER_PATH, or CMAKE_PROGRAM_PATH point to
#   the directory of an older lld, the build must use LTO with the ld.lld
#   next to the compiler.
# - When the older ld.lld is forced with CARTOGRAM_LTO_LD_LLD, the build must
#   succeed: with LTO and that ld.lld if the test link accepts it, and
#   otherwise without LTO.
#
# Usage: check_lto_linker.sh <clang++> <directory of an older ld.lld>
# Example: check_lto_linker.sh clang++-22 /usr/lib/llvm-21/bin

set -uo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <clang++> <directory of an older ld.lld>" >&2
  exit 2
fi
cxx=$(command -v "$1") || { echo "ERROR: $1 not found" >&2; exit 2; }
old_dir=$2
if [[ ! -x "$old_dir/ld.lld" ]]; then
  echo "ERROR: no ld.lld in $old_dir" >&2
  exit 2
fi

source_dir=$(cd "$(dirname "$0")" && pwd)
adjacent="$(dirname "$(realpath "$cxx")")/ld.lld"
work=$(mktemp -d)
failures=0

# check_case <name> <expected ld.lld, or "consistent"> [VAR=value ...] -- [cmake option ...]
check_case() {
  local name=$1 expected=$2
  shift 2
  local env_vars=()
  while [[ $# -gt 0 && $1 != "--" ]]; do
    env_vars+=("$1")
    shift
  done
  shift

  local build="$work/$name" log="$work/$name.log" status=0
  env ${env_vars[@]+"${env_vars[@]}"} cmake -S "$source_dir" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="$cxx" \
    -DCMAKE_EXE_LINKER_FLAGS=-v "$@" > "$log" 2>&1 &&
    env ${env_vars[@]+"${env_vars[@]}"} cmake --build "$build" --verbose >> "$log" 2>&1 &&
    "$build/lto_linker_check" >> "$log" 2>&1 || status=$?

  # CMake reports the selected ld.lld; with -v, clang prints each command
  # that it runs, and the link command starts with the quoted path of the
  # linker
  local selected used verdict=ok
  selected=$(sed -n 's/^-- LTO linker: //p' "$log")
  used=$(grep -o -E '^ "[^"]*ld\.lld"' "$log" | tr -d ' "' | sort -u)

  if [[ $status -ne 0 ]]; then
    verdict="FAILED (exit code $status)"
  elif [[ $expected == "consistent" ]]; then
    [[ $used == "$selected" ]] || verdict="FAILED (CMake selected '$selected', clang ran '$used')"
  elif [[ $selected != "$expected" || $used != "$expected" ]]; then
    verdict="FAILED (expected $expected; CMake selected '$selected', clang ran '$used')"
  fi

  printf '%-20s LTO linker: %-32s clang ran: %-32s %s\n' \
    "$name" "${selected:-none}" "${used:-no ld.lld}" "$verdict"
  if [[ $verdict != ok ]]; then
    failures=$((failures + 1))
    tail -n 30 "$log"
  fi
}

check_case default "$adjacent" --
check_case B_prefix "$adjacent" -- "-DCMAKE_CXX_FLAGS=-B$old_dir"
check_case COMPILER_PATH "$adjacent" "COMPILER_PATH=$old_dir" --
check_case CMAKE_PROGRAM_PATH "$adjacent" -- "-DCMAKE_PROGRAM_PATH=$old_dir"
check_case forced_older_lld consistent -- "-DCARTOGRAM_LTO_LD_LLD=$old_dir/ld.lld"

rm -rf "$work"
if [[ $failures -ne 0 ]]; then
  echo "ERROR: $failures of 5 cases failed" >&2
  exit 1
fi
echo "All 5 cases passed"
