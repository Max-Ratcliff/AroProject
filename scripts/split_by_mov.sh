#!/usr/bin/env bash
# Scipt made to split the pulses_tht into managable folders
# move this script into the directory with the TIFF files
# and run it to organize the files into subdirectories based on their names.
set -euo pipefail

# Enable nullglob so that *.tif expands to nothing if there are no matches
shopt -s nullglob

for f in *.tif; do
  # If somehow it's not a regular file, skip
  [[ -f "$f" ]] || continue

  # strip off the last "_<frame>.tif" to get the directory name
  dir="${f%_*}"

  # make the directory (no-op if it already exists)
  mkdir -p "$dir"

  # move the TIFF into it, preserving the filename
  mv -- "$f" "$dir/$f"
done