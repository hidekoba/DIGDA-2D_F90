#!/bin/bash

ncks=${NCKS:-ncks}  # Install NCO and ensure ncks is in PATH, or set NCKS

# Check if the directory is passed as an argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 directory"
    exit 1
fi

dir=$1

# Check if the specified directory exists
if [ ! -d "$dir" ]; then
    echo "Directory $dir does not exist."
    exit 1
fi

# Create an array of files
shopt -s nullglob
file_list=("$dir"/SEDTRR[0-1]*.nc)

# Check if there are any files to process
if [ "${#file_list[@]}" -eq 0 ]; then
    echo "No files found matching the pattern SEDTRR[0-1]*.nc in $dir"
    exit 1
fi

# Combine the NetCDF files using ncks
output_file="$dir/SEDTRR.nc"
first_file="${file_list[0]}"
cp "$first_file" "$output_file"

# Append each file to the output file
for file in "${file_list[@]:1}"; do
    $ncks -A "$file" "$output_file"
    if [ $? -ne 0 ]; then
        echo "ncks command failed on $file"
        exit 1
    fi
done

echo "NetCDF files combined successfully into $output_file"
# The original script deleted its input files after combination.
# Keep inputs in the public example; remove them manually after validation.
