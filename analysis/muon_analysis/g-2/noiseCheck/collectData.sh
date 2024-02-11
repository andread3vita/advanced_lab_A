#!/bin/bash

# Specify the directory containing your files
directory="../../g-2/data/"

# Specify the output file name
output_file="../../g-2/noiseCheck/total.txt"

# Navigate to the directory
cd "$directory" || exit

# Combine all files into a single file
cat * > "$output_file"

echo "Files have been concatenated into $output_file"
