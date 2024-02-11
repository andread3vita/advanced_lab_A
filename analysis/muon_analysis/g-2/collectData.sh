#!/bin/bash

# Specify the directory containing your files
directory="../g-2/data"

# Specify the output file name
output_file="../total.txt"

# Navigate to the directory
cd "$directory" || exit

# Combine all files into a single file
cat * > "$output_file"

echo "Files have been concatenated into $output_file"
