#!/bin/bash

# Check if a directory is provided
if [ "$#" -ne 2 ];then
    echo "Usage: $0 <directory> <suffix>" 
    exit 1
fi

# Assigning the first argument as the directory.
DIR=$1
# Assigning the second argument as the suffix.
SUFFIX=$2

# Output file 
OUTPUT_FILE="${DIR}/sample_ids.txt"

# Check if the directory exists
if [ ! -d "$DIR" ]; then
    echo "Directory not found: $DIR"
    exit 1
fi

# Temporarily store unique sample names;
declare -A sample_names

# Process files 
for file in "$DIR"/*"$SUFFIX";do 
    if [ -f "$file" ]; then
        # basename is used to extract the filename without the path from a complete file path，
        # and optionally remove a specific suffix from the filename.
        base_name=$(basename "$file" "$SUFFIX")
        sample_name=${base_name%_[12]}
        sample_names["$sample_name"]=1
    fi
done

> "$OUTPUT_FILE"
for name in "${!sample_names[@]}"; do
    echo "$name" >> "$OUTPUT_FILE"
done

echo "Sample IDs have been written to $OUTPUT_FILE"
