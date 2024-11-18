#!/bin/bash

# Define the file path
file="Position_Scan_Data_by_GEM_V_and_DRIFT_V_Configuration.csv"

# Extract unique numbers after the slash in the header row
numbers=$(head -n 1 "$file" | tr ',' '\n' | grep -oE '/[0-9]+' | tr -d '/' | sort -u)

# Loop through each number and run the Python script
for number in $numbers; do
    python3 anal.py "$number"
done