#!/bin/bash

# Check if a file was provided
if [ $# -ne 1 ]; then
  echo "Usage: $0 <source_file.c>"
  exit 1
fi

SRC_FILE="$1"

# Check if the file exists
if [ ! -f "$SRC_FILE" ]; then
  echo "Error: File '$SRC_FILE' does not exist."
  exit 1
fi

# Extract the filename without extension
OUT_FILE="${SRC_FILE}.out"

# Compile the C file
echo " >> Compiling $OUT_FILE"
gcc "$SRC_FILE" -c -o "tmp.o" -O3 -g -fopenmp -Wall
g++ tmp.o parse_and_run_multiple.cpp -o "$OUT_FILE" -fopenmp -g -O3 -Wall

# Check if compilation was successful
if [ $? -ne 0 ]; then
  echo " >> COMPILATION FAILED."
  exit 1
fi

echo "Compilation successful. Running $OUT_FILE..."
echo "--------------------"
./"$OUT_FILE"
