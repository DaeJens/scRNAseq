#!/bin/bash

# Check if cellranger-8.0.0 is in the PATH
command -v cellranger-8.0.0 >/dev/null 2>&1

# If command exits with non-zero status code, cellranger is not found
if [ $? -ne 0 ]; then
  echo "Error: cellranger-8.0.0 is not found in your PATH."
  exit 1
fi

echo "Cellranger-8.0.0 is installed and available in your PATH."

# Take the directory with fastq files as user input from the command line
if [[ $# -ne 1 ]]; then
  echo "Error: Please provide the directory containing fastq files as the first argument."
  exit 1
fi

fq_dir = "$1"

# Check if the provided directory exists
if [[ ! -d "$fq_dir" ]]; then
  echo "Error: The specified directory containing your fastq files does not exist."
  exit 1
fi

# Check if any files with .fastq extension exist in the directory
if [[ ! $(find "$fq_dir" -type f -name '*.fastq*') ]]; then
  echo "Error: No files with .fastq extension found in the directory."
  exit 1
fi

# Ask the user if their scRNAseq reads are multiplexed
read -p "Are your scRNAseq reads multiplexed? (y/n)" is_multi

while [[ ! is_multi =~ [YyNn]$ ]]; do
    read -p "Invalid input. Are your reads multiplexed (y/n)? " is_multi
done

if [[ is_multi =~ [Yy]]$ ]; do
    echo "Multiplexed reads detected. Adjusting CellRanger pipeline for multiplexing."
    # cellranger multi
else
    echo "Non-multiplexed reads detected. Using default CellRanger pipeline."
    # cellranger count
fi