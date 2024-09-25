#!/bin/bash

# Updated Directories
BAM_DIR="your_BAM_files_directory"
OUT_DIR="your_output_directory"
TEMP_DIR="temp_directory"

# Ensure the output and temporary directories exist
mkdir -p "$OUT_DIR"
mkdir -p "$TEMP_DIR"

# Iterate over each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Check if the corresponding BAM index exists. If not, create it.
    if [[ ! -f "${BAM_FILE}.bai" ]]; then
        echo "Indexing $BAM_FILE..."
        samtools index "$BAM_FILE"
    fi

    # Extract the base name of the BAM file to name the output accordingly
    BASE_NAME=$(basename "$BAM_FILE" ".bam")

    # Create a temporary BAM file for primary alignments with MQ >= 60
    TEMP_BAM="$TEMP_DIR/${BASE_NAME}_primary_MQ30.bam"
    echo "Filtering for primary alignments with MQ >= 60 in $BAM_FILE..."
    samtools view -b -F 256 -F 2048 -q 60 "$BAM_FILE" > "$TEMP_BAM"

    # Index the temporary primary BAM file
    echo "Indexing temporary BAM file for $BASE_NAME..."
    samtools index "$TEMP_BAM"

    # Run idxstats on the temporary primary BAM file and save the output
    echo "Running idxstats on primary alignments for $BASE_NAME..."
    samtools idxstats "$TEMP_BAM" > "$OUT_DIR/${BASE_NAME}_primary_MQ30_idxstats.txt"

    # Remove the temporary primary BAM file and its index
    rm "$TEMP_BAM"
    rm "${TEMP_BAM}.bai"
done

echo "Done."
