#!/bin/bash
# Upload TSV files to Harvard FTP via SCP

set -e

# Check required env vars
: "${HARVARD_FTP_HOST:?HARVARD_FTP_HOST not set}"
: "${HARVARD_FTP_USER:?HARVARD_FTP_USER not set}"
HARVARD_FTP_PATH="${HARVARD_FTP_PATH:-~/}"
OUTPUT_DIR="${OUTPUT_DIR:-./tsvs/}"

# Count files
file_count=0
for file in "$OUTPUT_DIR"/*.tsv; do
    if [ -f "$file" ]; then
        ((file_count++))
    fi
done

if [ "$file_count" -eq 0 ]; then
    echo "No TSV files found in $OUTPUT_DIR"
    exit 1
fi

echo "Found $file_count TSV files to upload"

# Find and upload TSV files
for file in "$OUTPUT_DIR"/*.tsv; do
    if [ -f "$file" ]; then
        echo "Uploading $(basename "$file")..."
        scp "$file" "${HARVARD_FTP_USER}@${HARVARD_FTP_HOST}:${HARVARD_FTP_PATH}/"
        echo "  Done"
    fi
done

echo ""
echo "All $file_count files uploaded to ${HARVARD_FTP_HOST}:${HARVARD_FTP_PATH}"
