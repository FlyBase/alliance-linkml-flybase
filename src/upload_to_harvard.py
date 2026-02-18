#!/usr/bin/env python3
"""Upload TSV files to Harvard FTP via SCP."""

import os
import sys
import subprocess
import glob
from pathlib import Path


def get_env_or_exit(var_name):
    """Get environment variable or exit with error."""
    value = os.environ.get(var_name)
    if not value:
        print(f"Error: {var_name} environment variable not set")
        sys.exit(1)
    return value


def upload_files():
    """Upload all TSV files to Harvard FTP via SCP."""
    # Get credentials from environment
    host = get_env_or_exit('HARVARD_FTP_HOST')
    user = get_env_or_exit('HARVARD_FTP_USER')
    remote_path = os.environ.get('HARVARD_FTP_PATH', '~/')

    # Find TSV files
    output_dir = os.environ.get('OUTPUT_DIR', './tsvs/')
    tsv_files = glob.glob(os.path.join(output_dir, '*.tsv'))

    if not tsv_files:
        print(f"No TSV files found in {output_dir}")
        sys.exit(1)

    print(f"Found {len(tsv_files)} TSV files to upload")

    # Upload each file via SCP
    for tsv_file in tsv_files:
        filename = Path(tsv_file).name
        remote_dest = f"{user}@{host}:{remote_path}/{filename}"

        print(f"Uploading {filename}...")
        result = subprocess.run(
            ['scp', tsv_file, remote_dest],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"Error uploading {filename}: {result.stderr}")
            sys.exit(1)

        print(f"  Done: {filename} uploaded successfully")

    print(f"\nAll {len(tsv_files)} files uploaded to {host}:{remote_path}")


if __name__ == '__main__':
    upload_files()
