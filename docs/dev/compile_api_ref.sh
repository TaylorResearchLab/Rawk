#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")/../api_reference"

printf "Compile API reference\n"

python ../dev/docstrings_to_md.py

python ../dev/md_to_html_pdf.py API.md API.pdf API.html

printf "done\n"
