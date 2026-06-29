#!/usr/bin/env bash

set -euo pipefail

cd "$(dirname "$0")/../tutorials"

printf "Run construct_recon3d_mrn\n"

sed -n '/^```python/,/^```/ { /^```/!p }' construct_recon3d_mrn.md \
  > construct_recon3d_mrn.py

python construct_recon3d_mrn.py

python ../dev/md_to_html_pdf.py construct_recon3d_mrn.md \
  construct_recon3d_mrn.pdf construct_recon3d_mrn.html



printf "Run example_mouse_data_analysis\n"

sed -n '/^```python/,/^```/ { /^```/!p }' example_mouse_data_analysis.md \
  > example_mouse_data_analysis.py

python example_mouse_data_analysis.py

python ../dev/md_to_html_pdf.py example_mouse_data_analysis.md \
  example_mouse_data_analysis.pdf example_mouse_data_analysis.html



printf "Run example_human_data_analysis\n"

sed -n '/^```python/,/^```/ { /^```/!p }' example_human_data_analysis.md \
  > example_human_data_analysis.py

python example_human_data_analysis.py

python ../dev/md_to_html_pdf.py example_human_data_analysis.md \
  example_human_data_analysis.pdf example_human_data_analysis.html



printf "done\n"
