# Rawk: Metabolic reaction network local enrichment analysis

Rawk is a python package to identify and visualize locally enriched metabolic pathways in metabolic reaction networks.

## Install

**Note** that Rawk currently can only be installed from local source files.

Following is an example procedure for installing Rawk on a Linux computer using command line tools:

```bash
# Download rawk-0.0.5.tar.gz

tar -xvf rawk-0.0.5.tar.gz

cd rawk-0.0.5

pip install .
```

## Tutorials

The tutorials of Rawk are the python scripts in the `docs/tutorials` folder contained in `rawk-0.0.5.tar.gz`. The `docs/tutorials` folder contains the following tutorials:

- `construct_recon3d_met_net.py`: Construct a metabolic reaction network from a genome scale metabolic model.
- `example_mouse_data_analysis.py`: Run Rawk standard analysis workflow on an example mouse dataset.
- `example_human_data_analysis.py`: Run Rawk standard analysis workflow on an example human dataset.

Following is an example procedure for running the tutorials on a Linux computer using command line tools:

```bash
python construct_recon3d_met_net.py

python example_mouse_data_analysis.py

python example_human_data_analysis.py
```

## Documentation

Rawk package documentation can be accessed using `help` in python interpreter. For example, `help(rawk.Rawk)` shows the documentation of the `Rawk` class.

## Troubleshooting

If you encounter any error related to `tkinter` multi-threading, try rerunning with `workers=1` or removing `tkinter` from the environment.

## Notice about license

This project is released under a Non-Commercial Research License. For commercial use, please contact licensing@chop.edu for licensing terms.
