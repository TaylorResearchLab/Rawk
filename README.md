# Rawk

A python package for metabolic pathway local enrichment analysis via random
walks on metabolic reaction network.

## Install

```bash
pip install rawk
```

## Documentation

### Tutorials

The tutorials of Rawk are in the `docs/tutorials` folder,
which contains the following tutorials:

- `construct_recon3d_mrn.md`: Construct a metabolic reaction network from a
  genome scale metabolic model.
- `example_mouse_data_analysis.md`: Run Rawk standard analysis workflow on an
  example mouse dataset.
- `example_human_data_analysis.md`: Run Rawk standard analysis workflow on an
  example human dataset.

### API reference

The API reference files of Rawk are in
`docs/api_reference`. The API reference files were
generated from the package docstrings. The docstrings can also be accessed
using `help` in python interpreter. For example, `help(rawk.Rawk)` shows the
documentation of the `Rawk` class.

## Troubleshooting

If you encounter any error related to `tkinter` multi-threading, try rerunning
with parameters set to use only one CPU core.

## Notice about license

This project is released under a Non-Commercial Research License. For commercial use, please contact licensing@chop.edu for licensing terms.
