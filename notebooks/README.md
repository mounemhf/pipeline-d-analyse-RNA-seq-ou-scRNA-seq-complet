# Notebooks

Exploratory notebooks only — the reproducible pipeline logic lives in
`workflow/scripts/` (plain Python, run by Snakemake).

- `01_exploration.ipynb` — interactive exploration of the pipeline outputs
  (clusters, annotations, marker genes). Run the pipeline first, then start
  Jupyter from an environment with scanpy (e.g. `workflow/envs/preprocess.yml`
  plus `jupyterlab`).
