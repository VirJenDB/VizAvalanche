# VirJenDB Visualizations – R
**Author:** Mohit Kumar

All data is fetched live from the [VirJenDB public API](https://api2.virjendb.org/v2/search) — no CSV download needed.

## Plots

| # | Plot | Method |
|---|------|--------|
| 1 | GC Content by Cluster Role | diverse API sample, log Y |
| 2 | Sequence Length by Viral Family | per-family API sample, log Y |
| 3 | Submissions Over Time | exact count per year via API, log Y |
| 4 | World Map by Collection Country | exact count per country via API |
| 5 | Molecule Type Distribution | exact count per type via API, log X |

## Run

```bash
conda env create -f environment.yml
conda activate virjendb_r_viz
Rscript virjendb_plots.R
```

Outputs: five interactive `.html` files and one `virjendb_plots.pdf`, saved in this folder.
