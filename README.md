# KaplanMeier-for-OAK_POPLAR
KaplanMeier with input cleaning, QC and survival analysis

This script `Jieun_KM.R` runs Kaplan–Meier (KM) survival analyses on the OAK/POPLAR-ready TSV provided here in .gz form, producing:

- KM summary tables (median survival, N, events, survival at requested timepoints)
- Log-rank test statistics and p-values
- KM plots (PNG) with optional risk tables and p-values
- A per-endpoint MASTER TSV indexing all outputs

It is designed to be robust when column names are supplied as **strings** (e.g., `"OS"`, `"OS_event"`), and avoids using `.data` inside model formulas (which causes errors in `survival::survfit()` / `survival::survdiff()`).

---

## Inputs

### Required columns in the dataset

The script expects the input file to contain at least:

- `PtID` — patient identifier (one row per patient is assumed)
- Time-to-event columns (choose one endpoint per run):
  - `OS` and `OS.CNSR` for overall survival
  - `PFS` and `PFS.CNSR` for progression-free survival
- Grouping columns for stratified KM analyses (e.g. `Treatment`, `trial`, mutations, clinical variables)

**Censoring convention used by this dataset:**
- `*.CNSR == 1` means **censored**
- `*.CNSR == 0` means **event observed**

So the script derives:
- `OS_event  = 1 - OS.CNSR`
- `PFS_event = 1 - PFS.CNSR`

### File format

- Tab-separated values: `.tsv`
- Header row is required

Example:
```r
df <- read.table(
  gzfile("ML-ready-oak-poplar.tsv.gz"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)
