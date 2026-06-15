# Data inventory -- on the D: drive (NOT in this repo)

UK Biobank / Welsh-registry participant data stays on D: under the UKB DUA. This is the map of where it lives; the repo ships only the aggregate result tables in `results_public/`.

| Metric | Value |
|---|---|
| Data files inventoried | **12045** |
| Total size on D: | **7.0 GB** |
| PARTICIPANT (restricted, on D: only) | **739** |
| AGGREGATE (safe) | **11306** |
| Copied into results_public/ | **102** |

## By folder

| Folder | Files | Size | Participant | Aggregate |
|---|---|---|---|---|
| `assets/` | 6 | 14.2 MB | 6 | 0 |
| `data/` | 268 | 285.5 MB | 55 | 213 |
| `projects/` | 11771 | 6.7 GB | 678 | 11093 |

Full per-file map: `DATA_INVENTORY_ON_D.csv` (column `abs_path_on_D` = exact location on the D: drive).
