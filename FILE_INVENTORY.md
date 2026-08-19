# Repository file inventory

Snapshot date: 2026-08-19

This is a read-only inventory of the working tree before any directory
reorganization. It excludes `.git` and reports apparent file sizes. The figures
describe the local working copy, including tracked, untracked, and ignored files.

## Summary

- Files: 1,329
- Total apparent size: 718.23 GB
- Git-tracked files: 161
- Visible untracked files before this policy update: 346
- Ignored untracked files before this policy update: 824
- Git object database: about 836 MB, including a 633.65 MB packed history

The working tree is dominated by generated simulation trajectories:

| File type | Files | Apparent size | Share of working tree |
|---|---:|---:|---:|
| `.lammpstrj` | 250 | 700.71 GB | 97.56% |
| `.data` | 387 | 13.50 GB | 1.88% |
| `.txt` | 168 | 2.16 GB | 0.30% |
| `.dat` | 18 | 1.12 GB | 0.16% |
| other files | 506 | about 0.74 GB | 0.10% |

Percentages are rounded, so they may not sum exactly to 100%.

## Largest top-level areas

| Area | Files | Apparent size |
|---|---:|---:|
| `adhesion/` | 619 | 460.70 GB |
| `pre_adhesion/` | 207 | 165.84 GB |
| `modulus/` | 84 | 14.46 GB |
| `friction/` | 65 | 12.53 GB |
| `melt/` | 57 | 11.38 GB |
| `ca_pdms/` | 28 | 10.09 GB |
| `900K_melt/` | 27 | 8.05 GB |
| `dif_is_melt/` | 29 | 8.01 GB |
| `600K_melt/` | 27 | 8.01 GB |
| `rp_melt/` | 32 | 7.42 GB |
| `solvent/` | 35 | 6.44 GB |

`adhesion/` and `pre_adhesion/` together account for approximately 626.54 GB,
or 87.2% of the working tree.

## Largest individual files

The largest files are generated `dump.*.lammpstrj` trajectories. The five
largest are approximately 21.37 GB, 19.86 GB, 16.55 GB, 14.83 GB, and 14.69 GB.
They occur under `pre_adhesion/N32bursh/` and `adhesion/brush_tersoff/`.

Large non-trajectory outputs also exist:

- `.data`: 13.50 GB across 387 structure files; these mix required inputs with
  generated snapshots and therefore cannot safely be classified by extension alone.
- `.txt`: 2.16 GB; examples include force maps around 25 MB each and a
  `modulus/melt_50/minp100.txt` file around 357 MB.
- `.dat`: 1.12 GB; most are generated `msd.dat` time series around 75 MB each.

## Interpretation

The repository itself is small enough for source-oriented Git collaboration, but
the local research workspace is not. Git should describe how simulations are
configured and analyzed; it should not become the primary store for full
trajectories or repeated intermediate snapshots.

This inventory is a point-in-time snapshot, not a generated manifest. No source
file, simulation directory, or result was moved or deleted while producing it.
