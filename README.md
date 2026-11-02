# Gyrobody geometry code

Custom code for the computational analyses in:

Bykowski, M. et al. *Mature plant chloroplasts can form gyroid cubic membranes.* Nature Communications (2026).

The repository holds only code. The measurements it runs on (TEM and light-microscopy measurements, GRANA output) are in the dataset https://doi.org/10.58132/VNNZV4.

## Scripts

| Script | What it calculates | Paper |
| --- | --- | --- |
| `gyrobody_prevalence.py` | Monte-Carlo calibration of the gyrobody-to-chloroplast size ratio κ, with the gyrobody either central or randomly placed in the chloroplast, and the detection probability p of a gyrobody in a 1 µm section | Supplementary Fig. 9; Methods "Gyrobody prevalence estimation" |
| `gyrobody_membrane_fraction.py` | Per chloroplast: chloroplast and gyrobody volumes, stroma thylakoid membrane area S_st, gyrobody membrane area S_gb and the gyrobody share f_gb | Supplementary Fig. 12e; Methods "Gyrobody volume and membrane-area estimation" |
| `gyrobody_membrane_area.py` | Unit-cell membrane area A_unit; total membrane area of an average gyrobody S_gb; membrane area of an average granum S_gs; number of grana equivalent to one gyrobody | Fig. 2c; Methods "SPIRE modeling and membrane geometry analysis" and "Gyrobody volume and membrane-area estimation" |
| `gaussian_curvature_tpms.py` | Gaussian curvature of nodal gyroid and diamond surfaces by the bordered-Hessian formula; curvature extremes and a per-vertex curvature mesh (PLY) | Fig. 1g; Methods "SPIRE modeling and membrane geometry analysis" |

The formulas are those given in the Methods sections listed above.

## System requirements

- Python 3.12 (tested with 3.12.14 on macOS 26; any OS with Python ≥ 3.10 should work)
- `numpy` (tested 2.3.5); `scikit-image` (tested 0.26.0, only for `gaussian_curvature_tpms.py`)
- No special hardware.

## Installation

```bash
python -m pip install -r requirements.txt
```

Installation takes under a minute.

## Input files

All inputs are CSV files with a header row. Lengths are in nm and areas in nm².

| Script | Option | Columns (one row per …) |
| --- | --- | --- |
| `gyrobody_prevalence.py` | `--tem-axes` | `genotype, chloroplast_major_nm, gyrobody_major_nm` (TEM chloroplast section) |
| | `--lm-counts` (optional) | `genotype, chloroplasts_with_gyrobody, chloroplasts_total` (light-microscopy cell) |
| `gyrobody_membrane_fraction.py` | `--chloroplasts` | `chloroplast_id, genotype, chloroplast_major_nm, chloroplast_minor_nm, chloroplast_area_nm2` (chloroplast section) |
| | `--gyrobodies` | `chloroplast_id, gyrobody_major_nm, gyrobody_minor_nm, gyrobody_area_nm2` (gyrobody profile) |
| | `--stroma-traces` | `chloroplast_id, length_nm` (traced stroma thylakoid profile) |
| `gyrobody_membrane_area.py` | `--gyrobodies` | `chloroplast_id, gyrobody_major_nm, gyrobody_minor_nm` (gyrobody profile; rows with the same `chloroplast_id` are summed) |
| | `--grana` | `granum_diameter_nm, n_thylakoids` (granum from GRANA; empty cells are skipped) |

`gaussian_curvature_tpms.py` needs no input.

## Options

| Script | Option | Default | Meaning |
| --- | --- | --- | --- |
| `gyrobody_prevalence.py` | `--position` | `central` | Gyrobody placement in the κ calibration: `central` (concentric with the chloroplast) or `random` (random position inside the chloroplast) |
| | `--thickness-nm` | 1000 | Section thickness t |
| | `--replicates` | 10 | Number of independent Monte-Carlo calibrations averaged; their SD is printed as `p_MC_SD` |
| | `--seed` | 1 | First random seed (replicate i uses seed + i) |
| `gyrobody_membrane_fraction.py` | `--unit-cell-nm` | 519 | Gyroid lattice parameter a; V_uc = a³ |
| | `--unit-cell-area-um2` | 1.62 | SPIRE membrane area of one unit cell, S_uc (both membranes) |
| `gyrobody_membrane_area.py` | `--lattice-nm` | 424 551 | Lattice parameters for which A_unit is printed |
| | `--lumen-nm` | 30 | Membrane–membrane spacing 2d |
| | `--unit-cell-nm`, `--unit-cell-area-um2` | 519, 1.62 | As above |
| `gaussian_curvature_tpms.py` | `--structure` | – | `gyroid` or `diamond` (nodal approximations) |
| | `--unit-cell-nm` | – | Lattice parameter a |
| | `--t-level` | 0 | Level-set value t (0 for the gyroid, 0.43 for the diamond PLB) |
| | `--out-dir` | `curvature_output` | Folder for the PLY mesh (vertex coordinates in nm, K in nm⁻², R = 1/√\|K\| in nm) |

Results are printed to the terminal as tab-separated tables; redirect them to a file to save them (`> results.tsv`).

## Running

```bash
# Supplementary Fig. 9 (both gyrobody placements)
python gyrobody_prevalence.py --tem-axes tem_axes.csv --lm-counts lm_counts.csv --position central
python gyrobody_prevalence.py --tem-axes tem_axes.csv --lm-counts lm_counts.csv --position random

# Supplementary Fig. 12e
python gyrobody_membrane_fraction.py --chloroplasts chloroplasts.csv \
    --gyrobodies gyrobodies.csv --stroma-traces stroma_traces.csv > membrane_fraction.tsv

# Fig. 2c and A_unit
python gyrobody_membrane_area.py --lattice-nm 424 551 \
    --gyrobodies stn7_gyrobodies.csv --grana stn7_grana.csv

# Fig. 1g (stn7-1 and aba1-6 gyroid unit cells, PLB diamond unit cell)
python gaussian_curvature_tpms.py --structure gyroid  --unit-cell-nm 551
python gaussian_curvature_tpms.py --structure gyroid  --unit-cell-nm 424
python gaussian_curvature_tpms.py --structure diamond --unit-cell-nm 90 --t-level 0.43
```

## Expected output and run time

Run on the dataset measurements, the scripts print the values reported in the paper:

| Script | Output | Run time* |
| --- | --- | --- |
| `gyrobody_prevalence.py --position central` | p = 0.632 (Col-0 6 dD), 0.754 (*stn7-1* 16 hD), 0.698 (*aba1-6* 16 hD); observed fractions 0.687, 0.730, 0.648 | ~2 s |
| `gyrobody_prevalence.py --position random` | p = 0.619, 0.747, 0.686 | ~5 s |
| `gyrobody_membrane_fraction.py` | f_gb for 11 *aba1-6* and 10 *stn7-1* chloroplasts (values in Supplementary Fig. 12e) | < 1 s |
| `gyrobody_membrane_area.py` | A_unit = 1.06 µm² (a = 424 nm), 1.81 µm² (a = 551 nm); mean V_gb = 25.96 µm³; S_gb = 300.9 µm²; S_gs = 8.60 µm²; 34.98 grana per gyrobody | < 1 s |
| `gaussian_curvature_tpms.py` | R at the most curved point of the nodal surface = 101.26 nm (gyroid, a = 551 nm), 77.92 nm (gyroid, a = 424 nm), 12.92 nm (diamond, a = 90 nm); PLY mesh with K and R per vertex | ~5 s per surface |

\* On an Apple-silicon laptop.

The Monte-Carlo SD of p over the 10 replicates is ≤ 0.0005.

## Licence

GNU General Public License v3.0 or later (GPL-3.0-or-later); see `LICENSE`.

## How to cite

Cite the paper above and this code: Bykowski, M. et al. Code for "Mature plant chloroplasts can form gyroid cubic membranes". Dane Badawcze UW https://doi.org/[DOI] (2026).
