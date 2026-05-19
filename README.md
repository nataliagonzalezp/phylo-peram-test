# phylo.peram.test

**Phylogenetic Peramorphosis Test with Bootstrap Recovery**

`phylo.peram.test` is an R script for inferring heterochronic changes across a phylogeny using geometric morphometric data. The script implements a phylogenetic extension of the peramorphosis test of Piras et al. (2011), allowing users to classify ancestor-descendant branches as **peramorphic**, **paedomorphic**, or **not heterochronically changed** based on distances in multivariate shape space.

The method is designed for comparative studies of ontogenetic shape evolution and can be used for reproducibility, teaching, and methodological development.

---

## Authors

**Natalia González-Piñeres**  
PhD scholar
Unidad Ejecutora Lillo (UEL-FML-CONICET)

**Santiago A. Catalano**  
Investigador Independiente CONICET  
Unidad Ejecutora Lillo (UEL-FML-CONICET)  
FCNeIML-UNT

---

## Version

**v1.0**

Developed in:

```r
R version 4.5.3 (2026-03-11) "Reassured Reassurer"
```

---

## Citation

If you use this script, please cite:

```text
González-Piñeres N, Catalano SA. 2026. phylo.peram.test:
a phylogenetic peramorphosis test with bootstrap recovery. In preparation.
```

This script is based on and extends the logic of the peramorphosis test proposed by Piras et al. (2011).

Please also cite the relevant R packages used in your analysis, especially `geomorph`, `ape`, and `phytools`.

---

## Purpose

This script infers heterochronic events along all ancestor-descendant branches of a phylogenetic tree using landmark-based geometric morphometric data.

For each branch, the script compares the developmental displacement of the ancestor and descendant from a shared juvenile reference shape. Branches are classified as:

- **PERAMORPHOSIS**: the descendant is farther from the juvenile reference shape than its ancestor.
- **PAEDOMORPHOSIS**: the descendant is closer to the juvenile reference shape than its ancestor.
- **NC**: no change; the difference between ancestor and descendant distances is below the user-defined threshold.

Bootstrap resampling is used to estimate the recovery frequency of each observed classification.

---

## Methodological overview

The script performs the following major steps:

1. Reads landmark coordinates from a `.tps` file.
2. Reads a phylogenetic tree from a Newick `.tre` file.
3. Uses either an external developmental indicator, such as age or stage, or log-transformed centroid size.
4. Standardizes developmental progress within each species.
5. Estimates a shared juvenile reference shape at developmental onset.
6. Estimates species-specific maximum developmental shapes using Procrustes regression.
7. Reconstructs ancestral shapes using phylogenetic PCA and Maximum Likelihood ancestral state estimation.
8. Compares ancestor and descendant distances from the juvenile reference shape.
9. Classifies each branch as peramorphic, paedomorphic, or no change.
10. Performs bootstrap resampling to estimate recovery percentages.
11. Exports result tables and annotated phylogenetic plots.

---

## Input files

The script requires two input files:

1. a `.tps` file containing landmark coordinates; and
2. a `.tre` file containing the phylogenetic tree in Newick format.

Both files can be placed in the same working directory as the script, or full file paths can be specified in the user-defined parameters.

---

## TPS file

The `.tps` file must contain landmark coordinates and specimen identifiers.

The TPS file may contain:

- raw, unaligned landmark coordinates; or
- pre-aligned Procrustes coordinates.

The file may also include a developmental indicator field, such as:

```text
AGE=0.75
```

Example TPS block:

```text
LM=10
0.123 0.456
0.789 0.012
...
ID=Species_name
AGE=0.75
```

Specimen IDs must match the species names in the phylogenetic tree, or they must be modifiable within the script to match the tree tip labels.

Each species must have at least two specimens. Each species must also show variation in the developmental indicator used for regression.

---

## Tree file

The tree file must be provided in Newick format:

```text
your_tree.tre
```

Requirements:

- Tip labels must match the species names in the TPS file.
- The tree should preferably include meaningful branch lengths.
- If the tree has no branch lengths, the script assigns unit branch lengths and issues a warning.
- The tree should be fully resolved for ancestral reconstruction.

---

## User-defined parameters

The main parameters are located at the top of the script in the section **USER-DEFINED PARAMETERS**.

### Input files

```r
TPS_file <- "yourfile.tps"
tree_file <- "yourtree.tre"
```

### TPS field identifiers

```r
tps_specID <- "ID"
tps_devtime_keyword <- "AGE"
```

Use `tps_devtime_keyword` to define the TPS field containing the developmental indicator. This can be `"AGE"`, `"STAGE"`, `"CS"`, `"DEVTIME"`, or another user-defined keyword.

---

## Input mode flags

The script supports different combinations of coordinate and developmental-indicator inputs.

### Pre-aligned coordinates with age or stage labels

```r
coords_aligned <- TRUE
has_age_label <- TRUE
```

The script reads the aligned coordinates directly and uses the provided developmental indicator.

---

### Raw coordinates without age or stage labels

```r
coords_aligned <- FALSE
has_age_label <- FALSE
```

The script performs Generalized Procrustes Analysis using `gpagen()`, extracts centroid size, applies a logarithmic transformation, and uses `log(centroid size)` as the developmental progress indicator.

---

### Raw coordinates with age or stage labels

```r
coords_aligned <- FALSE
has_age_label <- TRUE
```

The script performs Generalized Procrustes Analysis but keeps the developmental indicator provided in the TPS file.

---

### Unsupported mode

```r
coords_aligned <- TRUE
has_age_label <- FALSE
```

This combination is not supported because centroid size cannot be recovered from already aligned coordinates. In this case, the script stops and asks the user to provide a developmental indicator.

---

## Analytical parameters

```r
resampling_proportion <- 0.66
diff_ratio_threshold <- 0.10
n_iter <- 100
set.seed(5422)
```

### `resampling_proportion`

Proportion of specimens retained per species in each bootstrap replicate.

Default:

```r
resampling_proportion <- 0.66
```

Species with exactly two specimens are always retained completely.

---

### `diff_ratio_threshold`

Threshold used to classify a branch as heterochronic or no change.

Default:

```r
diff_ratio_threshold <- 0.10
```

The script calculates:

```text
diff_ratio = |dist_A - dist_B| / mean(dist_A, dist_B)
```

Classification rule:

```text
If diff_ratio < threshold:
    NC

If diff_ratio >= threshold and dist_A > dist_B:
    PAEDOMORPHOSIS

If diff_ratio >= threshold and dist_A < dist_B:
    PERAMORPHOSIS
```

Higher thresholds make the test more conservative. Lower thresholds make it more sensitive.

---

### `n_iter`

Number of valid bootstrap replicates.

Default:

```r
n_iter <- 100
```

The script attempts additional bootstrap draws if some replicates fail validity checks.

---

### `set.seed`

The random seed is set to ensure reproducibility.

```r
set.seed(5422)
```

Change this value to use a different random stream, or remove it to run the analysis without a fixed seed.

---

## Required R packages

The following R packages are required:

```r
install.packages(c(
  "abind",
  "ape",
  "dplyr",
  "geiger",
  "geomorph",
  "ggplot2",
  "ggpubr",
  "ggthemes",
  "gridExtra",
  "phytools",
  "readr",
  "tidyr",
  "tidyverse",
  "viridis"
))
```

Some installations may also require:

```r
install.packages("devtools")
install.packages("Morpho")
```

Depending on the user's system and package versions, additional dependencies may be required.

---

## How to run the script

1. Clone or download this repository.

2. Place the script, TPS file, and tree file in the same working directory.

3. Open `phylo.peram.test.R` in R or RStudio.

4. Edit the user-defined parameters at the top of the script:

```r
TPS_file <- "yourfile.tps"
tree_file <- "yourtree.tre"

coords_aligned <- TRUE
has_age_label <- TRUE

resampling_proportion <- 0.66
diff_ratio_threshold <- 0.10
n_iter <- 100
```

5. Run the script:

```r
source("phylo.peram.test.R")
```

6. Inspect the output files written to the working directory.

---

## Analytical workflow

The script executes the following workflow:

1. Cleans the R environment.
2. Loads the required libraries.
3. Reads user-defined parameters.
4. Extracts landmark coordinates and metadata from the TPS file.
5. Aligns coordinates when required.
6. Assigns the developmental progress indicator.
7. Reads the phylogenetic tree.
8. Checks consistency between species names and tree tip labels.
9. Standardizes developmental progress within each species.
10. Estimates the shared juvenile reference shape.
11. Estimates maximum developmental shapes for each species.
12. Reconstructs ancestral shapes.
13. Infers observed heterochronic trends for all tree branches.
14. Plots the observed heterochronic classifications.
15. Performs bootstrap resampling.
16. Aggregates bootstrap results.
17. Exports CSV tables and PDF plots.

---

## Output files

With default parameters, the script writes the following files:

```text
HeterochronicChanges_10_DThreshold.pdf
HeterochronicChanges_Bootstrap66_10_DThreshold.pdf
66Keep10Threshold_TotalSample_Heterochrony_Results_WithRecovery.csv
66Keep10Threshold_RawBootstrapTrends.csv
```

### Main results table

The main results CSV includes:

- ancestor node;
- descendant node or species;
- Procrustes distance from ancestor to juvenile reference shape;
- Procrustes distance from juvenile reference shape to descendant;
- distance-ratio value;
- inferred heterochronic trend;
- abbreviated trend label;
- bootstrap recovery percentage.

### Raw bootstrap table

The raw bootstrap CSV includes the inferred trend for every branch in every valid bootstrap replicate.

### PDF plots

The script produces two annotated phylogenetic plots:

1. observed heterochronic classifications; and
2. heterochronic classifications with bootstrap recovery percentages.

Color code:

```text
Blue  = PERAMORPHOSIS
Red   = PAEDOMORPHOSIS
Gray  = No Change
```

---

## Bootstrap recovery

Bootstrap recovery is calculated relative to the observed trend from the full dataset.

For each branch:

1. the observed trend is inferred from the complete dataset;
2. each valid bootstrap replicate is classified independently;
3. a replicate is counted as recovered if it produces the same trend as the observed analysis;
4. recovery percentage is calculated as:

```text
Recovery percentage = matching replicates / valid replicates * 100
```

Recovery measures the consistency of the classification under resampling. It does not directly measure the magnitude of the heterochronic change.

---

## Important notes

### Tree branch lengths

The original tree is used for all calculations. If the tree lacks branch lengths, the script assigns unit branch lengths and issues a warning.

Equal branch lengths are used only for plotting clarity and do not affect the calculations.

### Developmental indicator

The developmental indicator can be age, stage, log-transformed centroid size, or another biologically meaningful proxy of developmental progress.

When log-transformed centroid size is used, the input coordinates must be raw and unaligned, because centroid size cannot be recovered from Procrustes-aligned coordinates.

### Species names

Species names in the TPS file and tree tip labels must match exactly. If the TPS identifiers include specimen-level suffixes, the script can be modified to extract the species name.

Example:

```r
Species <- sub("_IND_[0-9]+$", "", ID)
```

### Bootstrap validity

A bootstrap replicate is discarded if any species loses developmental variation, if a regression fails, if ancestral reconstruction fails, or if the replicate does not produce classifications for all branches.

Species represented by exactly two specimens are not subsampled; both specimens are retained in every bootstrap replicate.

---

## Troubleshooting

### `coords_aligned = TRUE but has_age_label = FALSE`

This combination is not supported. Provide an age, stage, or developmental indicator in the TPS file, or use raw coordinates with:

```r
coords_aligned <- FALSE
has_age_label <- FALSE
```

---

### `The number of AGE values does not match the number of specimens`

Some specimens are missing the developmental indicator field, or the keyword defined in `tps_devtime_keyword` does not match the TPS file.

---

### `The following tree tips are missing from the morphometric data`

Tree tip labels and TPS species names do not match. Check spelling, underscores, capitalization, and suffixes.

---

### `Species X has fewer than 2 specimens`

Each species must have at least two specimens for the regression-based workflow.

---

### `Species X has no variation in Devtime`

Each species must include specimens with variation in developmental progress.

---

### `No valid bootstrap replicates were obtained`

This usually occurs when species sample sizes are small or when resampling removes developmental variation. Try increasing:

```r
resampling_proportion <- 0.80
```

or use a smaller exploratory value of `n_iter`.

---

### Bootstrap runs slowly

Large datasets or high bootstrap failure rates can slow down the analysis. For exploratory runs, try:

```r
n_iter <- 10
```

Once the workflow is running correctly, increase `n_iter` for the final analysis.

---

## Suggested repository structure

```text
phylo.peram.test/
├── README.md
├── LICENSE
├── phylo.peram.test.R
└── docs/
    └── phylo.peram.test.v1.0_DOCUMENTATION.txt
```

If example data are later distributed with the repository, the following structure can be used:

```text
phylo.peram.test/
├── README.md
├── LICENSE
├── phylo.peram.test.R
├── docs/
│   └── phylo.peram.test.v1.0_DOCUMENTATION.txt
└── example/
    ├── example_data.tps
    └── example_tree.tre
```

---

## License

This repository is distributed under the MIT License.

```text
MIT License

Copyright (c) 2026 Natalia González-Piñeres and Santiago A. Catalano
```

This software is provided as a research tool to support reproducibility, teaching, and methodological development. Although care has been taken to ensure that the code performs as described, the authors provide it without any express or implied warranty. Users are responsible for verifying the suitability of the software for their own data, analyses, and research purposes. The authors are not liable for any errors, misuse, data loss, or consequences arising from the use, modification, or distribution of this software.

See the `LICENSE` file for the full license text.

---

## Contact

For questions, bug reports, or suggestions, please use the GitHub Issues section of this repository.
