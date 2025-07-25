# Thesis can be located here: 
https://etd.ohiolink.edu/acprod/odb_etd/etd/r/1501/10?clear=10&p10_accession_num=ohiou1628858064984242


# Population Genetic Analyses with DAPC in R

This repository contains R code for conducting population structure and isolation-by-distance analyses using **Discriminant Analysis of Principal Components (DAPC)** and **Nei’s genetic distance**. These analyses support a thesis project examining population differentiation in bobcats.

## Methods Included

- **Clustering using `find.clusters()` from `adegenet`**
- **DAPC analysis with cross-validation via `optim.a.score()`**
- **Visualization with `scatter()` and `compoplot()`**
- **Genetic distance calculation using Nei's distance (`poppr`)**
- **Isolation by Distance testing with Mantel tests**

## File Structure

── Bobcat_2pop.txt.DAT # Genetic data input
├── Bobcat_3pop.DAT # Alternate dataset for IBD analysis
├── samp_names.csv # Sample name metadata
├── nei_coord.csv # Coordinate data for spatial analysis
├── dapc2.2.csv # DAPC assignment probabilities (output)
├── dapc3.2.csv # DAPC assignment probabilities (output)
└── DAPC_analysis.R # Main script (this repo)



## Requirements

Make sure the following R packages are installed:

```r
install.packages("adegenet")
install.packages("poppr")
```

Running the Analysis
Clone this repo or download the files.

Set your working directory in DAPC_analysis.R to the folder containing the data files.

Run the R script using RStudio or R console.

## Results include:

Clustering summaries

DAPC scatter plots

Assignment plots

Nei's distance matrices

Mantel test for IBD

## References
Jombart, T., Devillard, S., & Balloux, F. (2010). Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genetics, 11(1), 94.

Kamvar, Z. N., Tabima, J. F., & Grünwald, N. J. (2014). Poppr: an R package for genetic analysis of populations with clonal, partially clonal, and/or sexual reproduction. PeerJ, 2, e281.

## License
This project is released under the MIT License.
