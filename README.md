# GC-matched Randomization for Genomic Regions

[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/Elsa720/UcscGenomeRegexScanner/blob/main/LICENSE)

A custom randomization function for generating genomic controls that match target regions by **length** and **GC content**. This powerful R function is used as the randomization function part of the [`RegioneR`](https://www.bioconductor.org/packages/release/bioc/vignettes/regioneR/inst/doc/regioneR.html) packages for the Permutation test.

## Features

- **GC content matching**: Ensures the randomized control shares the same GC distribution as your target peaks/regions.
- **Support for Blacklists (Masking)**: Easily exclude gaps, repetitive elements, or specific genomic regions from the sampling pool.
- **Customizable Redundancy**: Choose whether random regions should overlap with original regions (`non.overlapping`, default=FALSE) or with each other (`random_overlap`, default=TRUE).
- **Dynamic Tolerance**: Automatically widens the GC tolerance if matching is difficult for certain regions, preventing infinite loops.

## Installation

```r
# download the gc_matched_randomize.R on the current working dir and then load the function locally
download.file("https://raw.githubusercontent.com/Elsa720/GCMatchedRandomization/main/gc_matched_randomize.R", 
              destfile = "gc_matched_randomize.R")
source("gc_matched_randomize.R")

#or load the function online

source("https://github.com/Elsa720/GCMatchedRandomization/raw/refs/heads/main/gc_matched_randomize.R")

```
Since this function depends on `GenomicRanges` and `Biostrings`, please ensure that you have these Bioconductor packages installed:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "Biostrings"))
```

After loading, you can run the following command to check if the function exists:

```r
print(gc_matched_randomize)
```

If the function's code is output, it means you have successfully deployed the function in your R environment!

**Notes:**
You can also load the function online via: 
1. using `devtools` or `remotes` via `source_url()`；
2. using the terminal to clone the function using `git clone` and then loading it via `source()`.

## Usage Example

```r
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

# Define your target regions (e.g., ChIP-seq peaks) and mask regions (e.g., blacklist regions)
query_regions <- toGRanges("path/to/my_peaks.bed")
blacklist_gr <- toGRanges("path/to/my_blacklist.bed")

# Generate GC-matched control regions
random_control <- gc_matched_randomize(
  A               = query_regions,
  genome          = genome,
  mask            = blacklist_gr,  # Optional: regions to exclude
  non.overlapping = TRUE,          # Do not overlap original A, defalut = FALSE
  random_overlap  = FALSE,         # Tolerable overlapping of random regions, default = TRUE
  per.chromosome  = FALSE,         # Sample from the same chr distrubution as A, default = FALSE
  gc.tol          = 0.05,          # Initial ±5% GC difference
  max.iter        = 1000,          # max iternation, default = 500
  verbose         = TRUE           # Whether to print progress messages, default = FALSE
)

# Use it as a custom randomization function in regioneR::permTest
# for example, the permutation test with the query_regions and another self-define region

self_defined_regions <- toGRanges("path/to/my_interested_regions.bed")

pt <- permTest(
    A=query_regions, 
    B=self_defined_regions, 
    randomize.function=random_control, 
    ntimes = 1000,
    evaluate.function = xxx,
    ...)
```

## Parameters

- `A`: Input `GRanges` object.
- `genome`: `BSgenome` object for the corresponding assembly.
- `mask`: `GRanges` object of regions to avoid (e.g., gap regions).
- `random_overlap`: If `FALSE`, it ensures that the generated set of random regions have no internal overlaps.
- `non.overlapping`: If `TRUE`, it ensures that the generated set of random regions have no overlaps with the input `A` region set.
- `per.chromosome`: If `TRUE`, it ensures that the generated set of random regions have the same chromosome distrubtion with the input `A` region set.
- `gc.tol`: Initial GC tolerance (absolute difference), default is 0.05.
- `max.iter`: Maximum iterations per region before widening GC tolerance, default is 500.


## Acknowledgment and Contribution

Thanks for the contribution of the function of the permutation test of [regioneR](https://www.bioconductor.org/packages/release/bioc/vignettes/regioneR/inst/doc/regioneR.html).

If you have any other questions or suggestions, please feel free to open an issue or contact 24111510029@m.fudan.edu.cn. We appreciate everyone's contribution!