# Faster-R-cor-function
Get correlation analysis done faster when data.frame is BIG by using parallelization :) 
I started working on this when I painfully experienced that the in-built R cor() function took me an entire day to complete.

# Installation
```R
devtools::install_github("utnesp/Easy-bioMart")
devtools::install_github("utnesp/Faster-R-cor-function")
library(corParallell)
```

# Example
Co-expression analysis of MYCN (ENSG00000134323) to find genes with strong correlation in a dataset containing 171 samples and 10 000 genes.

```R

> counts[1:6,1:6]
                sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
ENSG00000227232      109      117       63      179       67       81
ENSG00000268903       18        3        4       16        1       18
ENSG00000269981       14        1        5       12        1        7
ENSG00000134323      463      472      261      247      350      475
ENSG00000228463        0       14       12       65       21        8
ENSG00000237094      135       47       30       32       47       99

## the output will be saved in MYCN.cor.txt
cor.parallell(counts, "ENSG00000134323", file = "/path/to/file/MYCN.cor.txt")

```

You can also set:
```R
correlation_type = "spearman"   ## use "pearson", "kendall", or "spearman" (default "pearson")
annotate = T                    ## annotate file with gene names and biotype
read.file = T                   ## will read in file, and assign it to global environment with name MYCN.cor
no_cores = 5                    ## default uses all cores - 1
```

Good luck! :)

