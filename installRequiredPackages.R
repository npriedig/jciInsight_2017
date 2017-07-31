# prerequisite packages
requiredPackages.cran <- c(
  "devtools", 
  "reshape",  
  "RColorBrewer", 
  "ggplot2", 
  "gplots",
  "corrplot",
  "lattice",
  "survminer")

requiredPackages.bioconductor <- c(
  "genefu",
  "ComplexHeatmap",
  "edgeR",
  "biomaRt",
  "DESeq2")

# install packages if not installed on your system
# CRAN
install.packages(requiredPackages.cran)
# bioconductor 
source("https://bioconductor.org/biocLite.R")
biocLite(requiredPackages.bioconductor)

# QoRTs install, version used for jciInsight_2017 v1.1.8 (https://hartleys.github.io/QoRTs/)
install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",
                   repos=NULL, 
                   type="source")
