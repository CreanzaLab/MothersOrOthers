# MothersOrOthers
## Required packages
This code runs in R (written for 4.3.1) with several packages. The following can be installed with `install.packages`: `maps, ggplot2, magrittr, dplyr`. To generate the mtDNA distance matrix, you will also need to install the `ape` package.

Instructions for installing the `rethinking` package can be found here: https://github.com/rmcelreath/rethinking
## Running the code
In the `code` directory, either run all the code together with `Rscript run_all_together.R` or source each file separately in R or Rstudio:

> "collect_data_for_analysis.R"   : uses the files in the 'data' directory to generate objects that are used for analysis. Must be run before other scripts.
>
> "worldwide_threshold_coef.R"    : analysis of worlwide trends, and creation of line plots for nuclear and mtDNA analyses.
> 
> "per_pop_coef.R"                : analysis of coefficient values per population and the creation of maps to plot those values.
> 
> "worldwide_endoexogamy_coef.R"  : extended analysis that measures the effects of endogamy -- rather than those of postmarital residence and descent patterns, which are included in the previous two files.
