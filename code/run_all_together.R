library(rethinking)
library(ggplot2)
library(maps)
library(magrittr)
library(dplyr)


# Generate two lists, nuc_list and mt_data that contain geographic, linguistic,
#  genetic and cultural data within and between populations. 
source("collect_data_for_analysis.R")

# Create line plots that describe the worldwide relationship between genetics 
#  (nuclear and mtDNA), geography, and language (and, for mtDNA, ethnographic
#  traits)
source("worldwide_threshold_coef.R")

# Create maps that describe the localized relationships between genetics (
#  nuclear and mtDNA), geography and language at multiple geographic thresholds.
source("per_pop_coef.R")

# Similar to "worldwide_threshold_coef.R", but only for endogamy/exogamy, which
#  had different distributions than matriliny/matrilocality used elsewhere.
source("worldwide_endoexogamy_coef.R")