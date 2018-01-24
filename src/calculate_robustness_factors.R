#!/bin/env Rscript
#
# Author: Alex Eng
# Date: 1/18/2018

###############################################################################

suppressPackageStartupMessages(library(argparse))

# Defining arguments and parser for command line arguments
parser = ArgumentParser(description = "Calculate community robustness factors.")

# Define required positional arguments
parser$add_argument("measure_table", help = "Table of measures for perturbation dissimilarities")
parser$add_argument("output_robustness_factor_table", help = "The location to save the table of community robustness coefficients")

# Define options
parser$add_argument("--file_handling", default = "src/file_handling.R", help = "Location of custom library for reading files")

# Parse command line arguments
args = parser$parse_args()

# Load and/or source necessary libraries
library(data.table)
library(methods)
library(stringr)
source(args$file_handling)

###############################################################################

### Read table

measure_table = custom_read(args$measure_table)

###############################################################################

### Calculating robustness factors

# Supply the minimum possible float to add to all dissimilarities since we are fitting in log space and we can't take the log of zero
measure_table$min_val = rep(.Machine$double.eps, nrow(measure_table))

# Calculate attenuation first
attenuation_table = measure_table[,coef(nls(log(exp(pathway_cosine_dissimilarity) + min_val) ~ -a + (b * log(exp(taxon_weighted_unifrac) + min_val)), data = .SD, start = c(a = 1, b = 2), control=nls.control(maxiter = 1000, warnOnly = TRUE)))[1], by = "original_sample"]
colnames(attenuation_table)[which(colnames(attenuation_table) == "V1")] = "attenuation"

# Next calculate buffering
buffering_table = measure_table[,coef(nls(log(exp(pathway_cosine_dissimilarity) + min_val) ~ -a + (b * log(exp(taxon_weighted_unifrac) + min_val)), data = .SD, start = c(a = 1, b = 2), control=nls.control(maxiter = 1000, warnOnly = TRUE)))[2], by = "original_sample"]
colnames(buffering_table)[which(colnames(buffering_table) == "V1")] = "buffering"

# Merge attenuation and buffering tables
robustness_factor_table = merge(attenuation_table, buffering_table, by = "original_sample")

# Write the table of robustness factors
write.table(robustness_factor_table, args$output_robustness_factor_table, sep="\t", quote=F, row.names=F, col.names=T)