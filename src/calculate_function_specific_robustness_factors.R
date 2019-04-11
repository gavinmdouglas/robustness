#!/bin/env Rscript
#
# Author: Alex Eng
# Date: 1/18/2018

###############################################################################

suppressPackageStartupMessages(library(argparse))

# Determine relative path to optional input files based on full path to script.
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# Defining arguments and parser for command line arguments
parser = ArgumentParser(description = "Calculate function-specific robustness factors.")

# Define required positional arguments
parser$add_argument("measure_table",
                    help = "Table of measures for perturbated function relative differences")

parser$add_argument("output_robustness_factor_table",
                    help = "The location to save the table of function-specific robustness coefficients")

# Define options
parser$add_argument("--file_handling",
                    default = paste(script.basename, "/file_handling.R", sep=""),
                    help = "Location of custom library for reading files")

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

# Calculate attenuation for each function first
attenuation_table = measure_table[,lapply(.SD, function(col){
  subset_data = data.table(func_shift = col, tax_shift = measure_table$taxon_weighted_unifrac, min_val = rep(.Machine$double.eps, nrow(measure_table)))
  return(coef(nls(log(exp(func_shift) + min_val) ~ -a + (b * log(exp(tax_shift) + min_val)), data = subset_data, start = c(a = 1, b = 2), control=nls.control(maxiter = 1000, warnOnly = TRUE)))[1])}), .SDcols = which(grepl("^ko[0-9]*$", colnames(measure_table)) | grepl("^Superpathway_.*$", colnames(measure_table))), by = "original_sample"]

# Next calculate buffering
buffering_table = measure_table[,lapply(.SD, function(col){
  subset_data = data.table(func_shift = col, tax_shift = measure_table$taxon_weighted_unifrac, min_val = rep(.Machine$double.eps, nrow(measure_table)))
  return(coef(nls(log(exp(func_shift) + min_val) ~ -a + (b * log(exp(tax_shift) + min_val)), data = subset_data, start = c(a = 1, b = 2), control=nls.control(maxiter = 1000, warnOnly = TRUE)))[2])
}), .SDcols = which(grepl("^ko[0-9]*$", colnames(measure_table)) | grepl("^Superpathway_.*$", colnames(measure_table))), by = "original_sample"]

# Merge robustness factor tables
melted_attenuation_table = melt(attenuation_table, id.vars = "original_sample", variable.vars = colnames(attenuation_table)[which(!(colnames(attenuation_table) %in% c("original_sample")))], variable.name = "function", value.name = "attenuation")
melted_buffering_table = melt(buffering_table, id.vars = "original_sample", variable.vars = colnames(buffering_table)[which(!(colnames(buffering_table) %in% c("original_sample")))], variable.name = "function", value.name = "buffering")
robustness_factor_table = merge(melted_attenuation_table, melted_buffering_table, by = c("original_sample", "function"))

# Write the table of function-specific robustness factors
write.table(robustness_factor_table, args$output_robustness_factor_table, sep="\t", quote=F, row.names=F, col.names=T)