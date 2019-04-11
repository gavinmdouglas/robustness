#!/bin/env Rscript
#
# Author: Alex Eng
# Date: 1/18/2018
# Edited by Gavin Douglas (4/11/2019).

###############################################################################

suppressPackageStartupMessages(library(argparse))

# Determine relative path to optional input files based on full path to script.
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# Defining arguments and parser for command line arguments
parser = ArgumentParser(description = "Generate perturbed taxonmic and functional profile tables for a given community.")

# Define required positional arguments
parser$add_argument("taxonomic_profile_table", help = "Table defining original taxonomic profile to perturb")
parser$add_argument("output_taxonomic_profile_table", help = "The location to save the table of taxonoimc profile perturbations")
parser$add_argument("output_function_table", help = "The location to save the function perturbation table")
parser$add_argument("output_weighted_unifrac_table", help = "The location to save the weighted UniFrac dissimilarities between the original community and the perturbed community compositions")

# Define options
parser$add_argument("--fast_matrix_multiply",
                    default = paste(script.basename, "/fast_matrix_multiplication.cpp", sep=""),
                    help = "Location of the Rcpp implementation of matrix multiplaction (default: %(default)s)")

parser$add_argument("--file_handling",
                    default = paste(script.basename, "/file_handling.R", sep=""),
                    help = "Location of custom library for reading files")

parser$add_argument("--perturbation_list",
                    type = "double",
                    default = seq(1.2, 10, 0.2),
                    nargs = "+",
                    help = "The list of maximum perturbation magnitudes to generate perturbations for (default: %(default)s)")

parser$add_argument("--num_replicates",
                    type = "integer",
                    default = 100,
                    help = "The number of perturbations to generate for each maximum perturbation magnitude (default: %(default)s)")

parser$add_argument("--copy_number_norm_table",
                    default = paste(script.basename, "/../data/16S_13_5_precalculated.tab.gz", sep=""),
                    help = "The table of 16S copy numbers for each taxon to normalize taxon counts with (default: %(default)s)")

parser$add_argument("--genome_content_table",
                    default = paste(script.basename, "/../data/ko_13_5_precalculated.tab.gz", sep=""),
                    help = "The genome content table for each taxon (default: %(default)s)")

parser$add_argument("--bacterial_function_table",
                    default = paste(script.basename, "/../data/KO_BACTERIAL_KEGG_2013_07_15.lst", sep=""),
                    help = "The table of bacterial functions (default: %(default)s)")

parser$add_argument("--tree",
                    default = paste(script.basename, "/../data/gg_13_5_otus_99_annotated.tree", sep=""),
                    help = "The tree to use for calculating weighted UniFrac (default: %(default)s)")

parser$add_argument("--normal_scale",
                    action = "store_true",
                    help = "Flag to indicate that filtering should be done on normal scale (default filtering is done on log scale)")

parser$add_argument("--range_max",
                    type = "double",
                    default = 0.4,
                    help = "Set the maximum of the range of taxonomic profile dissimilarities, any perturbations with a greater taxonomic profile dissimilarity will be filtered out (default: %(default)s)")

parser$add_argument("--range_min",
                    type = "double",
                    default = 0.0001,
                    help = "Set the minimum of the range of taxonomic profile dissimilarities, any perturbations with a lower taxonomic profile dissimilarity will be filtered out (default: %(default)s)")

parser$add_argument("--num_windows",
                    type = "integer",
                    default = 50,
                    help = "The number of windows to use when filtering (default: %(default)s)")

parser$add_argument("--num_perturbations",
                    type = "integer",
                    default = 50,
                    help = "The number of perturbations to subsample to within each window (default: %(default)s)")

parser$add_argument("--drop_low_perturbation_windows",
                    action = "store_true",
                    help = "Flag to indicate that perturbations in windows with fewer total perturbations than num_perturbations should be ignored")

# Parse command line arguments
args = parser$parse_args()

# Load and/or source necessary libraries
library(data.table)
library(Rcpp)
library(foreach)
library(phyloseq)
library(ape)
sourceCpp(args$fast_matrix_multiply)
source(args$file_handling)

###############################################################################

### Read tables

taxonomic_profile_table = custom_read(args$taxonomic_profile_table)
norm_16S_table = custom_read(args$copy_number_norm_table)
genome_content_table = custom_read(args$genome_content_table, fill = T)
bacterial_function_table = custom_read(args$bacterial_function_table, header = F)

###############################################################################

### Format column names of tables

colnames(taxonomic_profile_table)[1] = c("taxon")
colnames(norm_16S_table) = c("taxon", "copy_number")
colnames(genome_content_table)[1] = "taxon"

###############################################################################

### Process tables

# Filter for taxa that are present
nonzero_rows = which(unlist(taxonomic_profile_table[,2,with=F]) > 0)
taxonomic_profile_table = taxonomic_profile_table[nonzero_rows]
num_taxa = nrow(taxonomic_profile_table)

# Save the number of reads to keep read counts constant between original profile and perturbations
num_reads = sum(unlist(taxonomic_profile_table[,2,with=F]))

# Remove functions from the genome content table that have been labeled on non-bacterial
bacterial_functions = unlist(bacterial_function_table)
bacterial_function_filtered_content_cols = c(1, which(colnames(genome_content_table) %in% bacterial_functions))
genome_content_table = genome_content_table[,bacterial_function_filtered_content_cols,with=F]

# Remove rows without NSTI data from the genome content table (either KO descriptions or blank rows)
genome_content_table = na.omit(genome_content_table, ncol(genome_content_table))

###############################################################################

### Correct taxonomic abundances by 16S copy number

# Merge normalization factor table with taxonomic profile table to match normalization factors to taxa
taxonomic_profile_table = merge(taxonomic_profile_table, norm_16S_table, by="taxon")

# Calculated corrected abundances and then readjust to same read count as before
taxonomic_profile_table[,corrected_abundances:=get(colnames(taxonomic_profile_table)[2]) / copy_number]
taxonomic_profile_table[,(colnames(taxonomic_profile_table)[2]):=(corrected_abundances / sum(corrected_abundances)) * num_reads]
taxonomic_profile_table[,copy_number:=NULL]
taxonomic_profile_table[,corrected_abundances:=NULL]

# Convert taxon ids to character for indexing by taxon
taxonomic_profile_table[,taxon:=as.character(taxon)]
genome_content_table[,taxon:=as.character(taxon)]

###############################################################################

### Generate perturbed compositions

invisible(foreach(new_col = 1:(length(args$perturbation_list) * args$num_replicates)) %do% {
  
  # Use the same set of 0-indexed replicate IDs for each perturbation size in the perturbation list
  replicate_num = (new_col - 1) %% args$num_replicates
  
  # Grab the perturbation size for the current perturbation
  perturbation_sd = args$perturbation_list[ceiling(new_col / args$num_replicates)]
  
  # Generate the vector of perturbation coefficients
  perturbation_vec = runif(num_taxa, min = 1, max = perturbation_sd)
  
  # For each taxon, randomly select whether its abundance will increased or decreased
  direction_vec = rbinom(num_taxa, 1, 0.5)
  
  # If the direction value is 1, then the taxon's abundance will be multiplied by the perturbation coefficient, otherwise it will be divided by the perturbation coefficient
  fixed_perturbation_vec = ifelse(direction_vec == 1, perturbation_vec, 1/perturbation_vec)
  perturbed_vec = unlist(taxonomic_profile_table[,2,with=F]) * fixed_perturbation_vec
  
  # Renormalize the perturbed abundance vector
  normalized_perturbed_vec = (perturbed_vec / sum(perturbed_vec)) * num_reads
  
  # Add the perturbed abundance vector to the taxonomic profile table
  taxonomic_profile_table[,(paste(colnames(taxonomic_profile_table)[2], "_perturbation_", replicate_num, "_sd_", perturbation_sd, sep="")):=normalized_perturbed_vec]
})

###############################################################################

### Filter perturbations to create more uniform distribution over range of perturbation sizes

# Read and filter tree for calculating weighted UniFrac distances
tree = read_tree_greengenes(args$tree)
tree = drop.tip(tree, tree$tip.label[which(!(tree$tip.label %in% taxonomic_profile_table$taxon))])

# Convert taxonomic profile table to a matrix for conversion to phyloseq object
taxonomic_profile_table_mat = as.matrix(taxonomic_profile_table[,2:ncol(taxonomic_profile_table),with=F])
rownames(taxonomic_profile_table_mat) = as.character(taxonomic_profile_table$taxon)

# Generate table of weighted UniFrac distances between the original sample and the perturbed samples
weighted_unifrac_vec = unlist(foreach(col = 2:ncol(taxonomic_profile_table_mat)) %do% {
  return(UniFrac(phyloseq(otu_table(taxonomic_profile_table_mat[,c(1,col)], taxa_are_rows=T), tree), weighted=T))
})

# Remove perturbations above the maximum range or below the minimum range desired
taxonomic_profile_table = taxonomic_profile_table[,c(1,2,which(weighted_unifrac_vec <= args$range_max & weighted_unifrac_vec >= args$range_min) + 2),with=F]
weighted_unifrac_vec = weighted_unifrac_vec[weighted_unifrac_vec <= args$range_max & weighted_unifrac_vec >= args$range_min]

# If not operating on normal scale, convert weighted UniFrac values and range to log scale
if (!args$normal_scale){
  weighted_unifrac_vec = log(weighted_unifrac_vec)
  args$range_min = log(args$range_min)
  args$range_max = log(args$range_max)
}

# Calculate window size and windows
sliding_window_size = (args$range_max - args$range_min) / args$num_windows
windows = seq(args$range_min, args$range_max, sliding_window_size)[1:args$num_windows]

# Subsample each window
for (window_start in windows){
  
  # Determine which sperturbations fall within the current window
  window_profile = window_start <= weighted_unifrac_vec & weighted_unifrac_vec <= window_start + sliding_window_size
  perturbations_inside_window = which(window_profile)
  perturbations_outside_window = which(!window_profile)
  
  # Determine which perturbations to keep
  perturbations_inside_window_to_keep = perturbations_inside_window
  if (length(perturbations_inside_window) >= args$num_perturbations){
    perturbations_inside_window_to_keep = sample(perturbations_inside_window, args$num_perturbations)
  }
  
  # Remove samples not kept from window
  taxonomic_profile_table = taxonomic_profile_table[,c(1,2,c(perturbations_outside_window, perturbations_inside_window_to_keep) + 2),with=F]
  weighted_unifrac_vec = weighted_unifrac_vec[c(perturbations_outside_window, perturbations_inside_window_to_keep)]
}

weighted_unifrac_table = data.table(original_sample = rep(colnames(taxonomic_profile_table)[2], ncol(taxonomic_profile_table) - 2), perturbation = colnames(taxonomic_profile_table)[3:ncol(taxonomic_profile_table)])

# If working on normal scale, convert weighted UniFrac values to log scale to increase precision
if (args$normal_scale){
  weighted_unifrac_vec = log(weighted_unifrac_vec)
}

weighted_unifrac_table$taxon_weighted_unifrac = weighted_unifrac_vec

###############################################################################

### Generate the functional profile table associated with the table of original and perturbed taxonomic profiles

# Remove taxa from the genome content table that are not present in the taxonomic abundnace table and then sort to match
genome_content_table = genome_content_table[which(taxon %in% taxonomic_profile_table$taxon)]
taxonomic_profile_table = taxonomic_profile_table[order(taxon)]
genome_content_table = genome_content_table[order(taxon)]

# Convert taxonomic profile table and function table to matrices for multiplication
taxonomic_profile_table_mat = as.matrix(taxonomic_profile_table[,2:ncol(taxonomic_profile_table),with=F])
rownames(taxonomic_profile_table_mat) = as.character(taxonomic_profile_table$taxon)
genome_content_table_mat = matrix(as.numeric(as.matrix(genome_content_table[,2:ncol(genome_content_table),with=F])), nrow = nrow(genome_content_table))
colnames(genome_content_table_mat) = colnames(genome_content_table)[2:ncol(genome_content_table)]
rownames(genome_content_table_mat) = as.character(genome_content_table$taxon)

# Multiply the matrices together to get functional profiles
function_table_mat = fast_matrix_multiply(t(genome_content_table_mat), taxonomic_profile_table_mat)

# Convert back into a datatable
function_table = as.data.table(function_table_mat)
colnames(function_table) = colnames(taxonomic_profile_table)[2:ncol(taxonomic_profile_table)]

# Add back in column of function names and reorder to put at front of table
function_table[["function"]] = colnames(genome_content_table_mat)
function_table = function_table[,c(ncol(function_table), 1:(ncol(function_table) - 1)),with=F]

# Remove functions with zero abundance
nonzero_rows = which(rowSums(function_table[,2:ncol(function_table),with=F]) > 0)
function_table = function_table[nonzero_rows]

###############################################################################

### Write resulting tables

write.table(taxonomic_profile_table, args$output_taxonomic_profile_table, quote=F, sep="\t", row.names=F, col.names=T)
write.table(function_table, args$output_function_table, quote=F, sep="\t", row.names=F, col.names=T)
write.table(weighted_unifrac_table, args$output_weighted_unifrac_table, quote=F, sep="\t", row.names=F, col.names=T)
