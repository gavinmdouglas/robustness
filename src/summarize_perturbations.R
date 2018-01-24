#!/bin/env Rscript
#
# Author: Alex Eng
# Date: 1/18/2018

###############################################################################

suppressPackageStartupMessages(library(argparse))

# Defining arguments and parser for command line arguments
parser = ArgumentParser(description = "Generate tables of measures summarizing the original community and its relation to the perturbed community compositions.")

# Define required positional arguments
parser$add_argument("taxonomic_profile_table", help = "Taxonoimc profile table with perturbed profiles")
parser$add_argument("function_table", help = "Function table with perturbed functional capacities")
parser$add_argument("weighted_unifrac_table", help = "Table of weighted UniFrac dissimilarities between the original sample and perturbations")
parser$add_argument("output_function_distribution_feature_table", help = "The location to save the table of community-level function distribution feature values")
parser$add_argument("output_measure_table", help = "The location to save the table of community-level difference measures")
parser$add_argument("output_function_relative_difference_table", help = "The location to save the table of relative function differences")

# Define options
parser$add_argument("--fast_matrix_multiply", default = "src/fast_matrix_multiplication.cpp", help = "Location of the Rcpp implementation of matrix multiplaction (default: %(default)s)")
parser$add_argument("--file_handling", default = "src/file_handling.R", help = "Location of custom library for reading files")
parser$add_argument("--function_distribution_features", default = "src/function_distribution_features.R", help = "Location of custom library for function distribution feature calculation (default: %(default)s)")
parser$add_argument("--genome_content_table", default = "data/ko_13_5_precalculated.tab.gz", help = "The genome content table for each taxon (default: %(default)s)")
parser$add_argument("--bacterial_function_table", default = "data/KO_BACTERIAL_KEGG_2013_07_15.lst", help = "The table of bacterial functions (default: %(default)s)")
parser$add_argument("--pathway_mapping_table", default = "data/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15", help = "The table mapping functions to pathways (default: %(default)s)")
parser$add_argument("--pathway_label_table", default = "data/map_title.tab", help = "Table with names for pathways (default: %(default)s)")
parser$add_argument("--pathway_label_to_superpathway_label_mapping_table", default = "data/pathway_to_superpathway_mapping.tab", help = "Table mapping human readable pathway labels to human readable superpathway labels (default: %(default)s)")

# Parse command line arguments
args = parser$parse_args()

# Load and/or source necessary libraries
library(data.table)
library(Rcpp)
library(foreach)
library(coop)
sourceCpp(args$fast_matrix_multiply)
source(args$file_handling)
source(args$function_distribution_features)

###############################################################################

### Read tables

taxonomic_profile_table = custom_read(args$taxonomic_profile_table)
function_table = custom_read(args$function_table)
weighted_unifrac_table = custom_read(args$weighted_unifrac_table)
genome_content_table = custom_read(args$genome_content_table)
bacterial_function_table = custom_read(args$bacterial_function_table, header = F)
pathway_mapping_table = custom_read(args$pathway_mapping_table)
pathway_label_table = custom_read(args$pathway_label_table, header=F, colClasses = c("string", "string"))
pathway_label_to_superpathway_label_mapping_table = custom_read(args$pathway_label_to_superpathway_label_mapping_table)

###############################################################################

### Format column names of tables

colnames(taxonomic_profile_table)[1] = "taxon"
colnames(function_table)[1] = "function_name"
colnames(genome_content_table)[1] = "taxon"
colnames(pathway_mapping_table)[1] = "function_name"
colnames(pathway_label_table) = c("pathway", "pathway_label")
colnames(pathway_label_to_superpathway_label_mapping_table)[3:4] = c("superpathway_label", "pathway_label")

###############################################################################

### Process tables

# Convert taxon ids to character for indexing by taxon and processing taxon ids
taxonomic_profile_table[,taxon:=as.character(taxon)]
genome_content_table[,taxon:=as.character(taxon)]

# Remove taxa from the genome content table that are not present in the taxonomic profile table and then sort to match
genome_content_table = genome_content_table[which(taxon %in% taxonomic_profile_table$taxon)]
taxonomic_profile_table = taxonomic_profile_table[order(taxon)]
genome_content_table = genome_content_table[order(taxon)]

# Remove functions from the genome content table that have not been labeled bacterial
bacterial_functions = unlist(bacterial_function_table)
bacterial_function_filtered_content_cols = c(1, which(colnames(genome_content_table) %in% bacterial_functions))
genome_content_table = genome_content_table[,bacterial_function_filtered_content_cols,with=F]

# Normalize pathway contributions for each function
pathway_mapping_table[pathway_mapping_table == 0] = NA
melted_pathway_mapping_table = melt(pathway_mapping_table, id.vars = "function_name", na.rm=T, variable.name = "pathway", value.name = "contribution")
pathway_function_counts = melted_pathway_mapping_table[,.N,by="function_name"]
melted_pathway_mapping_table = merge(melted_pathway_mapping_table, pathway_function_counts, by = "function_name")
melted_pathway_mapping_table[,normalized_contribution:=contribution/N]
pathway_mapping_table = dcast(melted_pathway_mapping_table, function_name ~ pathway, value.var="normalized_contribution")
pathway_mapping_table[is.na(pathway_mapping_table)] = 0

# Generate table mapping pathways in koXXXXX format to tagged superpathway names
colnames(pathway_label_table) = c("pathway", "pathway_label")
pathway_label_table[,pathway := paste("ko", pathway, sep="")]
pathway_label_to_superpathway_label_mapping_table = pathway_label_to_superpathway_label_mapping_table[,3:4,with=F]
pathway_label_to_superpathway_label_mapping_table = pathway_label_to_superpathway_label_mapping_table[!duplicated(pathway_label_to_superpathway_label_mapping_table)]
pathway_to_superpathway_mapping_table = merge(pathway_label_table, pathway_label_to_superpathway_label_mapping_table, by="pathway_label")
pathway_to_superpathway_mapping_table[,superpathway := paste("Superpathway_", superpathway_label, sep="")]

# Generate matrix version of genome content table for matrix multiplication
genome_content_mat = matrix(as.numeric(as.matrix(genome_content_table[,2:ncol(genome_content_table),with=F])), nrow = nrow(genome_content_table))
colnames(genome_content_mat) = colnames(genome_content_table)[2:ncol(genome_content_table)]
rownames(genome_content_mat) = as.character(genome_content_table$taxon)

###############################################################################

### Generate table of taxon pathway content

# Create pathway mapping matrix for matrix multiplication with only functions that exist in the genome content table
existing_functions = which(pathway_mapping_table$function_name %in% colnames(genome_content_table))
pathway_mapping_mat = as.matrix(pathway_mapping_table[existing_functions,2:ncol(pathway_mapping_table),with=F])
rownames(pathway_mapping_mat) = as.character(pathway_mapping_table[existing_functions]$function_name)

# Create copy of genome content matrix to match the pathway mapping matrix
temp_genome_content_mat = genome_content_mat[,which(colnames(genome_content_mat) %in% rownames(pathway_mapping_mat))]

# Order rows/columns to match between matrices
pathway_mapping_mat_row_order = sapply(colnames(temp_genome_content_mat), function(function_name){
  return(which(rownames(pathway_mapping_mat) == function_name))
})
pathway_mapping_mat = pathway_mapping_mat[pathway_mapping_mat_row_order,]

# Multiply the matrices together to get a taxon-to-pathway mapping table
pathway_content_table = as.data.table(fast_matrix_multiply(temp_genome_content_mat, pathway_mapping_mat))

# Add correct row and column labels
colnames(pathway_content_table) = colnames(pathway_mapping_table)[2:ncol(pathway_mapping_table)]
pathway_content_table$taxon = genome_content_table$taxon
pathway_content_table = pathway_content_table[,c(ncol(pathway_content_table), 1:(ncol(pathway_content_table) - 1)),with=F]

# Remove columns for pathways with no contributions from the present taxa
nonzero_cols = which(colSums(pathway_content_table[,2:ncol(pathway_content_table)]) > 0)
pathway_content_table = pathway_content_table[,c(1,nonzero_cols + 1),with=F]

###############################################################################

### Generate table of taxon genome presence

genome_presence_table = genome_content_table
for (column_name in colnames(genome_content_table)[2:ncol(genome_content_table)]){
  genome_presence_table[[column_name]] = ifelse(genome_content_table[[column_name]] > 0, 1, 0)
}

###############################################################################

### Generate table of taxon pathway presence

pathway_presence_table = pathway_content_table
for (column_name in colnames(pathway_content_table)[2:ncol(pathway_content_table)]){
  pathway_presence_table[[column_name]] = ifelse(pathway_content_table[[column_name]] > 0, 1, 0)
}

###############################################################################

### Generate pathway-level functional profile table for original and perturbed community compositions

# Remove functions that can't be mapped to pathways and convert to a matrix for matrix multiplication
existing_functions = which(function_table$function_name %in% pathway_mapping_table$function_name)
function_mat = as.matrix(function_table[existing_functions,2:ncol(function_table),with=F])
rownames(function_mat) = as.character(function_table[existing_functions]$function_name)

# Create a pathway mapping matrix with functions matched to the function table
existing_functions = which(pathway_mapping_table$function_name %in% function_table$function_name)
pathway_mapping_mat = as.matrix(pathway_mapping_table[existing_functions,2:ncol(pathway_mapping_table),with=F])
rownames(pathway_mapping_mat) = as.character(pathway_mapping_table[existing_functions]$function_name)

# Reorder functions in pathway mapping matrix to match up with functions in function matrix
pathway_mapping_mat_row_order = sapply(rownames(function_mat), function(function_name){
  return(which(rownames(pathway_mapping_mat) == function_name))
})
pathway_mapping_mat = pathway_mapping_mat[pathway_mapping_mat_row_order,]

# Multiply matrices together to get community pathway-level profiles
pathway_table = as.data.table(fast_matrix_multiply(t(pathway_mapping_mat), function_mat))

# Add correct row and column labels
colnames(pathway_table) = colnames(function_table)[2:ncol(function_table)]
pathway_table$pathway = colnames(pathway_mapping_mat)
pathway_table = pathway_table[,c(ncol(pathway_table), 1:(ncol(pathway_table) - 1)),with=F]

# Remove rows for pathways with no abundance
nonzero_rows = which(rowSums(pathway_table[,2:ncol(pathway_table),with=F]) > 0)
pathway_table = pathway_table[nonzero_rows]

###############################################################################

### Generate a superpathway-level profile table

temp_pathway_table = pathway_table[which(pathway %in% pathway_to_superpathway_mapping_table$pathway)]
temp_pathway_table = merge(temp_pathway_table, pathway_to_superpathway_mapping_table, by = "pathway")
superpathway_table = temp_pathway_table[,lapply(.SD, sum),.SDcols = which(!(colnames(temp_pathway_table) %in% c("pathway", "pathway_label", "superpathway", "superpathway_label"))),by="superpathway"]

###############################################################################

### Create function and pathway profile tables for just original community

original_function_table = function_table[,c(1,2),with=F]
original_pathway_table = pathway_table[,c(1,2),with=F]

###############################################################################

### Calculate function distribution metrics

# Get results for each set of function distribution metrics
unique_pathway_abundance_features = get_unique_content_features(taxonomic_profile_table[,c(1,2),with=F], original_pathway_table, pathway_presence_table)
pathway_redundancy_features = get_entropy_redundancy_features(taxonomic_profile_table[,c(1,2),with=F], original_pathway_table, pathway_content_table)
taxon_pathway_similarity_features = get_taxon_similarity_features(taxonomic_profile_table[,c(1,2),with=F], pathway_content_table)
genome_size_features = get_genome_size_features(taxonomic_profile_table[,c(1,2),with=F], pathway_presence_table, pathway_content_table)
unique_content_taxon_abundance_features = get_number_and_abundance_taxons_with_unique_content(taxonomic_profile_table[,c(1,2),with=F], pathway_presence_table, pathway_content_table)
pathway_containing_taxon_features = get_num_and_abundance_taxons_containing_content(taxonomic_profile_table[,c(1,2),with=F], pathway_presence_table, pathway_content_table)
pathway_occurrence_similarity_features = get_function_occurrence_similarity(taxonomic_profile_table[,c(1,2),with=F], pathway_content_table)

# Merge results together
sample_col_name = colnames(unique_pathway_abundance_features)[1]
function_distribution_feature_table = merge(unique_pathway_abundance_features, pathway_redundancy_features, by=sample_col_name)
function_distribution_feature_table = merge(function_distribution_feature_table, taxon_pathway_similarity_features, by=sample_col_name)
function_distribution_feature_table = merge(function_distribution_feature_table, genome_size_features, by=sample_col_name)
function_distribution_feature_table = merge(function_distribution_feature_table, unique_content_taxon_abundance_features, by=sample_col_name)
function_distribution_feature_table = merge(function_distribution_feature_table, pathway_containing_taxon_features, by=sample_col_name)
function_distribution_feature_table = merge(function_distribution_feature_table, pathway_occurrence_similarity_features, by=sample_col_name)

# Write the table of function distribution features
write.table(function_distribution_feature_table, args$output_function_distribution_feature_table, quote=F, sep="\t", row.names=F, col.names=T)

###############################################################################

### Calculate functional dissimilarities between the original and perturbed compositions

# Generate table of pathway profile cosine dissimilarities between the original sample and the perturbed samples
pathway_table_mat = as.matrix(pathway_table[,2:ncol(pathway_table),with=F])
cosine_vec = 1 - as.numeric(sapply(2:ncol(pathway_table_mat), function(col){return(cosine(pathway_table_mat[,1], pathway_table_mat[,col]))}))
pathway_cosine_dist_table = data.table(original_sample = rep(colnames(pathway_table_mat)[1], ncol(pathway_table_mat) - 1), perturbation = colnames(pathway_table_mat)[2:ncol(pathway_table_mat)], pathway_cosine_dissimilarity = log(cosine_vec))

# Combine with the weighted unifrac dissimilarities
measure_table = merge(weighted_unifrac_table, pathway_cosine_dist_table, by = c("original_sample", "perturbation"))

# Write the table of taxonomic and functional dissimilarities between original and perturbed compositions
write.table(measure_table, args$output_measure_table, quote = F, sep = "\t", row.names=F, col.names=T)

###############################################################################

### Calculating individual pathway and superpathway shifts

# Convert pathway and superpathway tables to matrices for matrix math
pathway_mat = as.matrix(pathway_table[,2:ncol(pathway_table),with=F])
rownames(pathway_mat) = pathway_table$pathway
superpathway_mat = as.matrix(superpathway_table[,2:ncol(superpathway_table),with=F])
rownames(superpathway_mat) = superpathway_table$superpathway

# Calculate relative difference in pathway abundances and convert to log scale for saving with increased precision
pathway_relative_diff_mat = log(abs(pathway_mat[,2:ncol(pathway_mat)] - pathway_mat[,1])/pathway_mat[,1])

# Transpose and convert to a datatable
pathway_relative_diff_table = as.data.table(t(pathway_relative_diff_mat))

# Add in labels for original community and perturbations
pathway_relative_diff_table[,original_sample:=rep(colnames(pathway_mat)[1], nrow(pathway_relative_diff_table))]
pathway_relative_diff_table[,perturbation:=colnames(pathway_relative_diff_mat)]

# Calculate relative difference in pathway abundances and convert to log scale for saving with increased precision
superpathway_relative_diff_mat = log(abs(superpathway_mat[,2:ncol(superpathway_mat)] - superpathway_mat[,1])/superpathway_mat[,1])

# Transpose and convert to a datatable
superpathway_relative_diff_table = data.table(t(superpathway_relative_diff_mat))

# Add in labels for original community and perturbations
superpathway_relative_diff_table[,original_sample:=rep(colnames(superpathway_mat)[1], nrow(superpathway_relative_diff_table))]
superpathway_relative_diff_table[,perturbation:=colnames(superpathway_relative_diff_mat)]

# Merge pathway and superpathway relative difference tables
function_relative_diff_table = merge(pathway_relative_diff_table, superpathway_relative_diff_table, by = c("original_sample", "perturbation"))

# Write the table of individual function relative shifts in abundance
write.table(function_relative_diff_table, args$output_function_relative_difference_table, quote=F, sep="\t", row.names=F, col.names=T)