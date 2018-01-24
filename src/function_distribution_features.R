#!/bin/env Rscript

library(vegan)
library(foreach)
library(data.table)
library(coop)

write("Defining function distribution feature methods", stderr())

cosine_similarity = function(vec1, vec2){
  num = sum(vec1 * vec2)
  denom = sqrt(sum(vec1^2)) * sqrt(sum(vec2^2))
  value = num/denom
  return(value)
}

get_unique_content_features = function(taxon_table, function_table, presence_table){
  unique_cols = which(colSums(presence_table[,2:ncol(presence_table),with=F]) == 1)
  if (length(unique_cols) > 0){
    total_function_abundance = sum(function_table[,2,with=F])
    unique_function_rows = which(unlist(function_table[,1,with=F]) %in% colnames(presence_table)[unique_cols + 1])
    unique_function_abundance = sum(function_table[unique_function_rows, 2, with=F])
    return(data.table(sample = colnames(taxon_table)[2], unique_content_relative_abundance = unique_function_abundance/ total_function_abundance))
  } else {
    return(data.table(sample = colnames(taxon_table)[2], unique_content_relative_abundance = 0))
  }
}

get_entropy_redundancy_features = function(taxon_table, function_table, content_table){
  
  content_table_present_func_cols = which(colnames(content_table)[2:ncol(content_table)] %in% unlist(function_table[,1,with=F]))
  subset_present_func_content_table = content_table[,c(1, content_table_present_func_cols + 1),with=F]
  if (!all(unlist(function_table[,1,with=F]) %in% colnames(subset_present_func_content_table)[2:ncol(subset_present_func_content_table)])){
    stop(paste("Functions in function table not in content table: ", paste(unlist(function_table[,1,with=F])[which(!(unlist(function_table[,1,with=F]) %in% colnames(subset_present_func_content_table)[2:ncol(subset_present_func_content_table)]))], collapse=","), sep=""))
  }
  subset_present_func_content_table_col_order = c(1, order(colnames(subset_present_func_content_table)[2:ncol(subset_present_func_content_table)]) + 1)
  ordered_subset_present_func_content_table = subset_present_func_content_table[,subset_present_func_content_table_col_order,with=F]
  function_table_row_order = order(unlist(function_table[,1,with=F]))
  ordered_function_table = function_table[function_table_row_order]
  
  func_entropies = unlist(foreach(function_col = 2:ncol(ordered_subset_present_func_content_table)) %do% {
    subset_cols_content_table = content_table[,c(1, function_col),with=F]
    subset_cols_content_table_nonzero_rows = which(unlist(subset_cols_content_table[,2,with=F]) > 0)
    subset_cols_nonzero_content_table = subset_cols_content_table[subset_cols_content_table_nonzero_rows]
    subset_taxon_table_rows = which(unlist(taxon_table[,1,with=F]) %in% unlist(subset_cols_nonzero_content_table[,1,with=F]))
    subset_taxon_table = taxon_table[subset_taxon_table_rows]
    subset_cols_nonzero_content_table_present_rows = which(unlist(subset_cols_nonzero_content_table[,1,with=F]) %in% unlist(subset_taxon_table[,1,with=F]))
    subset_cols_present_nonzero_content_table = subset_cols_nonzero_content_table[subset_cols_nonzero_content_table_present_rows]
    subset_taxon_table_row_order = order(unlist(subset_taxon_table[,1,with=F]))
    ordered_subset_taxon_table = subset_taxon_table[subset_taxon_table_row_order]
    subset_cols_present_nonzero_content_table_row_order = order(unlist(subset_cols_present_nonzero_content_table[,1,with=F]))
    ordered_subset_cols_present_nonzero_content_table = subset_cols_present_nonzero_content_table[subset_cols_present_nonzero_content_table_row_order]
    func_abundance_vec = unlist(ordered_subset_taxon_table[,2,with=F]) * unlist(ordered_subset_cols_present_nonzero_content_table[,2,with=F])
    normalized_func_abundance_vec = func_abundance_vec/sum(func_abundance_vec)
    return(diversity(normalized_func_abundance_vec, index="shannon"))
  })
  
  func_abundances = unlist(ordered_function_table[,2,with=F])
  normalized_func_abundances = func_abundances/sum(func_abundances)
  weighted_func_entropies = func_entropies * normalized_func_abundances
  weighted_func_entropy_mean = sum(weighted_func_entropies)
  weighted_func_entropy_sd = sqrt(sum((normalized_func_abundances ^ 2) * ((func_entropies - weighted_func_entropy_mean) ^ 2)))
  weighted_func_entropy_cov = weighted_func_entropy_sd/weighted_func_entropy_mean
  return(data.table(sample = colnames(taxon_table)[2], mean_entropy_redundancy = weighted_func_entropy_mean, cov_entropy_redundancy = weighted_func_entropy_cov))
}

get_taxon_similarity_features = function(taxon_table, content_table){
  content_table_mat = as.matrix(content_table[,2:ncol(content_table),with=F])
  
  taxon_cosine_dist_mat = as.dist(cosine(t(content_table_mat)))
  taxon_euclidean_dist_mat = vegdist(content_table_mat, method = "euclidean")
  taxon_jaccard_dist_mat = vegdist(content_table_mat, method = "jaccard")
  taxon_bray_dist_mat = vegdist(content_table_mat, method = "bray")
  taxon_manhattan_dist_mat = vegdist(content_table_mat, method = "manhattan")
  taxon_pearson_correlation_dist_mat = as.dist(cor(content_table_mat, method = "pearson"))
  taxon_spearman_correlation_dist_mat = as.dist(cor(content_table_mat, method = "spearman"))

  mean_taxon_cosine_dist = mean(taxon_cosine_dist_mat, na.rm=T)
  cov_taxon_cosine_dist = sd(taxon_cosine_dist_mat, na.rm=T)/mean_taxon_cosine_dist
  mean_taxon_euclidean_dist = mean(taxon_euclidean_dist_mat, na.rm=T)
  cov_taxon_euclidean_dist = sd(taxon_euclidean_dist_mat, na.rm=T)/mean_taxon_euclidean_dist
  mean_taxon_jaccard_dist = mean(taxon_jaccard_dist_mat, na.rm=T)
  cov_taxon_jaccard_dist = sd(taxon_jaccard_dist_mat, na.rm=T)/mean_taxon_jaccard_dist
  mean_taxon_bray_dist = mean(taxon_bray_dist_mat, na.rm=T)
  cov_taxon_bray_dist = sd(taxon_bray_dist_mat, na.rm=T)/mean_taxon_bray_dist
  mean_taxon_manhattan_dist = mean(taxon_manhattan_dist_mat, na.rm=T)
  cov_taxon_manhattan_dist = sd(taxon_manhattan_dist_mat, na.rm=T)/mean_taxon_manhattan_dist
  mean_taxon_pearson_correlation_dist = mean(taxon_pearson_correlation_dist_mat, na.rm=T)
  cov_taxon_pearson_correlation_dist = sd(taxon_pearson_correlation_dist_mat, na.rm=T)/mean_taxon_pearson_correlation_dist
  mean_taxon_spearman_correlation_dist = mean(taxon_spearman_correlation_dist_mat, na.rm=T)
  cov_taxon_spearman_correlation_dist = sd(taxon_spearman_correlation_dist_mat, na.rm=T)/mean_taxon_spearman_correlation_dist
  
  
  return(data.table(sample = colnames(taxon_table)[2], mean_taxon_cosine_similarity = mean_taxon_cosine_dist, cov_taxon_cosine_similarity = cov_taxon_cosine_dist, mean_taxon_euclidean_dist = mean_taxon_euclidean_dist, cov_taxon_euclidean_dist = cov_taxon_euclidean_dist, mean_taxon_jaccard_dist = mean_taxon_jaccard_dist, cov_taxon_jaccard_dist = cov_taxon_jaccard_dist, mean_taxon_bray_dist = mean_taxon_bray_dist, cov_taxon_bray_dist = cov_taxon_bray_dist, mean_taxon_manhattan_dist = mean_taxon_manhattan_dist, cov_taxon_manhattan_dist = cov_taxon_manhattan_dist, mean_taxon_pearson_correlation_dist = mean_taxon_pearson_correlation_dist, cov_taxon_pearson_correlation_dist = cov_taxon_pearson_correlation_dist, mean_taxon_spearman_correlation_dist = mean_taxon_spearman_correlation_dist, cov_taxon_spearman_correlation_dist = cov_taxon_spearman_correlation_dist))
}

get_genome_size_features = function(taxon_table, presence_table, content_table){
  presence_table_row_sums = rowSums(presence_table[,2:ncol(content_table), with=F])
  mean_number_functions = mean(presence_table_row_sums, na.rm=T)
  sd_number_functions = sd(presence_table_row_sums, na.rm=T)
  cov_number_functions = sd_number_functions/mean_number_functions
  
  content_table_row_sums = rowSums(content_table[,2:ncol(content_table),with=F])
  mean_amount_content = mean(content_table_row_sums, na.rm=T)
  sd_amount_content = sd(content_table_row_sums, na.rm=T)
  cov_amount_content = sd_amount_content/mean_amount_content
  
  return(data.table(sample = colnames(taxon_table)[2], mean_number_functions_per_taxon = mean_number_functions, cov_number_functions_per_taxon = cov_number_functions, mean_copy_number_content_per_taxon = mean_amount_content, cov_copy_number_content_per_taxon = cov_amount_content))
}

get_number_and_abundance_taxons_with_unique_content = function(taxon_table, presence_table, content_table){
  
  unique_cols = which(colSums(presence_table[,2:ncol(presence_table),with=F]) == 1)
  if (length(unique_cols) > 0){
    unique_content_table = content_table[,c(1,unique_cols + 1),with=F]
    present_unique_rows = which(rowSums(unique_content_table[,2:ncol(unique_content_table),with=F]) > 0)
    unique_content_table = unique_content_table[present_unique_rows]
    num_taxons_with_unique_content = nrow(unique_content_table)
    taxons_with_unique_rows = which(taxon_table[[colnames(taxon_table)[1]]] %in% unique_content_table[[colnames(taxon_table)[1]]])
    abundance_taxons_with_unique_content = sum(taxon_table[taxons_with_unique_rows,2,with=F])
    
    return(data.table(sample = colnames(taxon_table)[2], num_taxons_with_unique_content = num_taxons_with_unique_content, abundance_taxons_with_unique_content = abundance_taxons_with_unique_content))
  } else {
    return(data.table(sample = colnames(taxon_table)[2], num_taxons_with_unique_content = 0, abundance_taxons_with_unique_content = 0))
  }
}

get_num_and_abundance_taxons_containing_content = function(taxon_table, presence_table, content_table){
  setkeyv(taxon_table, colnames(taxon_table)[1])
  presence_table_col_sums = colSums(presence_table[,2:ncol(presence_table),with=F])
  content_table_col_sums = colSums(content_table[,2:ncol(content_table),with=F])
  mean_num_containing_taxons_presence = mean(presence_table_col_sums, na.rm=T)
  mean_num_containing_taxons_copy_number = mean(content_table_col_sums, na.rm=T)
  sd_num_containing_taxons_presence = sd(presence_table_col_sums, na.rm=T)
  sd_num_containing_taxons_copy_number = sd(content_table_col_sums, na.rm=T)
  cov_num_containing_taxons_presence = sd_num_containing_taxons_presence/mean_num_containing_taxons_presence
  cov_num_containing_taxons_copy_number = sd_num_containing_taxons_copy_number/mean_num_containing_taxons_copy_number
  
  abundance_weighted_presence_table = presence_table[,lapply(.SD, function(col){
    if (typeof(col[1]) != "character"){
      return(col * unlist(taxon_table[presence_table[[colnames(taxon_table)[1]]],2,with=F]))
    } else {
      return(col)
    }
  }),]
  
  abundance_weighted_content_table = content_table[,lapply(.SD, function(col){
    if (typeof(col[1]) != "character"){
      return(col * unlist(taxon_table[content_table[[colnames(taxon_table)[1]]],2,with=F]))
    } else {
      return(col)
    }
  }),]
  
  abundance_weighted_presence_table_col_sums = colSums(abundance_weighted_presence_table[,2:ncol(abundance_weighted_presence_table),with=F])
  abundance_weighted_content_table_col_sums = colSums(abundance_weighted_content_table[,2:ncol(abundance_weighted_content_table),with=F])
  mean_abund_containing_taxons_presence = mean(abundance_weighted_presence_table_col_sums, na.rm=T)
  mean_abund_containing_taxons_copy_number = mean(abundance_weighted_content_table_col_sums, na.rm=T)
  sd_abund_containing_taxons_presence = sd(abundance_weighted_presence_table_col_sums, na.rm=T)
  sd_abund_containing_taxons_copy_number = sd(abundance_weighted_content_table_col_sums, na.rm=T)
  cov_abund_containing_taxons_presence = sd_abund_containing_taxons_presence/mean_abund_containing_taxons_presence
  cov_abund_containing_taxons_copy_number = sd_abund_containing_taxons_copy_number/mean_abund_containing_taxons_copy_number
  
  return(data.table(sample = colnames(taxon_table)[2], mean_num_containing_taxons = mean_num_containing_taxons_presence, mean_num_containing_taxons_copy_number_weighted = mean_num_containing_taxons_copy_number, cov_num_containing_taxons = cov_num_containing_taxons_presence, cov_num_containing_taxons_copy_number_weighted = cov_num_containing_taxons_copy_number, mean_abund_containing_taxons = mean_abund_containing_taxons_presence, mean_abund_containing_taxons_copy_number_weighted = mean_abund_containing_taxons_copy_number, cov_abund_containing_taxons = cov_abund_containing_taxons_presence, cov_abund_containing_taxons_copy_number_weighted = cov_abund_containing_taxons_copy_number))
}

get_function_occurrence_similarity = function(taxon_table, content_table){
  content_table_mat = as.matrix(content_table[,2:ncol(content_table),with=F])
  
  function_occurrence_cosine_dist_mat = as.dist(cosine(content_table_mat))
  function_occurrence_euclidean_dist_mat = vegdist(t(content_table_mat), method = "euclidean")
  function_occurrence_jaccard_dist_mat = vegdist(t(content_table_mat), method = "jaccard")
  function_occurrence_bray_dist_mat = vegdist(t(content_table_mat), method = "bray")
  function_occurrence_manhattan_dist_mat = vegdist(t(content_table_mat), method = "manhattan")
  function_occurrence_pearson_correlation_dist_mat = as.dist(cor(t(content_table_mat), method = "pearson"))
  function_occurrence_spearman_correlation_dist_mat = as.dist(cor(t(content_table_mat), method = "spearman"))
  
  mean_function_occurrence_cosine_dist = mean(function_occurrence_cosine_dist_mat, na.rm=T)
  cov_function_occurrence_cosine_dist = sd(function_occurrence_cosine_dist_mat, na.rm=T)/mean_function_occurrence_cosine_dist
  mean_function_occurrence_euclidean_dist = mean(function_occurrence_euclidean_dist_mat, na.rm=T)
  cov_function_occurrence_euclidean_dist = sd(function_occurrence_euclidean_dist_mat, na.rm=T)/mean_function_occurrence_euclidean_dist
  mean_function_occurrence_jaccard_dist = mean(function_occurrence_jaccard_dist_mat, na.rm=T)
  cov_function_occurrence_jaccard_dist = sd(function_occurrence_jaccard_dist_mat, na.rm=T)/mean_function_occurrence_jaccard_dist
  mean_function_occurrence_bray_dist = mean(function_occurrence_bray_dist_mat, na.rm=T)
  cov_function_occurrence_bray_dist = sd(function_occurrence_bray_dist_mat, na.rm=T)/mean_function_occurrence_bray_dist
  mean_function_occurrence_manhattan_dist = mean(function_occurrence_manhattan_dist_mat, na.rm=T)
  cov_function_occurrence_manhattan_dist = sd(function_occurrence_manhattan_dist_mat, na.rm=T)/mean_function_occurrence_manhattan_dist
  mean_function_occurrence_pearson_correlation_dist = mean(function_occurrence_pearson_correlation_dist_mat, na.rm=T)
  cov_function_occurrence_pearson_correlation_dist = sd(function_occurrence_pearson_correlation_dist_mat, na.rm=T)/mean_function_occurrence_pearson_correlation_dist
  mean_function_occurrence_spearman_correlation_dist = mean(function_occurrence_spearman_correlation_dist_mat, na.rm=T)
  cov_function_occurrence_spearman_correlation_dist = sd(function_occurrence_spearman_correlation_dist_mat, na.rm=T)/mean_function_occurrence_spearman_correlation_dist
  
  
  return(data.table(sample = colnames(taxon_table)[2], mean_function_occurrence_cosine_similarity = mean_function_occurrence_cosine_dist, cov_function_occurrence_cosine_similarity = cov_function_occurrence_cosine_dist, mean_function_occurrence_euclidean_dist = mean_function_occurrence_euclidean_dist, cov_function_occurrence_euclidean_dist = cov_function_occurrence_euclidean_dist, mean_function_occurrence_jaccard_dist = mean_function_occurrence_jaccard_dist, cov_function_occurrence_jaccard_dist = cov_function_occurrence_jaccard_dist, mean_function_occurrence_bray_dist = mean_function_occurrence_bray_dist, cov_function_occurrence_bray_dist = cov_function_occurrence_bray_dist, mean_function_occurrence_manhattan_dist = mean_function_occurrence_manhattan_dist, cov_function_occurrence_manhattan_dist = cov_function_occurrence_manhattan_dist, mean_function_occurrence_pearson_correlation_dist = mean_function_occurrence_pearson_correlation_dist, cov_function_occurrence_pearson_correlation_dist = cov_function_occurrence_pearson_correlation_dist, mean_function_occurrence_spearman_correlation_dist = mean_function_occurrence_spearman_correlation_dist, cov_function_occurrence_spearman_correlation_dist = cov_function_occurrence_spearman_correlation_dist))
}