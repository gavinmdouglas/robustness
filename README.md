This file contains instructions for obtaining and running this microbial community taxa-function robustness estimation pipeline. This pipeline takes as input a table of microbial read counts or abundances and produces an estimation of the community's taxa-function robustness. Additionally, various gene distribution feature values are calculated for the original community.

### DOWNLOAD:

The source code for this microbial community taxa-function robustness estimation pipeline can be obtained by cloning this git repository:

`git clone git@github.com:borenstein-lab/robustness.git`

##### Downloading required data files

PICRUSt and Greengenes files used by default by this pipeline should be downloaded from their original sources. To run the default pipeline, users should download the PICRUSt [16S rRNA gene copy number table](http://kronos.pharmacology.dal.ca/public_files/picrust/picrust_precalculated_v1.1.3/13_5/16S_13_5_precalculated.tab.gz) and [genomic content annotation table](http://kronos.pharmacology.dal.ca/public_files/picrust/picrust_precalculated_v1.1.3/13_5/ko_13_5_precalculated.tab.gz), saving them to the `data` directory. Additionally, the [Greengenes tree](https://s3.amazonaws.com/gg_sg_web/gg_13_5_otus_99_annotated.tree.gz?AWSAccessKeyId=AKIAIKZRXPOMF7SLT42A&Signature=esMW7qFmRflNJ8hpBrF4WJ%2FLlZU%3D&Expires=1519260808) should be downloaded, unzipped, and saved to the `data` directory.

### USE:

The pipeline follows three main steps:
1. Generate perturbed taxonomic and functional compositions for each initial community composition
2. Measure taxonomic and functional differences between initial and perturbed compositions
3. Calculate robustness metrics

##### 1. Perturbing community compositions

The first step is performed using the `generate_perturbations.R` script as follows:

`src/generate_perturbations.R initial_community_taxonomic_composition perturbed_taxonomic_compositions perturbed_function_compositions perturbation_weighted_unifrac_distances`

where the `initial_community_taxonomic_composition` file should have two columns, each with column headers. The header of the second column should be the name or ID for the initial community. The first column should contain the ID for a taxon while the second column should contain the corresponding taxon's relative abundance in the initial community.

There are several parameters that can be adjusted for this step, and all can be found by running

`src/generate_perturbations.R`

but the main parameters that may be of interest to modify are the set of maximum perturbation magnitudes to use (the `--perturbation_list` option), the number of replicates for each maximum perturbation magnitude (the `--num_replicates` option), the minimum and maximum range of Weighted UniFrac distances to generate perturbations for (the `--range_min` and `--range_max` options respectively), and the number of windows and perturbations per window to use when filtering the generated perturbations (the `--num_windows` and `--num_perturbations` options respectively).

#### 2. Measuring taxonomic and functional differences

The second step is performed using the `summarize_perturbations.R` script as follows:

`src/summarize_perturbations.R perturbed_taxonomic_compositions perturbed_function_compositions perturbation_weighted_unifrac_distances function_distribution_features perturbation_distances function_specific_differences`

where `perturbed_taxonomic_compositions`, `perturbed_function_compositions`, and `perturbation_weighted_unifrac_distances` should be the output from the previous step in the pipeline.

#### 3. Calculating robustness metrics

The third and final step is performed using the `calculate_robustness_factors.R` and `calcualte_function_specific_robustness_factors.R` scripts as follows:

```
src/calculate_robustness_factors.R perturbation_distances robustness_factors
src/calculate_function_specific_robustness_factors.R function_specific_differences function_specific_robustness_factors
```

where `perturbation_distances` and `function_specific_differences` should be the output from the previous step in the pipeline.