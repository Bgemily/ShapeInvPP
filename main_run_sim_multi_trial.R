#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)
library(fdapace)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)
library(parallel)

# User input setup --------------------------------------------------------

N_replicate_total = 20
N_split = 2

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Multi_trial_v2'

### Parameters' possible values:
timeshift_subj_max_vec_list = list( c(1/4, 1/16)*0.1, c(1/4, 1/16)*0.3,
                                    c(1/4, 1/16)*0.5,
                                    c(1/4, 1/16)*0.75, c(1/4, 1/16)*1 )
timeshift_trial_max_list = list(0, 1/16, 1/8, 3/16, 1/4, 5/16, 3/8, 1/2)

for (id_timeshift_subj_max_vec in 1:length(timeshift_subj_max_vec_list)) {
  timeshift_subj_max_vec = timeshift_subj_max_vec_list[[id_timeshift_subj_max_vec]]
  method = paste0('shape_inv_pp','_timeshift_subj_max_vec_',
                  paste0(timeshift_subj_max_vec, collapse = '_'))
  default_setting = 'N_spks_total=100,N_subj=100,N_trial=10,N_clus=4,clus_sep=1.3,N_comp=2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_timeshift_trial_max in 1:length(timeshift_trial_max_list)) {
      timeshift_trial_max = timeshift_trial_max_list[[id_timeshift_trial_max]]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_shapeinvpp(SEED = SEED, 
                             N_trial = 10,
                             N_subj = 100,
                             N_clus = 4, 
                             N_component_true = 2,
                             N_spks_total = 100,
                             timeshift_subj_max_vec = timeshift_subj_max_vec,
                             timeshift_trial_max = timeshift_trial_max,
                             t_vec = seq(-1,1,0.01),
                             clus_sep = 1.3,
                             ### Parameters for algorithms
                             freq_trun = 10,
                             N_component = 2,
                             key_times_vec = c(-1,0,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "timeshift_trial_max"
      param_value = timeshift_trial_max
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method, 
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
      rm(results)
    }
  }
}




