#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)
library(fdasrvf)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)
library(parallel)

# User input setup --------------------------------------------------------

N_replicate_total = 500
N_split = 50

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------
test_N_component_2 = TRUE
save_res_details = FALSE

top_level_folder = "../Results/Rdata"
setup = 'Compare_methods_Nclus4_v4.3.3'
method = 'kmeans_align_use_intensity_internal_smth'

### Parameters' possible values:
N_trial_list = list(2,3,4,5,6,7,8,9,10)
timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2, c(1/32/4, 1/32)*3, 
                                   c(1/32/4, 1/32)*4, c(1/32/4, 1/32)*5,
                                   c(1/32/4, 1/32)*6, c(1/32/4, 1/32)*7)
N_subj_list = list(40, 60, 80, 100, 120, 140)
key_times_vec_list = list(c(-1,0-0.2,1.5), c(-1,0.02-0.2,1.5), c(-1,0.04-0.2,1.5), 
                          c(-1,0.06-0.2,1.5), c(-1,0.08-0.2,1.5), c(-1,0.1-0.2,1.5))
clus_sep_list = list(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)

if (TRUE) {
  default_setting = 'N_trial=2,timeshift_trial_max=0.1,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.5,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_kmeans_align(SEED = SEED, 
                                   N_trial = 2, timeshift_trial_max = 0.1,
                                   N_subj = N_subj_list[[1]],
                                   N_clus = 4, 
                                   N_component_true = 2,
                                   N_spks_total = 150,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   t_vec = seq(-1,1.5,0.01),
                                   clus_sep = clus_sep,
                                   key_times_vec = key_times_vec_list[[1]],
                                   ### Parameters for algorithms
                                   use_intensity = TRUE,
                                   smooth_data = TRUE,
                                   N_component = 1,
                                   save_center_pdf_array = save_center_pdf_array),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "clus_sep"
      param_value = clus_sep
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

if (TRUE) {
  default_setting = 'N_trial=2,timeshift_trial_max=0.1,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.5,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_trial in 1:length(N_trial_list)) {
      N_trial = N_trial_list[[id_N_trial]]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_kmeans_align(SEED = SEED, 
                                   N_trial = N_trial, timeshift_trial_max = 0.1,
                                   N_subj = N_subj_list[[1]],
                                   N_clus = 4, 
                                   N_component_true = 2,
                                   N_spks_total = 150,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   t_vec = seq(-1,1.5,0.01),
                                   clus_sep = 0.5,
                                   key_times_vec = key_times_vec_list[[1]],
                                   ### Parameters for algorithms
                                   use_intensity = TRUE, 
                                   smooth_data = TRUE,
                                   N_component = 1,
                                   save_center_pdf_array = save_center_pdf_array),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "N_trial"
      param_value = N_trial
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

