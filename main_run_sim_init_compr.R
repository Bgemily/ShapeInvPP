#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
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
save_res_details = FALSE

top_level_folder = "../Results/Rdata"
setup = 'Init_compr_v4.1'
default_setting = 'timeshift_trial_max=0.05,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.7,key_time_comp2=-0.2'

### Parameters' possible values:
N_restart_algo_list = list(1, 3, 5)
N_trial_list = list(1,2,3,4,5,6,7,8,9,10)
timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2 )
N_subj_list = list(40 )
key_times_vec_list = list(c(-1,0-0.2,1.5) )


# Test proposed init scheme -----
for (id_N_restart in 1:length(N_restart_algo_list)){
  N_restart = N_restart_algo_list[[id_N_restart]]
  method = paste0('shape_inv_pp_', 'our_init_Nrestart_algo', as.character(N_restart) )
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_trial in 1:length(N_trial_list)) {
      N_trial = N_trial_list[[id_N_trial]]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = 100*id_N_split + j
        tryCatch(main_shapeinvpp(SEED = SEED, 
                                 N_trial = N_trial,
                                 N_subj = N_subj_list[[1]],
                                 N_clus = 4, 
                                 N_component_true = 2,
                                 N_spks_total = 150,
                                 timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                 timeshift_trial_max = 0.05,
                                 t_vec = seq(-1,1.5,0.01),
                                 ### params when N_clus==4:
                                 clus_sep = 0.7,
                                 ### Parameters for algorithms
                                 rand_init = FALSE,
                                 N_restart = N_restart,
                                 N_start_kmean = 5,
                                 freq_trun = 10,
                                 gamma = 1,
                                 N_component = 2,
                                 key_times_vec = key_times_vec_list[[1]],
                                 fix_timeshift = FALSE,
                                 fix_membership = FALSE,
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


# Test random init scheme -----
for (id_N_restart in 1:length(N_restart_algo_list)){
  N_restart = N_restart_algo_list[[id_N_restart]]
  method = paste0('shape_inv_pp_', 'rand_init_Nrestart_algo', as.character(N_restart) )
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_trial in 1:length(N_trial_list)) {
      N_trial = N_trial_list[[id_N_trial]]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = 100*id_N_split + j
        tryCatch(main_shapeinvpp(SEED = SEED, 
                                 N_trial = N_trial,
                                 N_subj = N_subj_list[[1]],
                                 N_clus = 4, 
                                 N_component_true = 2,
                                 N_spks_total = 150,
                                 timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                 timeshift_trial_max = 0.05,
                                 t_vec = seq(-1,1.5,0.01),
                                 ### params when N_clus==4:
                                 clus_sep = 0.7,
                                 ### Parameters for algorithms
                                 rand_init = TRUE,
                                 N_restart = N_restart,
                                 N_start_kmean = 5,
                                 freq_trun = 10,
                                 gamma = 1,
                                 N_component = 2,
                                 key_times_vec = key_times_vec_list[[1]],
                                 fix_timeshift = FALSE,
                                 fix_membership = FALSE,
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


