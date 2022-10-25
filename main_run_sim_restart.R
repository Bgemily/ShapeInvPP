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

N_xxxxxxxx_total = 50
split = 5

N_xxxxxxxx = N_xxxxxxxx_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------
test_random_restart = TRUE
test_algorithm_restart = FALSE
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Multi_restart_v2.1'

if (test_random_restart) {
  ### Parameters' possible values:
  N_replicate_list = list(1,2,3,4,5)
  N_restart_algo_list = list(1, 3, 5, 10)
  for (id_method in 1:length(N_restart_algo_list)) {
    N_restart = N_restart_algo_list[[id_method]]
    method = paste0('shape_inv_pp_v2_', 'Rand_init_v2_Nrestart_algo', as.character(N_restart) )
    default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=2,N_comp=2'
    for (id_split in 1:split) {
      if (save_res_details & (id_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_N_replicate in 1:length(N_replicate_list)) {
        N_replicate = N_replicate_list[[id_N_replicate]]
        results <- foreach(id_xxxxxxxx = 1:N_xxxxxxxx) %dopar% {
          SEED = id_split * 1000 + id_N_replicate * 100 + id_xxxxxxxx + 10
          print(paste0("SEED: ", SEED))
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_subj = 100,
                               N_clus = 4,
                               N_component_true = 2,
                               t_vec = seq(-1, 1, by=0.01),
                               timeshift_max_vec = c(1/4, 1/16),
                               ### params when N_clus==4:
                               N_spks_total = 100,
                               N_replicate = N_replicate,
                               clus_sep = 2,
                               ### Parameters for algorithms
                               rand_init = TRUE,
                               N_restart = N_restart, 
                               N_start_kmean = 1,
                               freq_trun = 10, 
                               N_component = 2,
                               key_times_vec = c(-1,0,1),
                               fix_timeshift = FALSE,
                               fix_membership = FALSE,
                               save_center_pdf_array = save_center_pdf_array ),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_replicate"
        param_value = N_replicate
        folder_path = paste0(top_level_folder, '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
        rm(results)
      }
    }
  }
  method = paste0('shape_inv_pp_v2_', 'Nrestart_algo1_kmean5' )
  default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=2,N_comp=2'
  for (id_split in 1:split) {
    if (save_res_details & (id_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_replicate in 1:length(N_replicate_list)) {
      N_replicate = N_replicate_list[[id_N_replicate]]
      results <- foreach(id_xxxxxxxx = 1:N_xxxxxxxx) %dopar% {
        SEED = id_split * 1000 + id_N_replicate * 100 + id_xxxxxxxx + 10
        print(paste0("SEED: ", SEED))
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_subj = 100,
                             N_clus = 4,
                             N_component_true = 2,
                             t_vec = seq(-1, 1, by=0.01),
                             timeshift_max_vec = c(1/4, 1/16),
                             ### params when N_clus==4:
                             N_spks_total = 100,
                             N_replicate = N_replicate,
                             clus_sep = 2,
                             ### Parameters for algorithms
                             N_restart = 1, 
                             N_start_kmean = 5,
                             freq_trun = 10, 
                             N_component = 2,
                             key_times_vec = c(-1,0,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "N_replicate"
      param_value = N_replicate
      folder_path = paste0(top_level_folder, '/', setup,
                           '/', method, 
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
      rm(results)
    }
  }
}


if (test_algorithm_restart) {
  ### Parameters' possible values:
  N_replicate_list = list(1,2,3,4,5)
  N_restart_algo_N_start_kmean_list = list(c(1,5), c(1,1), c(3,1), c(5,1), c(3,5), c(10,1) )
  for (id_method in 1:length(N_restart_algo_N_start_kmean_list)) {
    N_restart = N_restart_algo_N_start_kmean_list[[id_method]][1]
    N_start_kmean = N_restart_algo_N_start_kmean_list[[id_method]][2]
    method = paste0('Nrestart', '_algo', as.character(N_restart), 
                    '_kmean', as.character(N_start_kmean) )
    default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=1.3,N_comp=2'
    for (id_split in 1:split) {
      if (save_res_details & (id_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_N_replicate in 1:length(N_replicate_list)) {
        N_replicate = N_replicate_list[[id_N_replicate]]
        results <- foreach(j = 1:N_xxxxxxxx) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_subj = 100,
                               N_clus = 4,
                               N_component_true = 2,
                               t_vec = seq(-1, 1, by=0.01),
                               timeshift_max_vec = c(1/4, 1/16),
                               ### params when N_clus==4:
                               N_spks_total = 100,
                               N_replicate = N_replicate,
                               clus_sep = 1.3,
                               ### Parameters for algorithms
                               N_restart = N_restart, 
                               N_start_kmean = N_start_kmean,
                               freq_trun = 10, 
                               N_component = 2,
                               key_times_vec = c(-1,0,1),
                               fix_timeshift = FALSE,
                               fix_membership = FALSE,
                               save_center_pdf_array = save_center_pdf_array ),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_replicate"
        param_value = N_replicate
        folder_path = paste0(top_level_folder, '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
        rm(results)
      }
    }
  }
}






