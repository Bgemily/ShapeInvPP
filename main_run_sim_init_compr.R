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

N_replicate_total = 50
N_split = 5

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------
test_random_restart = TRUE
test_algorithm_restart = FALSE
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Multi_restart_v3.1.2'

if (test_random_restart) {
  method = paste0('shape_inv_pp_', 'Nrestart_algo1_kmean5' )
  default_setting = 'N_spks_total=70,N_subj=100,N_clus=4,clus_sep=1.4,N_comp=2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    results <- foreach(id_replicate = 1:N_replicate) %dopar% {
      SEED = id_N_split * 1000 + id_replicate + 10
      print(paste0("SEED: ", SEED))
      tryCatch(main_shapeinvpp(SEED = SEED,
                               N_subj = 100,
                               N_clus = 4,
                               N_component_true = 2,
                               t_vec = seq(-1,1.5,0.01),
                               timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
                               N_spks_total = 70,
                               N_trial = 1,
                               clus_sep = 1.4,
                               ### Parameters for algorithms
                               N_restart = 1, 
                               N_start_kmean = 5,
                               MaxIter = 10,
                               conv_thres = -Inf, 
                               freq_trun = 10, 
                               gamma = 1,
                               N_component = 2,
                               key_times_vec = c(-1,0-0.2,1.5),
                               fix_timeshift = FALSE,
                               fix_membership = FALSE,
                               save_center_pdf_array = save_center_pdf_array ),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "N_trial"
    param_value = 1
    folder_path = paste0(top_level_folder, '/', setup,
                         '/', method, 
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
    rm(results)
    
  }
  
  N_restart_algo_list = list(1, 3, 10)
  for (id_method in 1:length(N_restart_algo_list)) {
    N_restart = N_restart_algo_list[[id_method]]
    method = paste0('shape_inv_pp_', 'Rand_init_v2_Nrestart_algo', as.character(N_restart) )
    default_setting = 'N_spks_total=70,N_subj=100,N_clus=4,clus_sep=1.4,N_comp=2'
    for (id_N_split in 1:N_split) {
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      
      results <- foreach(id_replicate = 1:N_replicate) %dopar% {
        SEED = id_N_split * 1000 + id_replicate + 10
        print(paste0("SEED: ", SEED))
        tryCatch(main_shapeinvpp(SEED = SEED,
                                 N_subj = 100,
                                 N_clus = 4,
                                 N_component_true = 2,
                                 t_vec = seq(-1,1.5,0.01),
                                 timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
                                 ### params when N_clus==4:
                                 N_spks_total = 70,
                                 N_trial = 1,
                                 clus_sep = 1.4,
                                 ### Parameters for algorithms
                                 rand_init = TRUE,
                                 N_restart = N_restart, 
                                 N_start_kmean = 1, 
                                 MaxIter = 10,
                                 conv_thres = -Inf, 
                                 freq_trun = 10, 
                                 gamma = 1,
                                 N_component = 2,
                                 key_times_vec = c(-1,0-0.2,1.5),
                                 fix_timeshift = FALSE,
                                 fix_membership = FALSE,
                                 save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "N_trial"
      param_value = 1
      folder_path = paste0(top_level_folder, '/', setup,
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


if (test_algorithm_restart) {
  ### Parameters' possible values:
  N_trial_list = list(1,2,3,4,5)
  N_restart_algo_N_start_kmean_list = list(c(1,5), c(1,1), c(3,1), c(5,1), c(3,5), c(10,1) )
  for (id_method in 1:length(N_restart_algo_N_start_kmean_list)) {
    N_restart = N_restart_algo_N_start_kmean_list[[id_method]][1]
    N_start_kmean = N_restart_algo_N_start_kmean_list[[id_method]][2]
    method = paste0('Nrestart', '_algo', as.character(N_restart), 
                    '_kmean', as.character(N_start_kmean) )
    default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=1.3,N_comp=2'
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
          tryCatch(main_shapeinvpp(SEED = SEED,
                               N_subj = 100,
                               N_clus = 4,
                               N_component_true = 2,
                               t_vec = seq(-1, 1, by=0.01),
                               timeshift_subj_max_vec = c(1/4, 1/16),
                               ### params when N_clus==4:
                               N_spks_total = 100,
                               N_trial = N_trial,
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
        param_name = "N_trial"
        param_value = N_trial
        folder_path = paste0(top_level_folder, '/', setup,
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
}






