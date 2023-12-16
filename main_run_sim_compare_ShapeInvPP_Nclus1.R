#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)

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
setup = 'Compare_methods_Nclus1_v3.2.7'
method = 'shape_inv_pp'

### Parameters' possible values:
timeshift_trial_max_list = list(0, 0.05, 0.1, 0.2)
N_trial_list = list(1,2,3,4,5,6,7,8,9,10)
timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2 )
N_subj_list = list(25, 50, 75, 100, 125, 150, 175, 200, 225, 250)
key_times_vec_list = list(c(-1,0-0.2,1.5), c(-1,0.04-0.2,1.5), c(-1,0.08-0.2,1.5),
                          c(-1,0.12-0.2,1.5), c(-1,0.16-0.2,1.5), c(-1,0.2-0.2,1.5))
if (TRUE) {
  for (timeshift_trial_max in timeshift_trial_max_list){
    default_setting = paste0("timeshift_trial_max=",timeshift_trial_max,",", 
                             'N_spks_total=150,N_subj=25,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2')
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
                                   N_trial = N_trial,
                                   N_subj = N_subj_list[[1]],
                                   N_clus = 1, 
                                   N_component_true = 2,
                                   N_spks_total = 150,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   timeshift_trial_max = timeshift_trial_max,
                                   t_vec = seq(-1,1.5,0.01),
                                   ### params when N_clus==4:
                                   clus_sep = 1.4,
                                   ### Parameters for algorithms
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
  
  for (use_true_timeshift in c(FALSE, TRUE)){
    default_setting = paste0("use_true_timeshift_", use_true_timeshift,",", 
                             'N_spks_total=150,N_subj=25,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2')
    for (id_N_split in 1:N_split) {
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_N_subj in 1:length(N_subj_list)) {
        N_subj = N_subj_list[[id_N_subj]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_shapeinvpp(SEED = SEED, 
                                   N_trial = 1,
                                   N_subj = N_subj,
                                   N_clus = 1, 
                                   N_component_true = 2,
                                   N_spks_total = 150,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   timeshift_trial_max = 0,
                                   t_vec = seq(-1,1.5,0.01),
                                   ### params when N_clus==4:
                                   clus_sep = 1.4,
                                   ### Parameters for algorithms
                                   freq_trun = 10,
                                   gamma = 1,
                                   N_component = 2,
                                   key_times_vec = key_times_vec_list[[1]],
                                   fix_timeshift = use_true_timeshift,
                                   use_true_timeshift = use_true_timeshift,
                                   fix_membership = FALSE,
                                   save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_subj"
        param_value = N_subj
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
    default_setting = paste0('N_trial=2,N_spks_total=150,N_subj=25,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2')
    for (id_N_split in 1:N_split) {
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_N_subj in 1:length(N_subj_list)) {
        N_subj = N_subj_list[[id_N_subj]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_shapeinvpp(SEED = SEED, 
                                   N_trial = 2,
                                   N_subj = N_subj,
                                   N_clus = 1, 
                                   N_component_true = 2,
                                   N_spks_total = 150,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   timeshift_trial_max = 0.05,
                                   t_vec = seq(-1,1.5,0.01),
                                   ### params when N_clus==4:
                                   clus_sep = 1.4,
                                   ### Parameters for algorithms
                                   freq_trun = 10,
                                   gamma = 1,
                                   N_component = 2,
                                   key_times_vec = key_times_vec_list[[1]],
                                   fix_timeshift = FALSE,
                                   fix_membership = FALSE,
                                   save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_subj"
        param_value = N_subj
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
  
}

if(TRUE){
  timeshift_trial_max_list = list(0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2)
  timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2, c(1/32/4, 1/32)*3,
                                     c(1/32/4, 1/32)*4, c(1/32/4, 1/32)*5,
                                     c(1/32/4, 1/32)*6, c(1/32/4, 1/32)*7 )
  
  default_setting = 'N_trial=2,N_spks_total=150,N_subj=25,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2'
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
                                 N_trial = 2,
                                 N_subj = 25,
                                 N_clus = 1, 
                                 N_component_true = 2,
                                 N_spks_total = 150,
                                 timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                 timeshift_trial_max = timeshift_trial_max,
                                 t_vec = seq(-1,1.5,0.01),
                                 clus_sep = 1.4,
                                 ### Parameters for algorithms
                                 freq_trun = 10,
                                 gamma = 1,
                                 N_component = 2,
                                 key_times_vec = c(-1,0-0.2,1.5),
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


