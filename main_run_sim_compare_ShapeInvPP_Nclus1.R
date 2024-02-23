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

N_replicate_total = 500
N_split = 50

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)

# Run simulations ---------------------------------------------------------
save_res_details = FALSE

top_level_folder = "../Results/Rdata"
setup = 'Compare_methods_Nclus1_v4'
method = 'shape_inv_pp'

### Parameters' possible values:
timeshift_trial_max_list = list(0.1, 0.2, 0.3)
# N_trial_list = list(2,3,4,5,6,7,8,9,10)
N_trial_list = list(2)
timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2 )
N_subj_list = list(10,20,30,40,50,60,70,80,90,100)
key_times_vec_list = list(c(-1,0-0.2,1.5), c(-1,0.04-0.2,1.5), c(-1,0.08-0.2,1.5),
                          c(-1,0.12-0.2,1.5), c(-1,0.16-0.2,1.5), c(-1,0.2-0.2,1.5))
# MISE vs tau vs R
if (TRUE) {
  for (timeshift_trial_max in timeshift_trial_max_list){
    default_setting = paste0("timeshift_trial_max=",timeshift_trial_max,",", 
                             'N_spks_total=200,N_subj=10,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2')
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
                                   N_spks_total = 200,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   timeshift_trial_max = timeshift_trial_max,
                                   t_vec = seq(-1,1.5,0.01)+1,
                                   ### params when N_clus==4:
                                   clus_sep = 1.4,
                                   ### Parameters for algorithms
                                   N_restart = 1,
                                   freq_trun = 10,
                                   gamma = 1,
                                   N_component = 2,
                                   key_times_vec = key_times_vec_list[[1]]+1,
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
}

if (TRUE) {
  # MISE vs n vs use_true_timeshift
  if (TRUE) {
    for (use_true_timeshift in c(TRUE,FALSE)){
      default_setting = paste0("use_true_timeshift_", use_true_timeshift,",", 
                               'N_spks_total=200,N_subj=10,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2')
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
                                     N_trial = N_trial_list[[1]],
                                     N_subj = N_subj,
                                     N_clus = 1, 
                                     N_component_true = 2,
                                     N_spks_total = 200,
                                     timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                     timeshift_trial_max = 0.1,
                                     t_vec = seq(-1,1.5,0.01)+1,
                                     ### params when N_clus==4:
                                     clus_sep = 1.4,
                                     ### Parameters for algorithms
                                     N_restart = 1,
                                     freq_trun = 10,
                                     gamma = 1,
                                     N_component = 2,
                                     key_times_vec = key_times_vec_list[[1]]+1,
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
    
  }
  # MISE vs n, use_true_timeshift = FALSE, N_trial = 3
  if (FALSE) {
    default_setting = paste0('N_trial=3,timeshift_trial_max=0.1,N_spks_total=200,N_subj=10,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2')
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
                                   N_trial = 3,
                                   N_subj = N_subj,
                                   N_clus = 1, 
                                   N_component_true = 2,
                                   N_spks_total = 200,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   timeshift_trial_max = 0.1,
                                   t_vec = seq(-1,1.5,0.01),
                                   ### params when N_clus==4:
                                   clus_sep = 1.4,
                                   ### Parameters for algorithms
                                   N_restart = 1,
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


if(FALSE){
  timeshift_trial_max_list = list(0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3)
  timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2, c(1/32/4, 1/32)*3,
                                     c(1/32/4, 1/32)*4, c(1/32/4, 1/32)*5,
                                     c(1/32/4, 1/32)*6, c(1/32/4, 1/32)*7 )
  
  default_setting = 'N_trial=2,N_spks_total=200,N_subj=10,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2'
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
                                 N_subj = N_subj_list[[1]],
                                 N_clus = 1, 
                                 N_component_true = 2,
                                 N_spks_total = 200,
                                 timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                 timeshift_trial_max = timeshift_trial_max,
                                 t_vec = seq(-1,1.5,0.01),
                                 clus_sep = 1.4,
                                 ### Parameters for algorithms
                                 N_restart = 1,
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


