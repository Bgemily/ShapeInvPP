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


# User input setup --------------------------------------------------------

N_trial_total = 20
split = 2

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------

### ICL ###########
### Parameters' possible values:
N_spks_total_list = list(50, 100, 150, 200, 250)
N_replicate_list = list(1,2,3,4,5)
N_node_list = list(100, 200, 300, 400, 500)
clus_sep_list = list(1.7, 1.8, 1.9, 2.0)

top_level_folder = "../Results/Rdata"
setup = 'ICL_Nclus4_v1.2'
default_setting = 'N_spks_total=50,N_node=100,clus_sep=1.7'
### Save estimated densities
for (. in 1:1) {
  method = 'timeshifts_est_v1.2'
  for (freq_trun in c(10)){
    ### interaction(clus_sep, N_spks_total)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_spks_total in 1:length(N_spks_total_list)) {
        N_spks_total = N_spks_total_list[[id_N_spks_total]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = 100,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = N_spks_total,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=TRUE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_spks_total"
        param_value = N_spks_total
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
    ### interaction(clus_sep, N_replicate)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_replicate in 1:length(N_replicate_list)) {
        N_replicate = N_replicate_list[[id_N_replicate]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = 100,
                               N_replicate = N_replicate,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = 50,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=TRUE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_replicate"
        param_value = N_replicate
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
    ### interaction(clus_sep, N_replicate)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_node in 1:length(N_node_list)) {
        N_node = N_node_list[[id_N_node]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = N_node,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = 50,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=TRUE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_node"
        param_value = N_node
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
  }
}

### NOT save estimated densities
for (. in 1:split) {
  method = 'timeshifts_est_v1.2'
  for (freq_trun in c(10)){
    ### interaction(clus_sep, N_spks_total)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_spks_total in 1:length(N_spks_total_list)) {
        N_spks_total = N_spks_total_list[[id_N_spks_total]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = 100,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = N_spks_total,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=FALSE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_spks_total"
        param_value = N_spks_total
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
    ### interaction(clus_sep, N_replicate)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_replicate in 1:length(N_replicate_list)) {
        N_replicate = N_replicate_list[[id_N_replicate]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = 100,
                               N_replicate = N_replicate,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = 50,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=FALSE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_replicate"
        param_value = N_replicate
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
    ### interaction(clus_sep, N_replicate)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_node in 1:length(N_node_list)) {
        N_node = N_node_list[[id_N_node]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = N_node,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = 50,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=FALSE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_node"
        param_value = N_node
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
  }
}



