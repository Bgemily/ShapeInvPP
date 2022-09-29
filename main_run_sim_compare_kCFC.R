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


# User input setup --------------------------------------------------------

N_trial_total = 20
split = 2

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------

### Compare with KCFC ###########
### Parameters' possible values:
clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5)

top_level_folder = "../Results/Rdata"
setup = 'Compare_methods_v1.4'
default_setting = 'N_spks_total=100,N_node=100,N_clus=4,N_comp=1'

### Save estimated densities
for (. in 1:1) {
  method = 'shape_inv_pp'
  for (id_clus_sep in 1:length(clus_sep_list)) {
    clus_sep = clus_sep_list[[id_clus_sep]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED, 
                           N_node = 100,
                           N_clus = 4, 
                           N_component_true = 1,
                           N_spks_total = 100,
                           timeshift_max_vec = c(1/8)*2,
                           ### params when N_clus==4:
                           clus_sep = clus_sep,
                           ### Parameters for algorithms
                           freq_trun = 10,
                           N_component = 1,
                           key_times_vec = c(-1,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = TRUE),
               error = function(x) print(SEED))
    }
    param_name = "clus_sep"
    param_value = clus_sep
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method, 
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  
  method = 'kcfc'
  for (id_clus_sep in 1:length(clus_sep_list)) {
    clus_sep = clus_sep_list[[id_clus_sep]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_kcfc(SEED = SEED,
                         N_node = 100,
                         N_clus = 4,
                         N_component_true = 1,
                         N_spks_total = 100,
                         timeshift_max_vec = c(1/8)*2,
                         ### params when N_clus==4:
                         clus_sep = clus_sep,
                         ### Parameters for algorithms
                         bw = 'SJ',
                         N_component = 1,
                         save_center_pdf_array = TRUE),
               error = function(x) print(SEED))
    }
    param_name = "clus_sep"
    param_value = clus_sep
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method,
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
}


### NOT save estimated densities
for (. in 1:split) {
  method = 'shape_inv_pp'
  for (id_clus_sep in 1:length(clus_sep_list)) {
    clus_sep = clus_sep_list[[id_clus_sep]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED, 
                           N_node = 100,
                           N_clus = 4, 
                           N_component_true = 1,
                           N_spks_total = 100,
                           timeshift_max_vec = c(1/8)*2,
                           ### params when N_clus==4:
                           clus_sep = clus_sep,
                           ### Parameters for algorithms
                           freq_trun = 10,
                           N_component = 1,
                           key_times_vec = c(-1,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = FALSE),
               error = function(x) print(SEED))
    }
    param_name = "clus_sep"
    param_value = clus_sep
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method, 
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  
  
  method = 'kcfc'
  for (id_clus_sep in 1:length(clus_sep_list)) {
    clus_sep = clus_sep_list[[id_clus_sep]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_kcfc(SEED = SEED,
                         N_node = 100,
                         N_clus = 4,
                         N_component_true = 1,
                         N_spks_total = 100,
                         timeshift_max_vec = c(1/8)*2,
                         ### params when N_clus==4:
                         clus_sep = clus_sep,
                         ### Parameters for algorithms
                         bw = 'SJ',
                         N_component = 1,
                         save_center_pdf_array = FALSE),
               error = function(x) print(SEED))
    }
    param_name = "clus_sep"
    param_value = clus_sep
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method,
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  
}


