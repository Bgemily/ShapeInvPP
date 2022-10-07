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

N_trial_total = 50
split = 5

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------

### Parameters' possible values:
timeshift_max_vec_list = list(c(1/4, 1/16), c(1/4, 1/16)*0.75,
                              c(1/4, 1/16)*0.5, c(1/4, 1/16)*0.4, c(1/4, 1/16)*0.3, 
                              c(1/4, 1/16)*0.25,
                              c(1/4, 1/16)*1.25, c(1/4, 1/16)*1.5,
                              c(1/4, 1/16)*1.75, c(1/4, 1/16)*2)

top_level_folder = "../Results/Rdata"
setup = 'Non_identifiability_v2'
default_setting = 'N_spks_total=1000,N_node=100,N_clus=1,N_comp=2'

### Save estimated densities
for (. in 1:1) {
  method = 'shape_inv_pp'
  for (id_timeshift_max_vec in 1:length(timeshift_max_vec_list)) {
    timeshift_max_vec = timeshift_max_vec_list[[id_timeshift_max_vec]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED, 
                           N_node = 100,
                           N_clus = 1, 
                           N_component_true = 2,
                           N_spks_total = 1000,
                           timeshift_max_vec = timeshift_max_vec,
                           t_vec = seq(-1,1,0.005),
                           ### Parameters for algorithms
                           freq_trun = 10,
                           step_size = 5e-5,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = TRUE, use_true_timeshift = TRUE,
                           fix_membership = FALSE,
                           save_center_pdf_array = TRUE),
               error = function(x) print(SEED))
    }
    param_name = "timeshift_max_vec"
    param_value = paste0(timeshift_max_vec, collapse = '_')
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
  for (id_timeshift_max_vec in 1:length(timeshift_max_vec_list)) {
    timeshift_max_vec = timeshift_max_vec_list[[id_timeshift_max_vec]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED, 
                           N_node = 100,
                           N_clus = 1, 
                           N_component_true = 2,
                           N_spks_total = 1000,
                           timeshift_max_vec = timeshift_max_vec,
                           t_vec = seq(-1,1,0.005),
                           ### Parameters for algorithms
                           freq_trun = 10,
                           step_size = 5e-5,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = TRUE, use_true_timeshift = TRUE,
                           fix_membership = FALSE,
                           save_center_pdf_array = FALSE),
               error = function(x) print(SEED))
    }
    param_name = "timeshift_max_vec"
    param_value = paste0(timeshift_max_vec, collapse = '_')
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

