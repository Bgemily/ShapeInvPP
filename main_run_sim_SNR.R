#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(matrixcalc)
library(fdapace)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)


# User input setup --------------------------------------------------------

N_trial_total = 20
split = 2

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 5
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Nclus4_v2.3'
method = 'timeshifts_est_v1.4'

### Parameters' possible values:
N_spks_total_list = list(50, 100, 150, 200, 250)
N_replicate_list = list(1,2,3,4,5)
N_node_list = list(100, 200, 300, 400, 500)
clus_sep_list = list(1.5, 1.6, 1.7, 1.8, 1.9, 2.0)

default_setting = 'N_spks_total=50,N_node=100,clus_sep=1.5'
for (id_split in 1:split) {
  if (save_res_details & (id_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = save_center_pdf_array
  }
  ### N_spks_total
  for (id_N_spks_total in 1:length(N_spks_total_list)) {
    N_spks_total = N_spks_total_list[[id_N_spks_total]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED,
                           N_node = 100,
                           N_clus = 4,
                           N_component_true = 2,
                           t_vec = seq(-1, 1, by=0.01),
                           timeshift_max_vec = c(1/4, 1/16),
                           ### params when N_clus==4:
                           N_spks_total = N_spks_total,
                           clus_sep = 1.5,
                           ### Parameters for algorithms
                           freq_trun = 10,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = save_center_pdf_array ),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "N_spks_total"
    param_value = N_spks_total
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
  
  ### N_replicate
  for (id_N_replicate in 1:length(N_replicate_list)) {
    N_replicate = N_replicate_list[[id_N_replicate]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED,
                           N_node = 100,
                           N_clus = 4,
                           N_component_true = 2,
                           t_vec = seq(-1, 1, by=0.01),
                           timeshift_max_vec = c(1/4, 1/16),
                           ### params when N_clus==4:
                           N_spks_total = 50,
                           N_replicate = N_replicate,
                           clus_sep = 1.5,
                           ### Parameters for algorithms
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
    
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)
  }
  
  ### N_node
  for (id_N_node in 1:length(N_node_list)) {
    N_node = N_node_list[[id_N_node]]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED,
                           N_node = N_node,
                           N_clus = 4,
                           N_component_true = 2,
                           t_vec = seq(-1, 1, by=0.01),
                           timeshift_max_vec = c(1/4, 1/16),
                           ### params when N_clus==4:
                           N_spks_total = 50,
                           clus_sep = 1.5,
                           ### Parameters for algorithms
                           freq_trun = 10,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = save_center_pdf_array ),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "N_node"
    param_value = N_node
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


