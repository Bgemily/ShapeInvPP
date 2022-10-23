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

N_trial_total = 50
split = 5

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'sensitivity_anal_v2.5'
method = 'ShapeInvPP'

### Parameters' possible values:
freq_trun_vec = c(1,5,10,15,20,25,30)

default_setting = 'N_spks_total=30,N_node=100,clus_sep=1.6'
for (id_split in 1:split) {
  if (save_res_details & (id_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = FALSE
  }
  for (id_freq_trun in 1:length(freq_trun_vec)){
    freq_trun = freq_trun_vec[id_freq_trun]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED,
                           N_node = 100,
                           N_clus = 4,
                           N_component_true = 2,
                           t_vec = seq(-1, 1, by=0.01),
                           timeshift_max_vec = c(1/4, 1/16),
                           ### params when N_clus==4:
                           N_spks_total = 30,
                           clus_sep = 1.6,
                           ### Parameters for algorithms
                           freq_trun = freq_trun,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = save_center_pdf_array ),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "freq_trun"
    param_value = freq_trun
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


default_setting = 'N_spks_total=50,N_node=100,clus_sep=1.6'
for (id_split in 1:split) {
  if (save_res_details & (id_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = FALSE
  }
  for (id_freq_trun in 1:length(freq_trun_vec)){
    freq_trun = freq_trun_vec[id_freq_trun]
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
                           clus_sep = 1.6,
                           ### Parameters for algorithms
                           freq_trun = freq_trun,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = save_center_pdf_array ),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "freq_trun"
    param_value = freq_trun
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


default_setting = 'N_spks_total=100,N_node=100,clus_sep=1.6'
for (id_split in 1:split) {
  if (save_res_details & (id_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = FALSE
  }
  for (id_freq_trun in 1:length(freq_trun_vec)){
    freq_trun = freq_trun_vec[id_freq_trun]
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED,
                           N_node = 100,
                           N_clus = 4,
                           N_component_true = 2,
                           t_vec = seq(-1, 1, by=0.01),
                           timeshift_max_vec = c(1/4, 1/16),
                           ### params when N_clus==4:
                           N_spks_total = 100,
                           clus_sep = 1.6,
                           ### Parameters for algorithms
                           freq_trun = freq_trun,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = FALSE,
                           fix_membership = FALSE,
                           save_center_pdf_array = save_center_pdf_array ),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "freq_trun"
    param_value = freq_trun
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




