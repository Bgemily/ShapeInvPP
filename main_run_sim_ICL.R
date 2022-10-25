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

N_xxxxxxxx_total = 20
split = 2

N_xxxxxxxx = N_xxxxxxxx_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 5
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'ICL_Nclus4_v2.3'
method = 'timeshifts_est_v1.2'

### Parameters' possible values:
N_spks_total_list = list(50, 100, 150, 200, 250)
N_replicate_list = list(1,2,3,4,5)
N_subj_list = list(100, 200, 300, 400, 500)
clus_sep_list = list(1.7, 1.8, 1.9, 2.0)

default_setting = 'N_spks_total=50,N_subj=100,clus_sep=1.7'
for (id_split in 1:split) {
  if (save_res_details & (id_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = FALSE
  }
  ### interaction(clus_sep, N_spks_total)
  for (id_clus_sep in 1:length(clus_sep_list)){
    clus_sep = clus_sep_list[[id_clus_sep]]
    for (id_N_spks_total in 1:length(N_spks_total_list)) {
      N_spks_total = N_spks_total_list[[id_N_spks_total]]
      results <- foreach(j = 1:N_xxxxxxxx) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_subj = 100,
                             N_clus = 4,
                             N_component_true = 2,
                             t_vec = seq(-1, 1, by=0.01),
                             timeshift_max_vec = c(1/4, 1/16),
                             ### params when N_clus==4:
                             N_spks_total = N_spks_total,
                             clus_sep = clus_sep,
                             ### Parameters for algorithms
                             freq_trun = 10,
                             N_clus_min = 1,
                             N_clus_max = 6,
                             N_component = 2,
                             key_times_vec = c(-1,0,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name_0 = "clus_sep"
      param_value_0 = clus_sep
      param_name = "N_spks_total"
      param_value = N_spks_total
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method, 
                           '/', default_setting,
                           '/', param_name_0, '/', param_value_0,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
      rm(results)
    }
  }
  ### interaction(clus_sep, N_replicate)
  for (id_clus_sep in 1:length(clus_sep_list)){
    clus_sep = clus_sep_list[[id_clus_sep]]
    for (id_N_replicate in 1:length(N_replicate_list)) {
      N_replicate = N_replicate_list[[id_N_replicate]]
      results <- foreach(j = 1:N_xxxxxxxx) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_subj = 100,
                             N_replicate = N_replicate,
                             N_clus = 4,
                             N_component_true = 2,
                             t_vec = seq(-1, 1, by=0.01),
                             timeshift_max_vec = c(1/4, 1/16),
                             ### params when N_clus==4:
                             N_spks_total = 50,
                             clus_sep = clus_sep,
                             ### Parameters for algorithms
                             freq_trun = 10,
                             N_clus_min = 1,
                             N_clus_max = 6,
                             N_component = 2,
                             key_times_vec = c(-1,0,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name_0 = "clus_sep"
      param_value_0 = clus_sep
      param_name = "N_replicate"
      param_value = N_replicate
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method, 
                           '/', default_setting,
                           '/', param_name_0, '/', param_value_0,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
      rm(results)
    }
  }
  ### interaction(clus_sep, N_subj)
  for (id_clus_sep in 1:length(clus_sep_list)){
    clus_sep = clus_sep_list[[id_clus_sep]]
    for (id_N_subj in 1:length(N_subj_list)) {
      N_subj = N_subj_list[[id_N_subj]]
      results <- foreach(j = 1:N_xxxxxxxx) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_subj = N_subj,
                             N_clus = 4,
                             N_component_true = 2,
                             t_vec = seq(-1, 1, by=0.01),
                             timeshift_max_vec = c(1/4, 1/16),
                             ### params when N_clus==4:
                             N_spks_total = 50,
                             clus_sep = clus_sep,
                             ### Parameters for algorithms
                             freq_trun = 10,
                             N_clus_min = 1,
                             N_clus_max = 6,
                             N_component = 2,
                             key_times_vec = c(-1,0,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name_0 = "clus_sep"
      param_value_0 = clus_sep
      param_name = "N_subj"
      param_value = N_subj
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method, 
                           '/', default_setting,
                           '/', param_name_0, '/', param_value_0,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      
      now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
      rm(results)
    }
  }
}



