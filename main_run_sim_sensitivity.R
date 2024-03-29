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

N_replicate_total = 200
N_split = 20

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = FALSE

top_level_folder = "../Results/Rdata"
setup = 'sensitivity_anal_v4.0'
method = 'ShapeInvPP'

### Parameters' possible values:
freq_trun_vec = c(2,4,6,8,10,12,14,16)

# ARI/MISE vs freq_trun, clus_sep = 0.5
if (TRUE) {
  default_setting = 'N_trial=3,timeshift_trial_max=0.3,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.5,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    } 
    for (id_freq_trun in 1:length(freq_trun_vec)){
      freq_trun = freq_trun_vec[id_freq_trun]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_shapeinvpp(SEED = SEED,
                                 N_trial = 3,
                                 N_subj = 40,
                                 N_clus = 4,
                                 N_component_true = 2,
                                 t_vec = seq(-1, 1.5, by=0.01)+1,
                                 timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
                                 timeshift_trial_max = 0.3,
                                 ### params when N_clus==4:
                                 N_spks_total = 150,
                                 clus_sep = 0.5,
                                 ### Parameters for algorithms
                                 freq_trun = freq_trun,
                                 gamma = 1/100,
                                 N_component = 2,
                                 key_times_vec = c(-1,0-0.2,1.5)+1,
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
      
      now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
      rm(results)
    }
    
  }
}
# ARI/MISE vs freq_trun, clus_sep = 0.9
if (FALSE) {
  default_setting = 'N_trial=3,timeshift_trial_max=0.3,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.9,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_freq_trun in 1:length(freq_trun_vec)){
      freq_trun = freq_trun_vec[id_freq_trun]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_shapeinvpp(SEED = SEED,
                                 N_trial = 3,
                                 N_subj = 40,
                                 N_clus = 4,
                                 N_component_true = 2,
                                 t_vec = seq(-1, 1.5, by=0.01)+1,
                                 timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
                                 timeshift_trial_max = 0.3,
                                 ### params when N_clus==4:
                                 N_spks_total = 150,
                                 clus_sep = 0.9,
                                 ### Parameters for algorithms
                                 freq_trun = freq_trun,
                                 gamma = 1/100,
                                 N_component = 2,
                                 key_times_vec = c(-1,0-0.2,1.5)+1,
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
      
      now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
      rm(results)
    }
    
  }
}




