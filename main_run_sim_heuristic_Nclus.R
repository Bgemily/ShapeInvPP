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
setup = 'heuristic_Nclus_v1.0'
method = 'ShapeInvPP'

### Parameters' possible values:
gamma_vec = c(1/100)
N_clus_est_vec = c(2,3,4,5,6)
SEED_0 = (as.numeric(format(Sys.time(), "%OS4"))*10^4)*1000

# L1/L2/ARI vs N_clus_est, N_trial = 5
if (TRUE) {
  default_setting = 'N_trial=5,timeshift_trial_max=0.3,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.5,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_clus_est in 1:length(N_clus_est_vec)){
      N_clus_est = N_clus_est_vec[id_N_clus_est]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = SEED_0 + (id_N_split-1)*N_replicate + j
        tryCatch(main_shapeinvpp(SEED = SEED,
                                 N_trial = 5,
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
                                 N_clus_est = N_clus_est,
                                 freq_trun = 10,
                                 gamma = 1/100,
                                 N_component = 2,
                                 key_times_vec = c(-1,0-0.2,1.5)+1,
                                 fix_timeshift = FALSE,
                                 fix_membership = FALSE,
                                 save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "N_clus_est"
      param_value = N_clus_est
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

# N_trial = 10
if (TRUE) {
  default_setting = 'N_trial=10,timeshift_trial_max=0.3,N_spks_total=150,N_subj=40,N_clus=4,clus_sep=0.5,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_clus_est in 1:length(N_clus_est_vec)){
      N_clus_est = N_clus_est_vec[id_N_clus_est]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = SEED_0 + (id_N_split-1)*N_replicate + j
        tryCatch(main_shapeinvpp(SEED = SEED,
                                 N_trial = 10,
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
                                 N_clus_est = N_clus_est,
                                 freq_trun = 10,
                                 gamma = 1/100,
                                 N_component = 2,
                                 key_times_vec = c(-1,0-0.2,1.5)+1,
                                 fix_timeshift = FALSE,
                                 fix_membership = FALSE,
                                 save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "N_clus_est"
      param_value = N_clus_est
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
