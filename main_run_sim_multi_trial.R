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
library(parallel)

# User input setup --------------------------------------------------------

N_replicate_total = 2000
N_split = 200

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Multi_trial_v2.2'
method = 'shape_inv_pp'

### Parameters' possible values:
timeshift_trial_max_list = list(0, 0.1, 0.2, 0.3, 0.4, 0.5)

for (N_trial in c(1,3,10)){
  default_setting = paste0('N_spks_total=70,N_subj=100,N_clus=4,N_trial=', N_trial, ',clus_sep=1.4,N_comp=2')
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
                                 N_trial = N_trial,
                                 N_subj = 100,
                                 N_clus = 4, 
                                 N_component_true = 2,
                                 N_spks_total = 70,
                                 timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
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







