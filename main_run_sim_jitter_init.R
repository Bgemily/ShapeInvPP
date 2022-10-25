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

N_replicate_total = 20
N_split = 2

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Jitter_init'

### Parameters' possible values:
N_trial_list = list(1,2,3,4,5)

for (jitter_level in c(0, 0.1, 0.2, 0.3, 0.5, 0.8)) {
  method = paste0('jitter_level_', as.character(jitter_level))
  default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=1.3,N_comp=2'
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
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_subj = 100,
                             N_clus = 4,
                             N_component_true = 2,
                             t_vec = seq(-1, 1, by=0.01),
                             timeshift_max_vec = c(1/4, 1/16),
                             ### params when N_clus==4:
                             N_spks_total = 100,
                             N_trial = N_trial,
                             clus_sep = 1.3,
                             ### Parameters for algorithms
                             jitter_level = jitter_level,
                             freq_trun = 10, 
                             N_component = 2,
                             key_times_vec = c(-1,0,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array ),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "N_trial"
      param_value = N_trial
      folder_path = paste0(top_level_folder, '/', setup,
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




