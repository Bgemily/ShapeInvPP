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

N_cores = 5
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'ICL_Nclus4_v3'
method = 'shape_inv_pp'

### Parameters' possible values:
gamma_vec = c(0.001, 0.1, 1,  10, 100)


default_setting = 'N_spks_total=70,N_subj=100,N_clus=4,clus_sep=1.4,key_time_comp2=-0.2'
for (id_N_split in 1:N_split) {
  if (save_res_details & (id_N_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = FALSE
  }
  for (id_gamma in 1:length(gamma_vec)) {
    gamma = gamma_vec[id_gamma]
    results <- foreach(j = 1:N_replicate) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_shapeinvpp(SEED = SEED, 
                               N_subj = 100,
                               N_clus = 4, 
                               N_component_true = 2,
                               N_spks_total = 70,
                               timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
                               t_vec = seq(-1,1.5,0.01),
                               clus_sep = 1.4,
                               ### Parameters for algorithms
                               freq_trun = 10,
                               gamma = 1,
                               N_clus_min = 2, N_clus_max = 6,
                               N_component = 2,
                               key_times_vec = c(-1,0.1-0.2,1.5),
                               fix_timeshift = FALSE,
                               fix_membership = FALSE,
                               save_center_pdf_array = save_center_pdf_array),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "gamma"
    param_value = gamma
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



