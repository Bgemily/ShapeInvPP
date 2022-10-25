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

N_replicate_total = 500
split = 50

N_replicate = N_replicate_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores = N_cores)


# Run simulations ---------------------------------------------------------
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Non_identifiability_v2.2'
method = 'shape_inv_pp'

### Parameters' possible values:
timeshift_max_vec_list = list(c(1/4, 1/16), c(1/4, 1/16)*0.75,
                              c(1/4, 1/16)*0.5, c(1/4, 1/16)*0.45,
                              c(1/4, 1/16)*0.4, c(1/4, 1/16)*0.35, 
                              c(1/4, 1/16)*0.325, c(1/4, 1/16)*0.3, 
                              c(1/4, 1/16)*0.275, c(1/4, 1/16)*0.25, 
                              c(1/4, 1/16)*1.25, c(1/4, 1/16)*1.5,
                              c(1/4, 1/16)*1.75, c(1/4, 1/16)*2)

default_setting = 'N_spks_total=100,N_subj=100,N_clus=1,N_comp=2'
for (id_split in 1:split) {
  if (save_res_details & (id_split == 1)) {
    save_center_pdf_array = TRUE
  } else {
    save_center_pdf_array = FALSE
  }
  for (id_timeshift_max_vec in 1:length(timeshift_max_vec_list)) {
    timeshift_max_vec = timeshift_max_vec_list[[id_timeshift_max_vec]]
    results <- foreach(j = 1:N_xxxxxxxx) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(main_v5_pdf(SEED = SEED, 
                           N_subj = 100,
                           N_clus = 1, 
                           N_component_true = 2,
                           N_spks_total = 100,
                           timeshift_max_vec = timeshift_max_vec,
                           t_vec = seq(-1,1,0.01),
                           ### Parameters for algorithms
                           freq_trun = 10,
                           N_component = 2,
                           key_times_vec = c(-1,0,1),
                           fix_timeshift = TRUE, use_true_timeshift = TRUE,
                           fix_membership = FALSE,
                           save_center_pdf_array = save_center_pdf_array),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "timeshift_max_vec"
    param_value = paste0(timeshift_max_vec, collapse = '_')
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method, 
                         '/', default_setting,
                         '/', param_name, '/', param_value)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    
    now_xxxxxxxx = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_xxxxxxxx', N_xxxxxxxx, '_', now_xxxxxxxx, '.Rdata'))
    rm(results)
  }
  
  
}


