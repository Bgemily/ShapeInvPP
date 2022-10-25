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

N_simtrial_total = 20
split = 2

N_simtrial = N_simtrial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
test_var_timeshift = FALSE
test_N_spks = TRUE
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Verify_lower_bound'

### Parameters' possible values:
timeshift_max_vec_list = list(c(1/4, 1/16), c(1/4, 1/16)*1.5, c(1/4, 1/16)*2,
                              c(1/4, 1/16)*2.25, c(1/4, 1/16)*2.5, c(1/4, 1/16)*2.75, c(1/4, 1/16)*3,
                              c(1/4, 1/16)*3.25, c(1/4, 1/16)*3.5, c(1/4, 1/16)*3.75,
                              c(1/4, 1/16)*0.5, c(1/4, 1/16)*0.75, 
                              c(1/4, 1/16)*0.25, c(1/4, 1/16)*0.125,
                              c(1/4, 1/16)*1.25, c(1/4, 1/16)*1.75)
N_spks_list = list(50, 75, 100, 125, 150, 175, 200)

if (test_var_timeshift) {
  for (method in c('shape_inv_pp', 'v_eq_0')) {
    fix_timeshift = ifelse(method == 'v_eq_0', yes = TRUE, no = FALSE)
    default_setting = 'N_spks_total=30,N_subj=100,N_clus=1,N_comp=1'
    for (id_split in 1:split) {
      if (save_res_details & (id_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_timeshift_max_vec in 1:length(timeshift_max_vec_list)) {
        timeshift_max_vec = timeshift_max_vec_list[[id_timeshift_max_vec]]
        results <- foreach(j = 1:N_simtrial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED, 
                               N_subj = 100,
                               N_clus = 1, 
                               N_component_true = 1,
                               N_spks_total = 30,
                               timeshift_max_vec = timeshift_max_vec,
                               t_vec = seq(-1,1,0.01),
                               ### Parameters for algorithms
                               freq_trun = 10,
                               N_component = 1,
                               key_times_vec = c(-1,1),
                               fix_timeshift = fix_timeshift,
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
        
        now_simtrial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_simtrial', N_simtrial, '_', now_simtrial, '.Rdata'))
        rm(results)
      }
    }
    
  }
  
}


if (test_N_spks) {
  method = 'shape_inv_pp'
  default_setting = 'N_spks_total_varies,N_subj=100,N_clus=1,N_comp=1'
  for (id_split in 1:split) {
    if (save_res_details & (id_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    for (id_N_spks in 1:length(N_spks_list)) {
      N_spks_total = N_spks_list[[id_N_spks]]
      results <- foreach(j = 1:N_simtrial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED, 
                             N_subj = 100,
                             N_clus = 1, 
                             N_component_true = 1,
                             N_spks_total = N_spks_total,
                             timeshift_max_vec = c(1/4, 1/16),
                             t_vec = seq(-1,1,0.01),
                             ### Parameters for algorithms
                             freq_trun = 10,
                             N_component = 1,
                             key_times_vec = c(-1,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = save_center_pdf_array),
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
      
      now_simtrial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_simtrial', N_simtrial, '_', now_simtrial, '.Rdata'))
      rm(results)
    }
  }
  
}

