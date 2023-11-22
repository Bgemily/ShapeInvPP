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

N_replicate_total = 200
N_split = 20

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------
test_N_component_2 = TRUE
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Compare_methods_Nclus1_v2.4.4'
method = 'fpca'

### Parameters' possible values:
timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2, c(1/32/4, 1/32)*3,
                                   c(1/32/4, 1/32)*4, c(1/32/4, 1/32)*5,
                                   c(1/32/4, 1/32)*6, c(1/32/4, 1/32)*7 )
# key_times_vec_list = list(c(-1,0.15-0.2,1.5), c(-1,0.17-0.2,1.5), c(-1,0.19-0.2,1.5), 
#                           c(-1,0.21-0.2,1.5), c(-1,0.23-0.2,1.5), c(-1,0.25-0.2,1.5))

if (test_N_component_2){
  default_setting = 'N_trial=2,timeshift_trial_max=0.05,use_true_density,N_subj=25,N_clus=1,clus_sep=1.4,key_time_comp2=-0.2'
  for (id_N_split in 1:N_split) {
    if (save_res_details & (id_N_split == 1)) {
      save_center_pdf_array = TRUE
    } else {
      save_center_pdf_array = FALSE
    }
    ### timeshift_subj_max_vec
    for (id_timeshift_subj_max_vec in 1:length(timeshift_subj_max_vec_list)) {
      timeshift_subj_max_vec = timeshift_subj_max_vec_list[[id_timeshift_subj_max_vec]]
      results <- foreach(j = 1:N_replicate) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_fpca(SEED = SEED,
                           N_trial = 2,
                           N_subj = 25,
                           N_clus = 1,
                           N_component_true = 2,
                           N_spks_total = 1000,
                           timeshift_subj_max_vec = timeshift_subj_max_vec,
                           timeshift_trial_max = 0.05,
                           t_vec = seq(-1,1.5,0.01),
                           clus_sep = 1.4,
                           key_times_vec = c(-1,0-0.2,1.5),
                           ### Parameters for algorithms
                           bw = 'SJ',
                           N_component = 2,
                           save_center_pdf_array = save_center_pdf_array),
                 error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
      }
      param_name = "timeshift_subj_max_vec"
      param_value = paste0(timeshift_subj_max_vec, collapse = '_')
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



