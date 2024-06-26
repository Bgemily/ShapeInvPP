#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)
library(parallel)

# User input setup --------------------------------------------------------

N_replicate_total = 200
N_split = 20

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)

# Run simulations ---------------------------------------------------------
test_N_component_2 = TRUE
save_res_details = FALSE

top_level_folder = "../Results/Rdata"
setup = 'LargeNtrial_Nclus4_v1'
method = 'shape_inv_pp_update_obj_func'

### Parameters' possible values:
timeshift_trial_max_list = list(0.1, 0.2, 0.3)
N_trial_list = list(10,20,30,40,50,60,70,80,90,100)
timeshift_subj_max_vec_list = list(c(1/32/4, 1/32)*2 )
N_subj_list = list(40, 60, 80, 100, 120, 140)
key_times_vec_list = list(c(-1,0-0.2,1.5), c(-1,0.04-0.2,1.5), c(-1,0.08-0.2,1.5),
                          c(-1,0.12-0.2,1.5), c(-1,0.16-0.2,1.5), c(-1,0.2-0.2,1.5))
clus_sep_list = list(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)

if (TRUE) {
  # ARI vs tau vs R
  for (timeshift_trial_max in timeshift_trial_max_list) {
    default_setting = paste0("timeshift_trial_max=",timeshift_trial_max,",", 
                             'N_spks_total=15,N_subj=40,N_clus=4,clus_sep=0.5,key_time_comp2=-0.2')
    for (id_N_split in 1:N_split){
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_N_trial in 1:length(N_trial_list)){
        N_trial = N_trial_list[[id_N_trial]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_shapeinvpp(SEED = SEED, 
                                   N_trial = N_trial,
                                   N_subj = N_subj_list[[1]],
                                   N_clus = 4, 
                                   N_component_true = 2,
                                   N_spks_total = 15,
                                   timeshift_subj_max_vec = timeshift_subj_max_vec_list[[1]],
                                   timeshift_trial_max = timeshift_trial_max,
                                   t_vec = seq(-1,1.5,0.01)+1,
                                   ### params when N_clus==4:
                                   clus_sep = 0.5,
                                   ### Parameters for algorithms
                                   freq_trun = 10,
                                   gamma = 1/100,
                                   N_component = 2,
                                   key_times_vec = key_times_vec_list[[1]]+1,
                                   fix_timeshift = FALSE,
                                   fix_membership = FALSE,
                                   save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_trial"
        param_value = N_trial
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
}
