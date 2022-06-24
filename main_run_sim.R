#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)


# User input setup --------------------------------------------------------

N_trial_total = 2000
split = 200

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------

### Nclus==1 ##############################
### Parameters' possible values:
N_node_list = list(100, 50, 30, 20, 10)
N_spks_ratio_list = list(3/2, 4/2, 5/2, 6/2, 8/2, 10/2)
sd_shrinkage_list = list(1, 1.5, 2.0, 2.5, 3.0)


top_level_folder = "../Results/Rdata"
setup = 'Nclus1'
default_setting = 'N_node=100,N_spks_ratio=1.5,sd_shrink=1'

for (. in 1:split) {
  method = 'timeshifts_est_v4'
  for (freq_trun in c(Inf)){
    ### N_node
    for (id_N_node in 1:length(N_node_list)) {
      N_node = N_node_list[[id_N_node]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = N_node,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = 3/2,
                             sd_shrinkage = 1,
                             ### Parameters for algorithms
                             fix_timeshift=FALSE,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "N_node"
      param_value = N_node
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### N_spks_ratio
    for (id_N_spks_ratio in 1:length(N_spks_ratio_list)) {
      N_spks_ratio = N_spks_ratio_list[[id_N_spks_ratio]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = 100,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = N_spks_ratio,
                             sd_shrinkage = 1,
                             ### Parameters for algorithms
                             fix_timeshift=FALSE,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "N_spks_ratio"
      param_value = N_spks_ratio
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### sd_shrinkage
    for (id_sd_shrinkage in 1:length(sd_shrinkage_list)) {
      sd_shrinkage = sd_shrinkage_list[[id_sd_shrinkage]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = 100,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = 3/2,
                             sd_shrinkage = sd_shrinkage,
                             ### Parameters for algorithms
                             fix_timeshift=FALSE,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "sd_shrinkage"
      param_value = sd_shrinkage
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
  }

  method = 'timeshifts_true_v4'
  for (freq_trun in c(Inf)){
    ### N_node
    for (id_N_node in 1:length(N_node_list)) {
      N_node = N_node_list[[id_N_node]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = N_node,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = 3/2,
                             sd_shrinkage = 1,
                             ### Parameters for algorithms
                             fix_timeshift=TRUE,
                             use_true_timeshift = TRUE,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "N_node"
      param_value = N_node
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### N_spks_ratio
    for (id_N_spks_ratio in 1:length(N_spks_ratio_list)) {
      N_spks_ratio = N_spks_ratio_list[[id_N_spks_ratio]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = 100,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = N_spks_ratio,
                             sd_shrinkage = 1,
                             ### Parameters for algorithms
                             fix_timeshift=TRUE,
                             use_true_timeshift = TRUE,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "N_spks_ratio"
      param_value = N_spks_ratio
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### sd_shrinkage
    for (id_sd_shrinkage in 1:length(sd_shrinkage_list)) {
      sd_shrinkage = sd_shrinkage_list[[id_sd_shrinkage]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = 100,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = 3/2,
                             sd_shrinkage = sd_shrinkage,
                             ### Parameters for algorithms
                             fix_timeshift=TRUE,
                             use_true_timeshift = TRUE,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "sd_shrinkage"
      param_value = sd_shrinkage
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
  }


  method = 'timeshifts_jitter_true_v4'
  for (freq_trun in c(Inf)){
    ### N_node
    for (id_N_node in 1:length(N_node_list)) {
      N_node = N_node_list[[id_N_node]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = N_node,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = 3/2,
                             sd_shrinkage = 1,
                             ### Parameters for algorithms
                             fix_timeshift=TRUE,
                             use_true_timeshift = TRUE,
                             jitter_prop_true_timeshift = 0.1,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "N_node"
      param_value = N_node
      folder_path = paste0(top_level_folder,
                           '/', setup,
                           '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### N_spks_ratio
    for (id_N_spks_ratio in 1:length(N_spks_ratio_list)) {
      N_spks_ratio = N_spks_ratio_list[[id_N_spks_ratio]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = 100,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = N_spks_ratio,
                             sd_shrinkage = 1,
                             ### Parameters for algorithms
                             fix_timeshift=TRUE,
                             use_true_timeshift = TRUE,
                             jitter_prop_true_timeshift = 0.1,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "N_spks_ratio"
      param_value = N_spks_ratio
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
    ### sd_shrinkage
    for (id_sd_shrinkage in 1:length(sd_shrinkage_list)) {
      sd_shrinkage = sd_shrinkage_list[[id_sd_shrinkage]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED,
                             N_node = 100,
                             N_clus = 1,
                             u_1 = 1, u_0 = 1,
                             ### params when N_clus==1:
                             N_spks_total = 60*10+40*10,
                             N_spks_ratio = 3/2,
                             sd_shrinkage = sd_shrinkage,
                             ### Parameters for algorithms
                             fix_timeshift=TRUE,
                             use_true_timeshift = TRUE,
                             jitter_prop_true_timeshift = 0.1,
                             save_center_pdf_array=FALSE ),
                 error = function(x) print(SEED))
      }
      param_name = "sd_shrinkage"
      param_value = sd_shrinkage
      folder_path = paste0(top_level_folder, '/', setup, '/', method,
                           '/', default_setting,
                           '/', param_name, '/', param_value)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

      now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
      save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
      rm(results)
    }
  }


}



# ### Nclus==2 ##############################
# ### Parameters' possible values:
# clus_mixture_list = list(0, 0.1, 0.2, 0.3, 0.4, 0.5)
# 
# 
# top_level_folder = "../Results/Rdata"
# setup = 'Nclus2'
# default_setting = 'N_node=100'
# 
# for (. in 1:split) {
#   method = 'timeshifts_est_v4'
#   for (freq_trun in c(Inf)){
#     ### clus_mixture
#     for (id_clus_mixture in 1:length(clus_mixture_list)) {
#       clus_mixture = clus_mixture_list[[id_clus_mixture]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100, 
#                              N_clus = 2, 
#                              u_1 = 1, u_0 = 1,
#                              ### params when N_clus==2:
#                              clus_mixture = clus_mixture,
#                              ### Parameters for algorithms
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_mixture"
#       param_value = clus_mixture
#       folder_path = paste0(top_level_folder, 
#                            '/', setup, 
#                            '/', method,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#   }
#   
#   method = 'timeshifts_true_v4'
#   for (freq_trun in c(Inf)){
#     ### clus_mixture
#     for (id_clus_mixture in 1:length(clus_mixture_list)) {
#       clus_mixture = clus_mixture_list[[id_clus_mixture]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100, 
#                              N_clus = 2, 
#                              u_1 = 1, u_0 = 1,
#                              ### params when N_clus==2:
#                              clus_mixture = clus_mixture,
#                              ### Parameters for algorithms
#                              fix_timeshift=TRUE, 
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_mixture"
#       param_value = clus_mixture
#       folder_path = paste0(top_level_folder, 
#                            '/', setup, 
#                            '/', method,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#   }
#   
#   
#   
#   
# }
# 
# 
# 
# 
# 
