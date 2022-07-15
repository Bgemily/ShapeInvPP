#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)


# User input setup --------------------------------------------------------

N_trial_total = 10
split = 1

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------

### ICL ###########
### Parameters' possible values:
N_node_list = list(100, 200, 300, 400, 500)
clus_sep_list = list(1.5, 1.6, 1.7, 1.8, 1.9, 2.0)

top_level_folder = "../Results/Rdata"
setup = 'ICL_Nclus4_v1.1'
default_setting = 'N_spks_total=50,N_node=100,clus_sep=1.5'
### Save estimated densities
for (. in 1:1) {
  method = 'timeshifts_est_v1.2'
  for (freq_trun in c(10)){
    ### interaction(clus_sep, N_node)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_node in 1:length(N_node_list)) {
        N_node = N_node_list[[id_N_node]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = N_node,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = 50,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=TRUE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_node"
        param_value = N_node
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
  }

}

### NOT save estimated densities
for (. in 1:split) {
  method = 'timeshifts_est_v1.2'
  for (freq_trun in c(10)){
    ### interaction(clus_sep, N_node)
    for (id_clus_sep in 1:length(clus_sep_list)){
      clus_sep = clus_sep_list[[id_clus_sep]]
      for (id_N_node in 1:length(N_node_list)) {
        N_node = N_node_list[[id_N_node]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED,
                               N_node = N_node,
                               N_clus = 4,
                               u_1 = 1, u_0 = 1,
                               t_vec = seq(-1, 1, by=0.01),
                               t_vec_extend = seq(-3/2, 1, by=0.01),
                               ### params when N_clus==4:
                               N_spks_total = 50,
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = freq_trun,
                               N_clus_min = 1,
                               N_clus_max = 6,
                               fix_timeshift=FALSE,
                               save_center_pdf_array=FALSE ),
                   error = function(x) print(SEED))
        }
        param_name_0 = "clus_sep"
        param_value_0 = clus_sep
        param_name = "N_node"
        param_value = N_node
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, '_freqtrun', freq_trun,
                             '/', default_setting,
                             '/', param_name_0, '/', param_value_0,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

        now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
        rm(results)
      }
    }
  }

}



# ### Frequency truncation ###########
# ### Parameters' possible values:
# freq_trun_vec = c(1,5,10,15,20,25,30)
# 
# top_level_folder = "../Results/Rdata"
# setup = 'sensitivity_anal'
# default_setting = 'N_spks_total=500,N_node=100,clus_sep=2'
# ### Save estimated densities
# for (. in 1:1) {
#   method = 'timeshifts_est_v1.2'
#   for (id_freq_trun in 1:length(freq_trun_vec)){
#     freq_trun = freq_trun_vec[id_freq_trun]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v5_pdf(SEED = SEED,
#                            N_node = 100,
#                            N_clus = 4,
#                            u_1 = 1, u_0 = 1,
#                            t_vec = seq(-1, 1, by=0.01),
#                            t_vec_extend = seq(-3/2, 1, by=0.01),
#                            ### params when N_clus==4:
#                            N_spks_total = 500,
#                            clus_sep = 2,
#                            ### Parameters for algorithms
#                            freq_trun = freq_trun,
#                            fix_timeshift=FALSE,
#                            save_center_pdf_array=TRUE ),
#                error = function(x) print(SEED))
#     }
#     param_name = "freq_trun"
#     param_value = freq_trun
#     folder_path = paste0(top_level_folder,
#                          '/', setup,
#                          '/', method,
#                          '/', default_setting,
#                          '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#     
#   }
# }
# ### NOT save estimated densities
# for (. in 1:split) {
#   method = 'timeshifts_est_v1.2'
#   for (id_freq_trun in 1:length(freq_trun_vec)){
#     freq_trun = freq_trun_vec[id_freq_trun]
#     results <- foreach(j = 1:N_trial) %dopar% {
#       SEED = sample(1:1e7,1)
#       tryCatch(main_v5_pdf(SEED = SEED,
#                            N_node = 100,
#                            N_clus = 4,
#                            u_1 = 1, u_0 = 1,
#                            t_vec = seq(-1, 1, by=0.01),
#                            t_vec_extend = seq(-3/2, 1, by=0.01),
#                            ### params when N_clus==4:
#                            N_spks_total = 500,
#                            clus_sep = 2,
#                            ### Parameters for algorithms
#                            freq_trun = freq_trun,
#                            fix_timeshift=FALSE,
#                            save_center_pdf_array=FALSE ),
#                error = function(x) print(SEED))
#     }
#     param_name = "freq_trun"
#     param_value = freq_trun
#     folder_path = paste0(top_level_folder,
#                          '/', setup,
#                          '/', method,
#                          '/', default_setting,
#                          '/', param_name, '/', param_value)
#     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     
#     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     rm(results)
#     
#   }
# }
# 
# 
# 
# 
# 
# ### Nclus==4. Varying signal strength. ###########
# ### Parameters' possible values:
# N_spks_total_list = list(50, 100, 150, 200, 250)
# N_replicate_list = list(1,2,3,4,5)
# N_node_list = list(100, 200, 300, 400, 500)
# clus_sep_list = list(1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
# 
# top_level_folder = "../Results/Rdata"
# setup = 'Nclus4'
# default_setting = 'N_spks_total=50,N_node=100,clus_sep=1.5'
# 
# ### Save estimated densities
# for (. in 1:1) {
#   method = 'timeshifts_est_v1.2'
#   for (freq_trun in c(10,20)){
#     ### N_spks_total
#     for (id_N_spks_total in 1:length(N_spks_total_list)) {
#       N_spks_total = N_spks_total_list[[id_N_spks_total]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = N_spks_total,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun=freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_spks_total"
#       param_value = N_spks_total
#       folder_path = paste0(top_level_folder, 
#                            '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_replicate
#     for (id_N_replicate in 1:length(N_replicate_list)) {
#       N_replicate = N_replicate_list[[id_N_replicate]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              N_replicate = N_replicate,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_replicate"
#       param_value = N_replicate
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_node
#     for (id_N_node in 1:length(N_node_list)) {
#       N_node = N_node_list[[id_N_node]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = N_node,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_node"
#       param_value = N_node
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### clus_sep
#     for (id_clus_sep in 1:length(clus_sep_list)) {
#       clus_sep = clus_sep_list[[id_clus_sep]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = clus_sep,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_sep"
#       param_value = clus_sep
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
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
#   method = 'timeshifts_true_v1.2'
#   for (freq_trun in c(10,20)){
#     ### N_spks_ratio
#     for (id_N_spks_total in 1:length(N_spks_total_list)) {
#       N_spks_total = N_spks_total_list[[id_N_spks_total]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = N_spks_total,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_spks_total"
#       param_value = N_spks_total
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_replicate
#     for (id_N_replicate in 1:length(N_replicate_list)) {
#       N_replicate = N_replicate_list[[id_N_replicate]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              N_replicate = N_replicate,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_replicate"
#       param_value = N_replicate
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     
#     ### N_node
#     for (id_N_node in 1:length(N_node_list)) {
#       N_node = N_node_list[[id_N_node]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = N_node,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_node"
#       param_value = N_node
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### clus_sep
#     for (id_clus_sep in 1:length(clus_sep_list)) {
#       clus_sep = clus_sep_list[[id_clus_sep]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = clus_sep,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_sep"
#       param_value = clus_sep
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
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
#   method = 'membership_true_v1.2'
#   for (freq_trun in c(10,20)){
#     ### N_spks_ratio
#     for (id_N_spks_total in 1:length(N_spks_total_list)) {
#       N_spks_total = N_spks_total_list[[id_N_spks_total]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = N_spks_total,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_spks_total"
#       param_value = N_spks_total
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_replicate
#     for (id_N_replicate in 1:length(N_replicate_list)) {
#       N_replicate = N_replicate_list[[id_N_replicate]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              N_replicate = N_replicate,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_replicate"
#       param_value = N_replicate
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_node
#     for (id_N_node in 1:length(N_node_list)) {
#       N_node = N_node_list[[id_N_node]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = N_node,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_node"
#       param_value = N_node
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### clus_sep
#     for (id_clus_sep in 1:length(clus_sep_list)) {
#       clus_sep = clus_sep_list[[id_clus_sep]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = clus_sep,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=TRUE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_sep"
#       param_value = clus_sep
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
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
# }
# 
# 
# ### NOT save estimated densities
# for (. in 1:split) {
#   method = 'timeshifts_est_v1.2'
#   for (freq_trun in c(10,20)){
#     ### N_spks_total
#     for (id_N_spks_total in 1:length(N_spks_total_list)) {
#       N_spks_total = N_spks_total_list[[id_N_spks_total]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = N_spks_total,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun=freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_spks_total"
#       param_value = N_spks_total
#       folder_path = paste0(top_level_folder, 
#                            '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
# 
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_replicate
#     for (id_N_replicate in 1:length(N_replicate_list)) {
#       N_replicate = N_replicate_list[[id_N_replicate]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              N_replicate = N_replicate,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_replicate"
#       param_value = N_replicate
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_node
#     for (id_N_node in 1:length(N_node_list)) {
#       N_node = N_node_list[[id_N_node]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = N_node,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_node"
#       param_value = N_node
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
# 
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### clus_sep
#     for (id_clus_sep in 1:length(clus_sep_list)) {
#       clus_sep = clus_sep_list[[id_clus_sep]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = clus_sep,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_sep"
#       param_value = clus_sep
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
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
#   method = 'timeshifts_true_v1.2'
#   for (freq_trun in c(10,20)){
#     ### N_spks_ratio
#     for (id_N_spks_total in 1:length(N_spks_total_list)) {
#       N_spks_total = N_spks_total_list[[id_N_spks_total]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = N_spks_total,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_spks_total"
#       param_value = N_spks_total
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_replicate
#     for (id_N_replicate in 1:length(N_replicate_list)) {
#       N_replicate = N_replicate_list[[id_N_replicate]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              N_replicate = N_replicate,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_replicate"
#       param_value = N_replicate
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     
#     ### N_node
#     for (id_N_node in 1:length(N_node_list)) {
#       N_node = N_node_list[[id_N_node]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = N_node,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_node"
#       param_value = N_node
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### clus_sep
#     for (id_clus_sep in 1:length(clus_sep_list)) {
#       clus_sep = clus_sep_list[[id_clus_sep]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = clus_sep,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=TRUE,
#                              use_true_timeshift = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_sep"
#       param_value = clus_sep
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
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
#   method = 'membership_true_v1.2'
#   for (freq_trun in c(10,20)){
#     ### N_spks_ratio
#     for (id_N_spks_total in 1:length(N_spks_total_list)) {
#       N_spks_total = N_spks_total_list[[id_N_spks_total]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = N_spks_total,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_spks_total"
#       param_value = N_spks_total
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_replicate
#     for (id_N_replicate in 1:length(N_replicate_list)) {
#       N_replicate = N_replicate_list[[id_N_replicate]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              N_replicate = N_replicate,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_replicate"
#       param_value = N_replicate
#       folder_path = paste0(top_level_folder, '/', setup, 
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### N_node
#     for (id_N_node in 1:length(N_node_list)) {
#       N_node = N_node_list[[id_N_node]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = N_node,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = 1.5,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "N_node"
#       param_value = N_node
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
#                            '/', default_setting,
#                            '/', param_name, '/', param_value)
#       dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#       
#       now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#       save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#       rm(results)
#     }
#     
#     ### clus_sep
#     for (id_clus_sep in 1:length(clus_sep_list)) {
#       clus_sep = clus_sep_list[[id_clus_sep]]
#       results <- foreach(j = 1:N_trial) %dopar% {
#         SEED = sample(1:1e7,1)
#         tryCatch(main_v5_pdf(SEED = SEED,
#                              N_node = 100,
#                              N_clus = 4,
#                              u_1 = 1, u_0 = 1,
#                              t_vec = seq(-1, 1, by=0.01),
#                              t_vec_extend = seq(-3/2, 1, by=0.01),
#                              ### params when N_clus==4:
#                              N_spks_total = 50,
#                              clus_sep = clus_sep,
#                              ### Parameters for algorithms
#                              freq_trun = freq_trun,
#                              fix_timeshift=FALSE,
#                              fix_membership = TRUE,
#                              save_center_pdf_array=FALSE ),
#                  error = function(x) print(SEED))
#       }
#       param_name = "clus_sep"
#       param_value = clus_sep
#       folder_path = paste0(top_level_folder,
#                            '/', setup,
#                            '/', method, '_freqtrun', freq_trun,
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
# }
# 
# 
# 
# ### Nclus==1 ##############################
# ### Parameters' possible values:
# N_spks_total_list = list(1000, 300, 100, 30, 10)
# N_node_list = list(200, 100, 50, 30, 20, 10)
# N_spks_ratio_list = list(0.2, 0.25, 0.33, 0.4, 0.5, 0.67, 1, 1.5, 2, 2.5, 3, 4, 5)
# sd_shrinkage_list = list(1, 1.5, 2.0, 2.5, 3.0)
# c_1_list = list(0.075, 0.15, 0.225, 0.3, 0.375, 0.45)
# delta_1_list = list(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
# c_2_list = list(0.075, 0.15, 0.225, 0.3, 0.375, 0.45)
# delta_2_list = list(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
# c_3_list = list(0.075, 0.15, 0.225, 0.3, 0.375, 0.45)
# delta_3_list = list(0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
# 
# 
# top_level_folder = "../Results/Rdata"
# setup = 'Nclus1'
# default_setting = 'N_node=100,N_spks_total=1000,N_spks_ratio=1.5,not_save_density,v2'
# 
# for (. in 1:split) {
#   method = 'timeshifts_est_v4.1.2'
#   for (freq_trun in c(Inf)){
#     ### interaction(c_1, delta_1)
#     for (id_c_1 in 1:length(c_1_list)) {
#       c_1 = c_1_list[[id_c_1]]
#       for (id_delta_1 in 1:length(delta_1_list)){
#         delta_1 = delta_1_list[[id_delta_1]]
#         results <- foreach(j = 1:N_trial) %dopar% {
#           SEED = sample(1:1e7,1)
#           tryCatch(main_v5_pdf(SEED = SEED,
#                                N_node = 100,
#                                N_clus = 1,
#                                u_1 = 1, u_0 = 1,
#                                t_vec = seq(-1, 1, by=0.01),
#                                t_vec_extend = seq(-1*3/2, 1, by=0.01),
#                                ### params when N_clus==1:
#                                N_spks_total = 1000,
#                                N_spks_ratio = 3/2,
#                                sd_shrinkage = 1, 
#                                c_1 = c_1, delta_1 = delta_1,
#                                ### Parameters for algorithms
#                                fix_timeshift=FALSE,
#                                save_center_pdf_array=FALSE ),
#                    error = function(x) print(SEED))
#         }
#         param_name_0 = "c_1"
#         param_value_0 = c_1
#         param_name = "delta_1"
#         param_value = delta_1
#         folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                              '/', default_setting,
#                              '/', param_name_0, '/', param_value_0,
#                              '/', param_name, '/', param_value)
#         dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#         
#         now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#         save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#         rm(results)
#       }
#       
#     }
#     
#     ### interaction(c_2, delta_2)
#     for (id_c_2 in 1:length(c_2_list)) {
#       c_2 = c_2_list[[id_c_2]]
#       for (id_delta_2 in 1:length(delta_2_list)){
#         delta_2 = delta_2_list[[id_delta_2]]
#         results <- foreach(j = 1:N_trial) %dopar% {
#           SEED = sample(1:1e7,1)
#           tryCatch(main_v5_pdf(SEED = SEED,
#                                N_node = 100,
#                                N_clus = 1,
#                                u_1 = 1, u_0 = 1,
#                                t_vec = seq(-1, 1, by=0.01),
#                                t_vec_extend = seq(-1*3/2, 1, by=0.01),
#                                ### params when N_clus==1:
#                                N_spks_total = 1000,
#                                N_spks_ratio = 3/2,
#                                sd_shrinkage = 1, 
#                                c_2 = c_2, delta_2 = delta_2,
#                                ### Parameters for algorithms
#                                fix_timeshift=FALSE,
#                                save_center_pdf_array=FALSE ),
#                    error = function(x) print(SEED))
#         }
#         param_name_0 = "c_2"
#         param_value_0 = c_2
#         param_name = "delta_2"
#         param_value = delta_2
#         folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                              '/', default_setting,
#                              '/', param_name_0, '/', param_value_0,
#                              '/', param_name, '/', param_value)
#         dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#         
#         now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#         save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#         rm(results)
#       }
#       
#     }
#     
#     ### interaction(c_3, delta_3)
#     for (id_c_3 in 1:length(c_3_list)) {
#       c_3 = c_3_list[[id_c_3]]
#       for (id_delta_3 in 1:length(delta_3_list)){
#         delta_3 = delta_3_list[[id_delta_3]]
#         results <- foreach(j = 1:N_trial) %dopar% {
#           SEED = sample(1:1e7,1)
#           tryCatch(main_v5_pdf(SEED = SEED,
#                                N_node = 100,
#                                N_clus = 1,
#                                u_1 = 1, u_0 = 1,
#                                t_vec = seq(-1, 1, by=0.01),
#                                t_vec_extend = seq(-1*3/2, 1, by=0.01),
#                                ### params when N_clus==1:
#                                N_spks_total = 1000,
#                                N_spks_ratio = 3/2,
#                                sd_shrinkage = 1, 
#                                c_3 = c_3, delta_3 = delta_3,
#                                ### Parameters for algorithms
#                                fix_timeshift=FALSE,
#                                save_center_pdf_array=FALSE ),
#                    error = function(x) print(SEED))
#         }
#         param_name_0 = "c_3"
#         param_value_0 = c_3
#         param_name = "delta_3"
#         param_value = delta_3
#         folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                              '/', default_setting,
#                              '/', param_name_0, '/', param_value_0,
#                              '/', param_name, '/', param_value)
#         dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#         
#         now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#         save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#         rm(results)
#       }
#       
#     }
#     
#     # ### N_spks_ratio
#     # for (id_N_spks_ratio in 1:length(N_spks_ratio_list)) {
#     #   N_spks_ratio = N_spks_ratio_list[[id_N_spks_ratio]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5_pdf(SEED = SEED,
#     #                          N_node = 10,
#     #                          N_clus = 1,
#     #                          u_1 = 1, u_0 = 1,
#     #                          ### params when N_clus==1:
#     #                          N_spks_total = 100,
#     #                          N_spks_ratio = N_spks_ratio,
#     #                          sd_shrinkage = 1,
#     #                          ### Parameters for algorithms
#     #                          fix_timeshift=FALSE,
#     #                          save_center_pdf_array=TRUE ),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "N_spks_ratio"
#     #   param_value = N_spks_ratio
#     #   folder_path = paste0(top_level_folder, '/', setup, '/', method,
#     #                        '/', default_setting,
#     #                        '/', param_name, '/', param_value)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     #   
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#     
#     # ### interact(N_spks_total, N_node) 
#     # for (id_N_spks_total in 1:length(N_spks_total_list)) {
#     #   N_spks_total = N_spks_total_list[[id_N_spks_total]]
#     #   for (id_N_node in 1:length(N_node_list)) {
#     #     N_node = N_node_list[[id_N_node]]
#     #     results <- foreach(j = 1:N_trial) %dopar% {
#     #       SEED = sample(1:1e7,1)
#     #       tryCatch(main_v5_pdf(SEED = SEED,
#     #                            N_node = N_node,
#     #                            N_clus = 1,
#     #                            u_1 = 1, u_0 = 1,
#     #                            ### params when N_clus==1:
#     #                            N_spks_total = N_spks_total,
#     #                            N_spks_ratio = 3/2,
#     #                            sd_shrinkage = 1,
#     #                            ### Parameters for algorithms
#     #                            fix_timeshift=FALSE,
#     #                            save_center_pdf_array=TRUE ),
#     #                error = function(x) print(SEED))
#     #     }
#     #     param_name_0 = "N_spks_total"
#     #     param_value_0 = N_spks_total
#     #     param_name = "N_node"
#     #     param_value = N_node
#     #     folder_path = paste0(top_level_folder,
#     #                          '/', setup,
#     #                          '/', method,
#     #                          # '/', default_setting,
#     #                          '/', param_name_0, '/', param_value_0,
#     #                          '/', param_name, '/', param_value)
#     #     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     #     
#     #     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #     rm(results)
#     #   }
#     #   
#     # }
#     
#     # ### N_node
#     # for (id_N_node in 1:length(N_node_list)) {
#     #   N_node = N_node_list[[id_N_node]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5_pdf(SEED = SEED,
#     #                          N_node = N_node,
#     #                          N_clus = 1,
#     #                          u_1 = 1, u_0 = 1,
#     #                          ### params when N_clus==1:
#     #                          N_spks_total = 60*10+40*10,
#     #                          N_spks_ratio = 3/2,
#     #                          sd_shrinkage = 1,
#     #                          ### Parameters for algorithms
#     #                          fix_timeshift=FALSE,
#     #                          save_center_pdf_array=FALSE ),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "N_node"
#     #   param_value = N_node
#     #   folder_path = paste0(top_level_folder,
#     #                        '/', setup,
#     #                        '/', method,
#     #                        '/', default_setting,
#     #                        '/', param_name, '/', param_value)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     # 
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#     
#     # ### sd_shrinkage
#     # for (id_sd_shrinkage in 1:length(sd_shrinkage_list)) {
#     #   sd_shrinkage = sd_shrinkage_list[[id_sd_shrinkage]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5_pdf(SEED = SEED,
#     #                          N_node = 100,
#     #                          N_clus = 1,
#     #                          u_1 = 1, u_0 = 1,
#     #                          ### params when N_clus==1:
#     #                          N_spks_total = 60*10+40*10,
#     #                          N_spks_ratio = 3/2,
#     #                          sd_shrinkage = sd_shrinkage,
#     #                          ### Parameters for algorithms
#     #                          fix_timeshift=FALSE,
#     #                          save_center_pdf_array=FALSE ),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "sd_shrinkage"
#     #   param_value = sd_shrinkage
#     #   folder_path = paste0(top_level_folder, '/', setup, '/', method,
#     #                        '/', default_setting,
#     #                        '/', param_name, '/', param_value)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     # 
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#   }
# 
#   method = 'timeshifts_true_v4.1.2'
#   for (freq_trun in c(Inf)){
#     ### interaction(c_1, delta_1)
#     for (id_c_1 in 1:length(c_1_list)) {
#       c_1 = c_1_list[[id_c_1]]
#       for (id_delta_1 in 1:length(delta_1_list)){
#         delta_1 = delta_1_list[[id_delta_1]]
#         results <- foreach(j = 1:N_trial) %dopar% {
#           SEED = sample(1:1e7,1)
#           tryCatch(main_v5_pdf(SEED = SEED,
#                                N_node = 100,
#                                N_clus = 1,
#                                u_1 = 1, u_0 = 1,
#                                t_vec = seq(-1, 1, by=0.01),
#                                t_vec_extend = seq(-1*3/2, 1, by=0.01),
#                                ### params when N_clus==1:
#                                N_spks_total = 1000,
#                                N_spks_ratio = 3/2,
#                                sd_shrinkage = 1, 
#                                c_1 = c_1, delta_1 = delta_1,
#                                ### Parameters for algorithms
#                                fix_timeshift=TRUE, 
#                                use_true_timeshift = TRUE,
#                                save_center_pdf_array=FALSE ),
#                    error = function(x) print(SEED))
#         }
#         param_name_0 = "c_1"
#         param_value_0 = c_1
#         param_name = "delta_1"
#         param_value = delta_1
#         folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                              '/', default_setting,
#                              '/', param_name_0, '/', param_value_0,
#                              '/', param_name, '/', param_value)
#         dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#         
#         now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#         save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#         rm(results)
#       }
#       
#     }
# 
#     ### interaction(c_2, delta_2)
#     for (id_c_2 in 1:length(c_2_list)) {
#       c_2 = c_2_list[[id_c_2]]
#       for (id_delta_2 in 1:length(delta_2_list)){
#         delta_2 = delta_2_list[[id_delta_2]]
#         results <- foreach(j = 1:N_trial) %dopar% {
#           SEED = sample(1:1e7,1)
#           tryCatch(main_v5_pdf(SEED = SEED,
#                                N_node = 100,
#                                N_clus = 1,
#                                u_1 = 1, u_0 = 1,
#                                t_vec = seq(-1, 1, by=0.01),
#                                t_vec_extend = seq(-1*3/2, 1, by=0.01),
#                                ### params when N_clus==1:
#                                N_spks_total = 1000,
#                                N_spks_ratio = 3/2,
#                                sd_shrinkage = 1, 
#                                c_2 = c_2, delta_2 = delta_2,
#                                ### Parameters for algorithms
#                                fix_timeshift=TRUE, 
#                                use_true_timeshift = TRUE,
#                                save_center_pdf_array=FALSE ),
#                    error = function(x) print(SEED))
#         }
#         param_name_0 = "c_2"
#         param_value_0 = c_2
#         param_name = "delta_2"
#         param_value = delta_2
#         folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                              '/', default_setting,
#                              '/', param_name_0, '/', param_value_0,
#                              '/', param_name, '/', param_value)
#         dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#         
#         now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#         save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#         rm(results)
#       }
#       
#     }
#     
#     ### interaction(c_3, delta_3)
#     for (id_c_3 in 1:length(c_3_list)) {
#       c_3 = c_3_list[[id_c_3]]
#       for (id_delta_3 in 1:length(delta_3_list)){
#         delta_3 = delta_3_list[[id_delta_3]]
#         results <- foreach(j = 1:N_trial) %dopar% {
#           SEED = sample(1:1e7,1)
#           tryCatch(main_v5_pdf(SEED = SEED,
#                                N_node = 100,
#                                N_clus = 1,
#                                u_1 = 1, u_0 = 1,
#                                t_vec = seq(-1, 1, by=0.01),
#                                t_vec_extend = seq(-1*3/2, 1, by=0.01),
#                                ### params when N_clus==1:
#                                N_spks_total = 1000,
#                                N_spks_ratio = 3/2,
#                                sd_shrinkage = 1, 
#                                c_3 = c_3, delta_3 = delta_3,
#                                ### Parameters for algorithms
#                                fix_timeshift=TRUE, 
#                                use_true_timeshift = TRUE,
#                                save_center_pdf_array=FALSE ),
#                    error = function(x) print(SEED))
#         }
#         param_name_0 = "c_3"
#         param_value_0 = c_3
#         param_name = "delta_3"
#         param_value = delta_3
#         folder_path = paste0(top_level_folder, '/', setup, '/', method,
#                              '/', default_setting,
#                              '/', param_name_0, '/', param_value_0,
#                              '/', param_name, '/', param_value)
#         dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#         
#         now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#         save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#         rm(results)
#       }
#       
#     }
#     
#     # ### N_spks_ratio
#     # for (id_N_spks_ratio in 1:length(N_spks_ratio_list)) {
#     #   N_spks_ratio = N_spks_ratio_list[[id_N_spks_ratio]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5_pdf(SEED = SEED,
#     #                          N_node = 10,
#     #                          N_clus = 1,
#     #                          u_1 = 1, u_0 = 1,
#     #                          ### params when N_clus==1:
#     #                          N_spks_total = 100,
#     #                          N_spks_ratio = N_spks_ratio,
#     #                          sd_shrinkage = 1,
#     #                          ### Parameters for algorithms
#     #                          fix_timeshift=TRUE,
#     #                          use_true_timeshift = TRUE,
#     #                          save_center_pdf_array=FALSE ),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "N_spks_ratio"
#     #   param_value = N_spks_ratio
#     #   folder_path = paste0(top_level_folder, '/', setup, '/', method,
#     #                        '/', default_setting,
#     #                        '/', param_name, '/', param_value)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     # 
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#     # ### interact(N_spks_total, N_node) 
#     # for (id_N_spks_total in 1:length(N_spks_total_list)) {
#     #   N_spks_total = N_spks_total_list[[id_N_spks_total]]
#     #   for (id_N_node in 1:length(N_node_list)) {
#     #     N_node = N_node_list[[id_N_node]]
#     #     results <- foreach(j = 1:N_trial) %dopar% {
#     #       SEED = sample(1:1e7,1)
#     #       tryCatch(main_v5_pdf(SEED = SEED,
#     #                            N_node = N_node,
#     #                            N_clus = 1,
#     #                            u_1 = 1, u_0 = 1,
#     #                            ### params when N_clus==1:
#     #                            N_spks_total = N_spks_total,
#     #                            N_spks_ratio = 3/2,
#     #                            sd_shrinkage = 1,
#     #                            ### Parameters for algorithms
#     #                            fix_timeshift=TRUE,
#     #                            use_true_timeshift = TRUE,
#     #                            save_center_pdf_array=FALSE ),
#     #                error = function(x) print(SEED))
#     #     }
#     #     param_name_0 = "N_spks_total"
#     #     param_value_0 = N_spks_total
#     #     param_name = "N_node"
#     #     param_value = N_node
#     #     folder_path = paste0(top_level_folder,
#     #                          '/', setup,
#     #                          '/', method,
#     #                          # '/', default_setting,
#     #                          '/', param_name_0, '/', param_value_0,
#     #                          '/', param_name, '/', param_value)
#     #     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     #     
#     #     now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #     save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #     rm(results)
#     #   }
#     #   
#     # }
#     
#     # ### N_node
#     # for (id_N_node in 1:length(N_node_list)) {
#     #   N_node = N_node_list[[id_N_node]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5_pdf(SEED = SEED,
#     #                          N_node = N_node,
#     #                          N_clus = 1,
#     #                          u_1 = 1, u_0 = 1,
#     #                          ### params when N_clus==1:
#     #                          N_spks_total = 60*10+40*10,
#     #                          N_spks_ratio = 3/2,
#     #                          sd_shrinkage = 1,
#     #                          ### Parameters for algorithms
#     #                          fix_timeshift=TRUE,
#     #                          use_true_timeshift = TRUE,
#     #                          save_center_pdf_array=FALSE ),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "N_node"
#     #   param_value = N_node
#     #   folder_path = paste0(top_level_folder,
#     #                        '/', setup,
#     #                        '/', method,
#     #                        '/', default_setting,
#     #                        '/', param_name, '/', param_value)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     # 
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#     # ### sd_shrinkage
#     # for (id_sd_shrinkage in 1:length(sd_shrinkage_list)) {
#     #   sd_shrinkage = sd_shrinkage_list[[id_sd_shrinkage]]
#     #   results <- foreach(j = 1:N_trial) %dopar% {
#     #     SEED = sample(1:1e7,1)
#     #     tryCatch(main_v5_pdf(SEED = SEED,
#     #                          N_node = 100,
#     #                          N_clus = 1,
#     #                          u_1 = 1, u_0 = 1,
#     #                          ### params when N_clus==1:
#     #                          N_spks_total = 60*10+40*10,
#     #                          N_spks_ratio = 3/2,
#     #                          sd_shrinkage = sd_shrinkage,
#     #                          ### Parameters for algorithms
#     #                          fix_timeshift=TRUE,
#     #                          use_true_timeshift = TRUE,
#     #                          save_center_pdf_array=FALSE ),
#     #              error = function(x) print(SEED))
#     #   }
#     #   param_name = "sd_shrinkage"
#     #   param_value = sd_shrinkage
#     #   folder_path = paste0(top_level_folder, '/', setup, '/', method,
#     #                        '/', default_setting,
#     #                        '/', param_name, '/', param_value)
#     #   dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
#     # 
#     #   now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
#     #   save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
#     #   rm(results)
#     # }
#   }
# 
# 
# }



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
