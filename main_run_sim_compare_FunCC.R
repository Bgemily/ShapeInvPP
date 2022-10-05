#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)
library(FunCC)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)


# User input setup --------------------------------------------------------

N_trial_total = 20
split = 2

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)



# Select tuning parameter -------------------------------------------------
N_component_true = 2
clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5)
res_delta_selection_list = list()
for (clus_sep in clus_sep_list) {
  # Generate synthetic data
  data_param = list(SEED = 1,
                    N_node = 100,
                    N_replicate = 1,
                    N_clus = 4, 
                    u_1 = 1, u_0 = 1,
                    t_vec = seq(-1, 1, by = 0.01),
                    t_vec_extend = seq(-1, 1, by = 0.01),
                    N_spks_total = 100,
                    timeshift_max_vec = c(1/4, 1/16),
                    clus_sep = clus_sep)
  if (N_component_true == 1) {
    data_generated = do.call(what = generate_data_Ncomp_1, args = data_param)
  } else if (N_component_true == 2) {
    data_generated = do.call(what = generate_data, args = data_param)
  }
  
  ### Prepare data for FunCC 
  key_times_vec = c(-1, 0, 1)
  spks_time_mlist = data_generated$spks_time_mlist
  density_array = array(dim = c(data_param$N_node, N_component_true, length(data_param$t_vec)))
  for (id_node in 1:data_param$N_node){
    res_smooth = density(spks_time_mlist[[id_node]], bw = 'SJ', 
                         from = min(data_param$t_vec), to = max(data_param$t_vec),
                         n = length(data_param$t_vec))
    for (id_component in 1:N_component_true) {
      y_curr_comp = res_smooth$y * I((data_param$t_vec >= key_times_vec[id_component]) & (data_param$t_vec <= key_times_vec[id_component+1]))
      density_array[id_node, id_component, ] = y_curr_comp
    }
  }
  
  # Get results for various delta
  res = FunCC_find_best_delta(fun_mat = density_array, 
                              delta_min = 0.01, delta_max = 0.2, num_delta = 2,
                              alpha = 0, beta = 0, theta = 1.25, shift.alignement = TRUE,
                              max.iter.align = 10, number = 10)
  res_delta_selection_list[[clus_sep]] = res
  
  # Save result
  top_level_folder = "../Results/Rdata"
  setup = 'Compare_methods_v1.8'
  method = 'FunCC'
  folder_path = paste0(top_level_folder,
                       '/', setup,
                       '/', method)
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
  save(res, file = paste0(folder_path, '/', 'tuning_parameter_selection', 
                          '_', 'clus_sep', clus_sep, '_', now_trial, '.Rdata'))
  
}

top_level_folder = "../Results/Rdata"
setup = 'Compare_methods_v1.8'
method = 'FunCC'
folder_path = paste0(top_level_folder,
                     '/', setup,
                     '/', method)
dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
save(res_delta_selection_list, file = paste0(folder_path, '/', 'tuning_parameter_selection', '_', now_trial, '.Rdata'))


# Run simulations ---------------------------------------------------------

if (FALSE) {
  ### Compare with FunCC ###########
  ### Parameters' possible values:
  clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5)
  
  top_level_folder = "../Results/Rdata"
  setup = 'Compare_methods_v1.6'
  default_setting = 'N_spks_total=100,N_node=100,N_clus=4,N_comp=1'
  
  ### Save estimated densities
  for (. in 1:1) {
    method = 'shape_inv_pp'
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED, 
                             N_node = 100,
                             N_clus = 4, 
                             N_component_true = 1,
                             N_spks_total = 100,
                             timeshift_max_vec = c(1/8)*2,
                             ### params when N_clus==4:
                             clus_sep = clus_sep,
                             ### Parameters for algorithms
                             freq_trun = 10,
                             step_size = 5e-5,
                             N_component = 1,
                             key_times_vec = c(-1,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = TRUE),
                 error = function(x) print(SEED))
      }
      param_name = "clus_sep"
      param_value = clus_sep
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
    
    method = 'kcfc'
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_kcfc(SEED = SEED,
                           N_node = 100,
                           N_clus = 4,
                           N_component_true = 1,
                           N_spks_total = 100,
                           timeshift_max_vec = c(1/8)*2,
                           ### params when N_clus==4:
                           clus_sep = clus_sep,
                           ### Parameters for algorithms
                           bw = 'SJ',
                           N_component = 1,
                           save_center_pdf_array = TRUE),
                 error = function(x) print(SEED))
      }
      param_name = "clus_sep"
      param_value = clus_sep
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
  }
  
  
  ### NOT save estimated densities
  for (. in 1:split) {
    method = 'shape_inv_pp'
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_v5_pdf(SEED = SEED, 
                             N_node = 100,
                             N_clus = 4, 
                             N_component_true = 1,
                             N_spks_total = 100,
                             timeshift_max_vec = c(1/8)*2,
                             ### params when N_clus==4:
                             clus_sep = clus_sep,
                             ### Parameters for algorithms
                             freq_trun = 10,
                             step_size = 5e-5,
                             N_component = 1,
                             key_times_vec = c(-1,1),
                             fix_timeshift = FALSE,
                             fix_membership = FALSE,
                             save_center_pdf_array = FALSE),
                 error = function(x) print(SEED))
      }
      param_name = "clus_sep"
      param_value = clus_sep
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
    
    
    method = 'kcfc'
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_kcfc(SEED = SEED,
                           N_node = 100,
                           N_clus = 4,
                           N_component_true = 1,
                           N_spks_total = 100,
                           timeshift_max_vec = c(1/8)*2,
                           ### params when N_clus==4:
                           clus_sep = clus_sep,
                           ### Parameters for algorithms
                           bw = 'SJ',
                           N_component = 1,
                           save_center_pdf_array = FALSE),
                 error = function(x) print(SEED))
      }
      param_name = "clus_sep"
      param_value = clus_sep
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
    
  }
  
  
  
}
