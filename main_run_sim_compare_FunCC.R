#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(matrixcalc)
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

N_cores = 5
registerDoParallel(cores=N_cores)


# Experiment type ---------------------------------------------------------

select_tuning_parameter = FALSE
test_algorithm_performance = TRUE

# Select tuning parameter -------------------------------------------------
if (select_tuning_parameter) {
  N_component_true = 2
  clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5)
  N_trial = 10
  
  # Define one replicate
  find_best_delta_tmp = function(SEED, clus_sep){
    # Generate synthetic data
    data_param = list(SEED = SEED,
                      N_node = 100,
                      N_replicate = 1,
                      N_clus = 4, 
                      u_1 = 1, u_0 = 1,
                      t_vec = seq(-1, 1, by = 0.01),
                      t_vec_extend = seq(-1, 1, by = 0.01),
                      N_spks_total = 100,
                      timeshift_max_vec = c(1/4, 1/16),
                      clus_sep = clus_sep)
    data_generated = do.call(what = generate_data, args = data_param)
    
    ### Prepare data for FunCC 
    N_component_true = 2
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
                                delta_min = 0.001, delta_max = 0.1, num_delta = 20,
                                alpha = 0, beta = 0, theta = 1, 
                                shift.alignement = TRUE, shift.max = 0.25,
                                max.iter.align = 10, number = 10)
    return(list(SEED = SEED, res = res))
  }
  
  # Get results of multiple replicates
  for (clus_sep in clus_sep_list) {
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(find_best_delta_tmp(SEED = SEED, clus_sep = clus_sep),
               error = function(x) print(SEED))
    }
    
    # Save results
    top_level_folder = "../Results/Rdata"
    setup = 'FunCC_tuning_param'
    method = 'FunCC_v4'
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method,
                         '/', 'clus_sep',
                         '/', as.character(clus_sep))
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'tuning_parameter_selection', 
                                '_', now_trial, '.Rdata'))
    rm(results)
  }
}


# Run simulations ---------------------------------------------------------

if (test_algorithm_performance) {
  ### Compare with FunCC ###########
  ### Parameters' possible values:
  clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5)
  
  top_level_folder = "../Results/Rdata"
  setup = 'Compare_methods_v2.5'
  default_setting = 'N_spks_total=100,N_node=100,N_clus=4,N_comp=2'
  
  ### Save estimated densities
  for (. in 1:1) {
    if (FALSE) {
      method = 'shape_inv_pp'
      for (id_clus_sep in 1:length(clus_sep_list)) {
        clus_sep = clus_sep_list[[id_clus_sep]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED, 
                               N_node = 100,
                               N_clus = 4, 
                               N_component_true = 2,
                               N_spks_total = 100,
                               timeshift_max_vec = c(1/4, 1/16),
                               t_vec = seq(-1,1,0.01),
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = 10,
                               N_component = 2,
                               key_times_vec = c(-1,0,1),
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
      
    }
    
    method = 'funcc'
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      if (clus_sep %in% c(2.0, 1.9, 1.8, 1.7)) {
        delta = 0.01
      } else if (clus_sep %in% c(1.6, 1.5)) {
        delta = 0.03
      }
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_funcc(SEED = SEED, 
                            N_node = 100,
                            N_clus = 4, 
                            N_component_true = 2,
                            t_vec = seq(-1,1,by=0.01),
                            N_spks_total = 100,
                            timeshift_max_vec = c(1/4, 1/16),
                            ### params when N_clus==4:
                            clus_sep = clus_sep,
                            ### Parameters for algorithms
                            delta = delta, 
                            theta = 1,
                            bw = 'SJ',
                            N_component = 2,
                            key_times_vec = c(-1, 0, 1),
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
    if (FALSE) {
      method = 'shape_inv_pp'
      for (id_clus_sep in 1:length(clus_sep_list)) {
        clus_sep = clus_sep_list[[id_clus_sep]]
        results <- foreach(j = 1:N_trial) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_v5_pdf(SEED = SEED, 
                               N_node = 100,
                               N_clus = 4, 
                               N_component_true = 2,
                               N_spks_total = 100,
                               timeshift_max_vec = c(1/4, 1/16),
                               t_vec = seq(-1,1,0.01),
                               clus_sep = clus_sep,
                               ### Parameters for algorithms
                               freq_trun = 10,
                               N_component = 2,
                               key_times_vec = c(-1,0,1),
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
      
    }    
    
    method = 'funcc'
    for (id_clus_sep in 1:length(clus_sep_list)) {
      clus_sep = clus_sep_list[[id_clus_sep]]
      if (clus_sep %in% c(2.0, 1.9, 1.8, 1.7)) {
        delta = 0.01
      } else if (clus_sep %in% c(1.6, 1.5)) {
        delta = 0.03
      }
      results <- foreach(j = 1:N_trial) %dopar% {
        SEED = sample(1:1e7,1)
        tryCatch(main_funcc(SEED = SEED, 
                            N_node = 100,
                            N_clus = 4, 
                            N_component_true = 2,
                            t_vec = seq(-1,1,by=0.01),
                            N_spks_total = 100,
                            timeshift_max_vec = c(1/4, 1/16),
                            ### params when N_clus==4:
                            clus_sep = clus_sep,
                            ### Parameters for algorithms
                            delta = delta, 
                            theta = 1,
                            bw = 'SJ',
                            N_component = 2,
                            key_times_vec = c(-1, 0, 1),
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
