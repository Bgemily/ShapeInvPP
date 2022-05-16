### Generate data, run our algorithm, and output measurements of errors.
 
main_v5_pdf = function(### Parameters for generative model
                        SEED, 
                        N_node = 100,
                        N_spks = 40,
                        N_clus=2, 
                        conn_patt_sep = 4,
                        time_shift_rad = 0.1,
                        v0 = 0.5, v1 = 0.5,
                        t_vec=seq(-v1, v0, length.out=400),
                        ### Parameters for algorithms
                        N_component=2,
                        gamma=0,
                        MaxIter=10,
                        N_clus_min=N_clus, N_clus_max=N_clus,
                        fix_timeshift=FALSE,
                        save_center_pdf_array=FALSE,
                        ### Unused
                        jitter_time_rad = 10, max_iter=50,
                        conv_thres=1e-2, 
                        freq_trun_vec=NULL,
                        opt_radius=total_time/2,
                        ...)
{
  
  # Generate data -------------------------------------------------------
  
  network_list = generate_data(SEED=SEED,
                               N_node=N_node, 
                               N_spks = N_spks,
                               conn_patt_sep = conn_patt_sep,
                               v0 = v0, v1 = v1,
                               t_vec = t_vec,
                               time_shift_rad = time_shift_rad )
  
  spks_time_mlist = network_list$spks_time_mlist
  stim_onset_vec = network_list$stim_onset_vec
  
  center_density_array_true = network_list$center_density_array_true
  mem_true_vec = network_list$mem_true_vec
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$v_vec_list
  

  # Fit model for various cluster number ------------------------------------
  
  res_list = list()
  for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
    res_list[[ind_N_clus]] = list()
    N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
    
    ### Get initialization -----------
    res = get_init(spks_time_mlist = spks_time_mlist, 
                   stim_onset_vec = stim_onset_vec,
                   N_clus = N_clus_tmp,
                   N_component = N_component,
                   v0 = v0, v1 = v1,
                   t_vec = t_vec, 
                   fix_timeshift = fix_timeshift)
    clusters_list_init = res$clusters_list
    v_vec_list_init = res$v_vec_list
    
    # Apply algorithm ---------
    
    clusters_list_init -> clusters_list_est
    v_vec_list_init -> v_vec_list_est

    time_start = Sys.time()
    ### Estimation z,v,f based on pdf
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist, 
                         stim_onset_vec = stim_onset_vec, 
                         clusters_list_init = clusters_list_est, 
                         v_vec_list_init = v_vec_list_est,
                         N_component = N_component, 
                         freq_trun = Inf,
                         MaxIter = MaxIter,
                         gamma=gamma,
                         t_vec=t_vec,
                         fix_timeshift = fix_timeshift,
                         ...)
    time_end = Sys.time()
    time_estimation = time_end - time_start
    time_estimation = as.numeric(time_estimation, units='secs')
    
    res$clusters_list -> clusters_list_est
    res$v_vec_list -> v_vec_list_est
    res$loss_history -> loss_history
    res$N_iteration -> N_iteration
    res$center_intensity_array -> center_intensity_array_est
    res$center_density_array -> center_density_array_est
    res$center_Nspks_mat -> center_Nspks_mat_est
    
    
    # Save results of N_clus_tmp ----------------------------------------------
    res_list[[ind_N_clus]] = res
    
  }
  
  

  # Select best cluster number using ICL ------------------------------------
  res_select_model = select_model(spks_time_mlist = spks_time_mlist, 
                                  stim_onset_vec = stim_onset_vec, 
                                  result_list = res_list)
  cand_N_clus_vec = N_clus_min:N_clus_max
  N_clus_est = cand_N_clus_vec[res_select_model$id_best_res]
  ICL_vec = res_select_model$ICL_vec 
  compl_log_lik_vec = res_select_model$compl_log_lik_vec 
  log_lik_vec = res_select_model$log_lik_vec
  penalty_vec = res_select_model$penalty_vec
  
  
  # Retrieve estimation results of the best cluster number ------------------
  res = res_list[[res_select_model$id_best_res]]
  
  res$clusters_list -> clusters_list_est
  res$v_vec_list -> v_vec_list_est
  res$center_intensity_array -> center_intensity_array_est
  res$center_density_array -> center_density_array_est
  res$center_Nspks_mat -> center_Nspks_mat_est
  res$loss_history -> loss_history
  res$N_iteration -> N_iteration
  
  
  # Compute estimation error ------------------------------------------------
  if (N_clus_est == N_clus) {
    # Compute errors of conn patts, i.e. F ---------
    
    ### Match clusters 
    res = find_permn(center_density_array_from = center_density_array_est,
                     center_density_array_to = center_density_array_true)
    permn = res$permn
    center_density_array_est_permn = center_density_array_est[permn, , ,drop=FALSE]
    
    
    ### Calculate distance 
    dist_mat = matrix(nrow=N_clus, ncol=N_component)
    for (id_clus in 1:N_clus) {
      for (id_component in 1:N_component) {
        dist_mat[id_clus,id_component] = sqrt(sum( (center_density_array_est_permn[id_clus,id_component,] - 
                                                      center_density_array_true[id_clus,id_component,])^2 * 
                                                     (t_vec[2]-t_vec[1]) ))
      }
    }
    F_mean_sq_err = mean(dist_mat^2)
    
    
    # Compute errors of clusters, i.e. Z ------------------------------------
    clusters_list_est_permn = clusters_list_est[permn]
    ARI = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est), 
                      memb_true_vec = mem_true_vec)
    

    # Compute error of time shifts, i.e. v --------------------------------------------
    v_mean_sq_err = mean((unlist(v_true_list)-unlist(v_vec_list_est))^2) 
    
    
    
  }
  else{
    ARI=NA 
    F_mean_sq_err=NA 
    v_mean_sq_err=NA
    clusters_list_est_permn=NA
    center_density_array_est_permn=NA
  }
  
  # Extract network related parameters -----------------------------------------
  data_param = list(SEED=SEED,
                    N_node=N_node, 
                    conn_patt_sep = conn_patt_sep,
                    v0 = v0, v1 = v1,
                    t_vec = t_vec,
                    time_shift_rad = time_shift_rad )
  
  
  # Output ------------------------------------------------------------------
  
  return(list(data_param=data_param, 
              # model selection result
              N_clus_est=N_clus_est, 
              correct_N_clus=I(N_clus_est==N_clus)*1, 
              ICL_vec=ICL_vec, 
              compl_log_lik_vec=compl_log_lik_vec, 
              log_lik_vec=log_lik_vec,
              penalty_vec=penalty_vec,
              # parameter estimates of best cluster number
              clusters_list_est=clusters_list_est,
              v_vec_list_est=v_vec_list_est,
              clusters_list_est_permn=clusters_list_est_permn,
              center_density_array_est_permn=switch(save_center_pdf_array, 
                                                "TRUE"=center_density_array_est_permn,
                                                "FALSE"=NULL),
              # estimation error
              ARI=ARI, 
              F_mean_sq_err=F_mean_sq_err, 
              v_mean_sq_err=v_mean_sq_err,
              # other
              time_estimation=time_estimation,
              N_iteration=N_iteration,
              loss_history=loss_history
              ))
  
}


