### Generate data, run our algorithm, and output measurements of errors.
 
main_v5_pdf = function(### Parameters for generative model
                        SEED, 
                        N_node = 100,
                        N_clus=2, 
                        conn_patt_sep = 4,
                        u_1 = 1, u_0 = 1,
                        t_vec = seq(-u_0,u_1,length.out=200),
                        ### params when N_clus==1:
                        N_spks_total = 60*10+40*10,
                        N_spks_ratio = 3/2,
                        sd_shrinkage = 1,
                        ### params when N_clus==2:
                        clus_mixture = 0,
                        ### Parameters for algorithms
                        N_component=2,
                        gamma=0,
                        MaxIter=10,
                        N_clus_min=N_clus, N_clus_max=N_clus,
                        fix_timeshift=FALSE,
                        use_true_timeshift=FALSE,
                        jitter_prop_true_timeshift=0,
                        save_center_pdf_array=FALSE,
                        ### Unused
                        jitter_time_rad = 10, max_iter=50,
                        conv_thres=5e-3, 
                        freq_trun_vec=NULL,
                        opt_radius=total_time/2,
                        ...)
{
  
  # Generate data -------------------------------------------------------
  
  network_list = generate_data(SEED=SEED,
                               N_node=N_node,
                               N_clus=N_clus, 
                               u_1=u_1, u_0=u_0,
                               t_vec=t_vec,
                               N_spks_total = N_spks_total,
                               N_spks_ratio = N_spks_ratio,
                               sd_shrinkage = sd_shrinkage,
                               clus_mixture = clus_mixture)
  
  spks_time_mlist = network_list$spks_time_mlist
  stim_onset_vec = network_list$stim_onset_vec
  
  center_density_array_true = network_list$center_density_array_true
  center_intensity_array_true = network_list$center_intensity_array_true
  mem_true_vec = network_list$mem_true_vec
  clus_true_list = network_list$clus_true_list
  v_true_list = network_list$v_vec_list
  t_vec = network_list$t_vec
  

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
                   v0 = u_1, v1 = u_0,
                   t_vec = t_vec, 
                   fix_timeshift = fix_timeshift, 
                   use_true_timeshift = use_true_timeshift, 
                   jitter_prop_true_timeshift = jitter_prop_true_timeshift, 
                   v_true_list = v_true_list)
    clusters_list_init = res$clusters_list
    v_vec_list_init = res$v_vec_list
    
    
    # Apply algorithm ---------
    
    time_start = Sys.time()
    ### Estimation z,v,f based on pdf
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                         stim_onset_vec = stim_onset_vec,
                         clusters_list_init = clusters_list_init,
                         v_vec_list_init = v_vec_list_init,
                         N_component = N_component, 
                         freq_trun = Inf,
                         MaxIter = MaxIter, 
                         gamma=gamma,
                         v0 = u_1, v1= u_0,
                         t_vec=t_vec, 
                         fix_timeshift = fix_timeshift,
                         conv_thres = conv_thres,
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
    if(N_clus>=2){
      res = find_permn(center_density_array_from = center_density_array_est,
                       center_density_array_to = center_density_array_true)
      permn = res$permn
    } else{
      permn = 1
    }
    
    center_density_array_est_permn = center_density_array_est[permn, , ,drop=FALSE]
    center_intensity_array_est_permn = center_intensity_array_est[permn, , ,drop=FALSE]
    center_Nspks_mat_est_permn = center_Nspks_mat_est[permn, ,drop=FALSE]
    
    
    ### Calculate distance 
    dist_mse_mat = matrix(nrow=N_clus, ncol=N_component)
    for (id_clus in 1:N_clus) {
      for (id_component in 1:N_component) {
        dist_mse_mat[id_clus,id_component] = sum( (center_density_array_est_permn[id_clus,id_component,] - 
                                                     center_density_array_true[id_clus,id_component,])^2 * 
                                                    (t_vec[2]-t_vec[1]) )
      }
    }
    F_mean_sq_err = mean(dist_mse_mat)
    F_mean_sq_err_vec = colMeans(dist_mse_mat)
    
    F_l2_squared_norm_mat = apply(center_density_array_true, 1:2, function(density){
      sum(density^2 * (t_vec[2]-t_vec[1]))
    })
    F_mse_squarel2_ratio_vec = colMeans( dist_mse_mat / F_l2_squared_norm_mat )
    
    
    # Compute errors of clusters, i.e. Z ------------------------------------
    clusters_list_est_permn = clusters_list_est[permn]
    ARI = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est), 
                      memb_true_vec = mem_true_vec)
    

    # Compute error of time shifts, i.e. w, v --------------------------------------------
    v_mean_sq_err = mean((unlist(v_true_list)-unlist(v_vec_list_est))^2) 
    v_mean_sq_err_vec = sapply(1:N_component, function(id_component){
      mean((unlist(v_true_list[[id_component]])-unlist(v_vec_list_est[[id_component]]))^2) 
    })
      
    
    
    
  }
  else{
    ARI=NA 
    F_mean_sq_err=F_mean_sq_err_vec=F_mse_squarel2_ratio_vec=NA 
    v_mean_sq_err=v_mean_sq_err_vec=NA
    clusters_list_est_permn=NA
    center_density_array_est_permn=NA
    center_intensity_array_est_permn=NA
    center_Nspks_mat_est_permn=NA
  }
  
  # Extract network related parameters -----------------------------------------
  data_param = list(SEED=SEED,
                    N_node=N_node,
                    N_clus=N_clus, 
                    u_1=u_1, u_0=u_0,
                    t_vec=t_vec,
                    N_spks_total = N_spks_total,
                    N_spks_ratio = N_spks_ratio,
                    sd_shrinkage = sd_shrinkage,
                    clus_mixture = clus_mixture)
  
  
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
              center_intensity_array_est_permn=switch(save_center_pdf_array, 
                                                      "TRUE"=center_intensity_array_est_permn,
                                                      "FALSE"=NULL),
              center_Nspks_mat_est_permn=switch(save_center_pdf_array, 
                                                "TRUE"=center_Nspks_mat_est_permn,
                                                "FALSE"=NULL),
              # estimation error
              ARI=ARI, 
              F_mean_sq_err=F_mean_sq_err, 
              F_mean_sq_err_vec=F_mean_sq_err_vec, 
              F_mse_squarel2_ratio_vec=F_mse_squarel2_ratio_vec,
              v_mean_sq_err=v_mean_sq_err,
              v_mean_sq_err_vec=v_mean_sq_err_vec,
              # other
              time_estimation=time_estimation,
              N_iteration=N_iteration,
              loss_history=loss_history
              ))
  
}


