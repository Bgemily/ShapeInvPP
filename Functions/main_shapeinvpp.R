### Generate data, run our algorithm, and output measurements of errors.
 
main_shapeinvpp = function(### Parameters for generative model
                        SEED, 
                        N_subj = 100,
                        N_trial = 1,
                        N_clus=2, 
                        N_component_true = 2,
                        t_vec = seq(-1,1,by=0.01),
                        N_spks_total = 1000,
                        timeshift_subj_max_vec = c(1/8, 1/32),
                        timeshift_trial_max = 1/8,
                        ### params when N_clus==4:
                        clus_sep = 2,
                        ### Parameters for algorithms
                        freq_trun = 5,
                        bw = 0,
                        N_component = 2,
                        v_subjwise_max = NULL,
                        key_times_vec = c(min(t_vec),0,max(t_vec)),
                        gamma = 0,
                        MaxIter = 10,
                        N_clus_est = N_clus,
                        fix_timeshift = FALSE,
                        use_true_timeshift = FALSE,
                        save_center_pdf_array = FALSE,
                        rand_init = FALSE,
                        N_restart = 1,
                        N_start_kmean = 5,
                        ### Unused
                        conv_thres = 5e-3, 
                        ...)
{
  t_unit = t_vec[2]-t_vec[1]
  # Generate data -------------------------------------------------------
  ### Extract network related parameters 
  data_param = list(SEED=SEED,
                    N_subj=N_subj,
                    N_trial=N_trial,
                    N_clus=N_clus, 
                    t_vec=t_vec,
                    key_times_vec = key_times_vec,
                    N_spks_total = N_spks_total,
                    timeshift_subj_max_vec = timeshift_subj_max_vec,
                    timeshift_trial_max = timeshift_trial_max,
                    clus_sep = clus_sep )
  data_generated = do.call(what = generate_data, args = data_param)
  
  
  spks_time_mlist = data_generated$spks_time_mlist

  center_density_array_true = data_generated$center_density_array_true
  center_intensity_array_true = data_generated$center_intensity_array_true
  mem_true_vec = data_generated$mem_true_vec
  clus_true_list = data_generated$clus_true_list
  v_true_mat_list = data_generated$v_mat_list
  v_trialwise_vec_list = data_generated$v_trialwise_vec_list
  
  # Calculate non-identifiability level --------
  if (N_component >= 2) {
    non_identifiability = get_non_identifiability(spks_time_mlist = spks_time_mlist, 
                                                    v_mat_list = v_true_mat_list, 
                                                    N_component = N_component, 
                                                    t_vec = t_vec, 
                                                  freq_trun = freq_trun)
  } else {
    non_identifiability = NA
  }
  
  # Fit model for various cluster number ------------------------------------
  
  res_list = list()
  for (ind_N_clus in 1:1) {
    res_list[[ind_N_clus]] = list()
    N_clus_tmp = N_clus_est
    
    # Restart -----------------------------------------------------------------
    res_best = NA
    l2_loss_best = Inf
    l2_loss_history = c()
    for (id_restart in 1:N_restart) {
      ### Get initialization -----------
      if (rand_init) {
        res = get_init_random(spks_time_mlist = spks_time_mlist, 
                              N_clus = N_clus_tmp,
                              N_component = N_component,
                              t_vec = t_vec, 
                              key_times_vec = key_times_vec,
                              N_start_kmean = N_start_kmean,
                              freq_trun = freq_trun,
                              bw = bw,
                              fix_timeshift = fix_timeshift, 
                              use_true_timeshift = use_true_timeshift, 
                              add_rand_to_init_timeshift = ifelse(id_restart>1, TRUE, FALSE),
                              v_true_mat_list = v_true_mat_list, 
                              v_trialwise_vec_list = v_trialwise_vec_list,
                              rmv_conn_prob = TRUE)
        center_density_array_init = res$center_density_array
        center_Nspks_mat_init = res$center_Nspks_mat
        clusters_list_init = res$clusters_list
        v_mat_list_init = res$v_mat_list
      } else {
        res = get_init(spks_time_mlist = spks_time_mlist, 
                       N_clus = N_clus_tmp,
                       N_component = N_component,
                       t_vec = t_vec, 
                       key_times_vec = key_times_vec,
                       N_start_kmean = N_start_kmean,
                       freq_trun = freq_trun,
                       bw = bw,
                       fix_timeshift = fix_timeshift, 
                       use_true_timeshift = use_true_timeshift, 
                       add_rand_to_init_timeshift = ifelse(id_restart>1, TRUE, FALSE),
                       v_true_mat_list = v_true_mat_list, 
                       v_trialwise_vec_list = v_trialwise_vec_list,
                       rmv_conn_prob = TRUE)
        center_density_array_init = res$center_density_array
        center_Nspks_mat_init = res$center_Nspks_mat
        clusters_list_init = res$clusters_list
        v_mat_list_init = res$v_mat_list
      }
      
      
      # Apply algorithm ---------
      time_start = Sys.time()
      ### Estimation z,v,f based on pdf
      res_new = apply_asimm(spks_time_mlist = spks_time_mlist,
                               v_trialwise_vec_list = v_trialwise_vec_list,
                               clusters_list_init = clusters_list_init,
                               v_mat_list_init = v_mat_list_init,
                               N_component = N_component, 
                               freq_trun = freq_trun,
                               bw = bw,
                               MaxIter = MaxIter, 
                               gamma=gamma,
                               t_vec=t_vec, 
                               v_subjwise_max = v_subjwise_max,
                               key_times_vec = key_times_vec,
                               fix_timeshift = fix_timeshift, 
                               rand_init = rand_init,
                               conv_thres = conv_thres )
      time_end = Sys.time()
      time_estimation = time_end - time_start
      time_estimation = as.numeric(time_estimation, units='secs')
      
      ### Extract loss function value
      l2_loss_new = tail(res_new$loss_history, 2)[1]
      
      ### Update best estimation
      if(l2_loss_new < l2_loss_best){
        l2_loss_best = l2_loss_new
        res_best = res_new
      }
      l2_loss_history[id_restart] = l2_loss_best
      
    }
    
    # Save results of N_clus_tmp ----------------------------------------------
    res_list[[ind_N_clus]] = res_best
  }
  
  # Retrieve estimation results of the best cluster number ------------------
  res = res_list[[1]]
  
  res$clusters_list -> clusters_list_est
  res$clusters_history -> clusters_history
  res$center_density_array_history -> center_density_array_history
  res$v_mat_list -> v_mat_list_est
  res$center_intensity_array -> center_intensity_array_est
  res$center_density_array -> center_density_array_est
  res$center_Nspks_mat -> center_Nspks_mat_est
  res$loss_history -> loss_history
  res$N_iteration -> N_iteration
  
  # Evaluate L2 loss and log-likelihood for the results of the best cluster number -----
  tmp = evaluate_model(spks_time_mlist = spks_time_mlist, 
                       v_trialwise_vec_list = v_trialwise_vec_list, 
                       N_component = N_component, 
                       key_times_vec = key_times_vec, 
                       model_fitted_list = res,
                       freq_trun = freq_trun)
  log_lik = tmp$log_lik
  log_lik_tmp_1 = tmp$log_lik_tmp_1
  log_lik_tmp_2 = tmp$log_lik_tmp_2
  clus_entropy = tmp$clus_entropy
  L2_loss_part_1 = tmp$L2_loss_part_1
  L2_loss_part_2 = tmp$L2_loss_part_2
  L2_loss_part_1_smoothdensity = tmp$L2_loss_part_1_smoothdensity
  compl_log_lik = tmp$compl_log_lik
  
  # Compute estimation error ------------------------------------------------
  # Compute errors of clusters, i.e. Z ------------------------------------
  ARI = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est), 
                    memb_true_vec = mem_true_vec)
  ARI_history = sapply(clusters_history, function(clusters_list_tmp)get_one_ARI(memb_est_vec = clus2mem(clusters_list_tmp), 
                                                                                memb_true_vec = mem_true_vec))
  
  # Compute error of time shifts, i.e. v --------------------------------------------
  if (N_component == N_component_true) {
    v_mean_sq_err_vec = sapply(1:N_component, function(id_component){
      mean((unlist(v_true_mat_list[[id_component]])-unlist(v_mat_list_est[[id_component]]))^2) 
    })
    v_mean_sq_err = mean(v_mean_sq_err_vec)
    
    v_align_mean_sq_err_vec = c()
    v_align_mat_list_est = v_mat_list_est
    for (id_component in 1:N_component) {
      for (id_clus in 1:N_clus_est) {
        v_align_mat_list_est[[id_component]][clusters_list_est[[id_clus]], ] = 
          v_mat_list_est[[id_component]][clusters_list_est[[id_clus]], ] - 
          mean(v_mat_list_est[[id_component]][clusters_list_est[[id_clus]], ]) + 
          mean(v_true_mat_list[[id_component]][clusters_list_est[[id_clus]], ])
      }
      mse_tmp = mean(( unlist(v_true_mat_list[[id_component]])-unlist(v_align_mat_list_est[[id_component]]) )^2) 
      v_align_mean_sq_err_vec[id_component] = mse_tmp
    }
    v_align_mean_sq_err = mean(v_align_mean_sq_err_vec)
    
  } else {
    v_mean_sq_err = v_mean_sq_err_vec = NA
    v_align_mean_sq_err_vec = v_align_mean_sq_err = NA
  }
  
  # Compute errors of conn patts, i.e. F ---------
  if (N_clus_est == N_clus) {
    ### Match clusters 
    if(N_clus>=2){
      permn = find_permn_clus_label(clusters_list_true = data_generated$clus_true_list, 
                                        clusters_list_est = clusters_list_est)
    } else{
      permn = 1
    }
    
    center_density_array_est_permn = center_density_array_est[permn, , ,drop=FALSE]
    center_intensity_array_est_permn = center_intensity_array_est[permn, , ,drop=FALSE]
    center_Nspks_mat_est_permn = center_Nspks_mat_est[permn, ,drop=FALSE]
    clusters_list_est_permn = clusters_list_est[permn]
    
    
    ### Calculate distance 
    calculate_F_mse_ratio = function(N_clus, N_component, 
                                     clusters_list_true, clusters_list_est, 
                                     center_density_array_true, center_density_array_est, 
                                     t_vec, t_unit) 
    {
      permn = find_permn_clus_label(clusters_list_true = clusters_list_true, 
                                    clusters_list_est = clusters_list_est)
      center_density_array_est_permn = center_density_array_est[permn, , ,drop=FALSE]
      clusters_list_est_permn = clusters_list_est[permn]
      
      dist_mse_mat = matrix(nrow=N_clus, ncol=N_component)
      for (id_clus in 1:N_clus) {
        for (id_component in 1:N_component) {
          f_target = center_density_array_true[id_clus,id_component,]
          density_est = center_density_array_est_permn[id_clus,id_component,]
          res_ccf = ccf(y = density_est, x = f_target, plot = FALSE, lag.max = length(t_vec)%/%2)
          n0_init = res_ccf$lag[which.max(res_ccf$acf)]
          if (length(n0_init) == 0) {
            n0_init = 0
          }
          f_target_array = array(data = f_target, dim = c(1,1,length(f_target)))
          f_origin_mat = matrix(density_est, nrow = 1)
          n0 = align_multi_components(f_target_array = f_target_array,
                                      f_origin_mat = f_origin_mat,
                                      n0_init_mat = as.matrix(n0_init),
                                      v_trialwise_vec_list = list(c(0)),
                                      N_spks_mat = as.matrix(c(1)),
                                      t_unit = t_unit, 
                                      n0_min_vec = -length(f_target) %/% 2,
                                      n0_max_vec = length(f_target) %/% 2 )$n0_mat
          n0 = round(n0) 
          if (n0 > 0) {
            if (FALSE) {
              density_est_shift = c(rep(0, n0), head(density_est, length(density_est) - n0) )
            } else {
              density_est_shift = c(rep(head(density_est,1), n0), head(density_est, length(density_est) - n0) )
            }
            
          } else if (n0 < 0) {
            if (FALSE) {
              density_est_shift = c(tail(density_est, length(density_est) - abs(n0)), rep(0, abs(n0)) )
            } else {
              density_est_shift = c(tail(density_est, length(density_est) - abs(n0)), rep(tail(density_est,1), abs(n0)) )
            }
            
          } else if (n0 == 0) {
            density_est_shift = density_est
          }
          dist_mse_mat[id_clus,id_component] = sum( (density_est_shift - 
                                                       center_density_array_true[id_clus,id_component,])^2 * 
                                                      t_unit )
        }
      }
      weight_vec = sapply(clusters_list_est_permn, length) / length(unlist(clusters_list_est_permn))
      F_mean_sq_err = sum( rowMeans(dist_mse_mat) * weight_vec )
      F_mean_sq_err_vec = colMeans(dist_mse_mat)
      F_l2_squared_norm_mat = apply(center_density_array_true, 1:2, function(density){
        sum(density^2 * t_unit)
      })
      F_mse_squarel2_ratio_mat =  dist_mse_mat / F_l2_squared_norm_mat 
      F_mse_squarel2_ratio = sum( rowMeans(F_mse_squarel2_ratio_mat) * weight_vec )
      F_l2_squared_norm_mat_2 = apply(center_density_array_est_permn, 1:2, function(density){
        sum(density^2 * t_unit)
      })
      F_mse_squarel2_ratio_mat_2 =  dist_mse_mat / F_l2_squared_norm_mat_2 
      F_mse_squarel2_ratio_2 = sum( rowMeans(F_mse_squarel2_ratio_mat_2) * weight_vec )
      
      return(list(dist_mse_mat = dist_mse_mat,
                  F_mean_sq_err = F_mean_sq_err,
                  F_mean_sq_err_vec = F_mean_sq_err_vec,
                  F_mse_squarel2_ratio = F_mse_squarel2_ratio,
                  F_mse_squarel2_ratio_mat = F_mse_squarel2_ratio_mat,
                  F_mse_squarel2_ratio_2 = F_mse_squarel2_ratio_2,
                  F_mse_squarel2_ratio_mat_2 = F_mse_squarel2_ratio_mat_2))
    }
    
    center_density_array_true_rmv_baseline = center_density_array_true
    for (id_clus in 1:N_clus) {
      id_component = 1
      baseline = head(center_density_array_true[id_clus, id_component, ],1)
      center_density_array_true_rmv_baseline[id_clus, id_component, ] = center_density_array_true[id_clus, id_component, ] - baseline
    }
    
    N_component_true = dim(center_density_array_true)[2]
    if (N_component == N_component_true) {
      tmp = calculate_F_mse_ratio(N_clus = N_clus, N_component = N_component, 
                                  clusters_list_true = data_generated$clus_true_list,
                                  clusters_list_est = clusters_list_est,
                                  center_density_array_true = center_density_array_true_rmv_baseline, 
                                  center_density_array_est = center_density_array_est, 
                                  t_vec = t_vec, 
                                  t_unit = t_unit)
      dist_mse_mat = tmp$dist_mse_mat
      F_mean_sq_err = tmp$F_mean_sq_err
      F_mean_sq_err_vec = tmp$F_mean_sq_err_vec
      F_mse_squarel2_ratio = tmp$F_mse_squarel2_ratio
      F_mse_squarel2_ratio_mat = tmp$F_mse_squarel2_ratio_mat
      F_mse_squarel2_ratio_2 = tmp$F_mse_squarel2_ratio_2
      F_mse_squarel2_ratio_mat_2 = tmp$F_mse_squarel2_ratio_mat_2
      
      F_mean_sq_err_history = c()
      F_mse_squarel2_ratio_history = c()
      F_mse_squarel2_ratio_history_2 = c()
      for (id_iteration in 1:length(center_density_array_history)){
        tmp = calculate_F_mse_ratio(N_clus = N_clus, N_component = N_component, 
                                    clusters_list_true = data_generated$clus_true_list,
                                    clusters_list_est = clusters_history[[id_iteration]],
                                    center_density_array_true = center_density_array_true_rmv_baseline, 
                                    center_density_array_est = center_density_array_history[[id_iteration]], 
                                    t_vec = t_vec, 
                                    t_unit = t_unit)
        F_mean_sq_err_history[id_iteration] = tmp$F_mean_sq_err
        F_mse_squarel2_ratio_history[id_iteration] = tmp$F_mse_squarel2_ratio
        F_mse_squarel2_ratio_history_2[id_iteration] = tmp$F_mse_squarel2_ratio_2
      }
      
    } else {
      F_mean_sq_err = F_mean_sq_err_vec = F_mse_squarel2_ratio = NA
      F_mse_squarel2_ratio_mat = dist_mse_mat = NA
      F_mean_sq_err_history = NA
      F_mse_squarel2_ratio_history = NA
      F_mse_squarel2_ratio_2 = NA
      F_mse_squarel2_ratio_mat_2 = NA
      F_mse_squarel2_ratio_history_2 = NA
    }
    
    
    
  }
  else{
    F_mean_sq_err=F_mean_sq_err_vec=F_mse_squarel2_ratio=NA
    F_mse_squarel2_ratio_mat = dist_mse_mat = NA
    clusters_list_est_permn=NA
    center_density_array_est_permn=NA
    center_intensity_array_est_permn=NA
    center_Nspks_mat_est_permn=NA
    F_mean_sq_err_history = NA
    F_mse_squarel2_ratio_history = NA
    F_mse_squarel2_ratio_2 = NA
    F_mse_squarel2_ratio_mat_2 = NA
    F_mse_squarel2_ratio_history_2 = NA
  }
  
  
  # Output ------------------------------------------------------------------
  
  return(list(data_param=data_param, 
              N_clus_est=N_clus_est, 
              correct_N_clus=I(N_clus_est==N_clus)*1, 
              clusters_list_est=clusters_list_est,
              v_mat_list_est=v_mat_list_est,
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
              # generated data
              data_generated=switch(save_center_pdf_array, 
                                  "TRUE"=data_generated,
                                  "FALSE"=NULL),
              non_identifiability = non_identifiability,
              # estimation error
              ARI=ARI, 
              F_mean_sq_err=F_mean_sq_err, 
              F_mean_sq_err_vec=F_mean_sq_err_vec, 
              F_mse_squarel2_ratio=F_mse_squarel2_ratio,
              F_mse_squarel2_ratio_mat =  F_mse_squarel2_ratio_mat,
              dist_mse_mat = dist_mse_mat,
              v_mean_sq_err=v_mean_sq_err,
              v_mean_sq_err_vec=v_mean_sq_err_vec,
              v_align_mean_sq_err = v_align_mean_sq_err,
              ARI_history = ARI_history,
              F_mean_sq_err_history = F_mean_sq_err_history,
              F_mse_squarel2_ratio_history = F_mse_squarel2_ratio_history,
              F_mse_squarel2_ratio_2 = F_mse_squarel2_ratio_2,
              F_mse_squarel2_ratio_mat_2 = F_mse_squarel2_ratio_mat_2,
              F_mse_squarel2_ratio_history_2 = F_mse_squarel2_ratio_history_2,
              # L2 loss and log-likelihood
              log_lik = log_lik,
              log_lik_tmp_1 = log_lik_tmp_1,
              log_lik_tmp_2 = log_lik_tmp_2,
              clus_entropy = clus_entropy,
              L2_loss_part_1 = L2_loss_part_1,
              L2_loss_part_2 = L2_loss_part_2,
              L2_loss_part_1_smoothdensity = L2_loss_part_1_smoothdensity,
              compl_log_lik = compl_log_lik,
              # other
              N_restart = N_restart,
              t_vec=t_vec, 
              time_estimation=time_estimation,
              N_iteration=N_iteration,
              loss_history=loss_history
              ))
  
}


