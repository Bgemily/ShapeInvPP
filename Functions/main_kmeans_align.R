### Generate data, apply kmeans_align(), and output measurements of errors.
# library(fdasrvf)

main_kmeans_align = function(### Parameters for generative model
  SEED, 
  N_node = 100,
  N_replicate = 1,
  N_clus=2, 
  N_component_true = 2,
  u_1 = 1, u_0 = 1,
  t_vec = seq(-u_0,u_1,by=0.01),
  t_vec_extend = t_vec,
  N_spks_total = 1000,
  timeshift_max_vec = c(1/8, 1/32),
  ### params when N_clus==4:
  clus_sep = 2,
  ### params when N_clus==1:
  N_spks_ratio = 3/2,
  sd_shrinkage = 1,
  c_1 = 0, delta_1 = 0,
  c_2 = 0, delta_2 = 0,
  c_3 = 0, delta_3 = 0,
  identical_components = FALSE,
  ### params when N_clus==2:
  clus_mixture = 0,
  ### Parameters for algorithms
  bw = 0,
  N_component = 2,
  save_center_pdf_array = FALSE)
{
  t_unit = t_vec[2]-t_vec[1]
  # Generate data -------------------------------------------------------
  ### Extract network related parameters 
  data_param = list(SEED=SEED,
                    N_node=N_node,
                    N_replicate=N_replicate,
                    N_clus=N_clus, 
                    u_1=u_1, u_0=u_0,
                    t_vec=t_vec,
                    t_vec_extend=t_vec_extend,
                    N_spks_total = N_spks_total,
                    timeshift_max_vec = timeshift_max_vec,
                    clus_sep = clus_sep,
                    N_spks_ratio = N_spks_ratio,
                    sd_shrinkage = sd_shrinkage,
                    c_1 = c_1, delta_1 = delta_1,
                    c_2 = c_2, delta_2 = delta_2,
                    c_3 = c_3, delta_3 = delta_3,
                    identical_components = identical_components,
                    clus_mixture = clus_mixture)
  
  if (N_component_true == 1) {
    data_generated = do.call(what = generate_data_Ncomp_1, args = data_param)
  } else if (N_component_true == 2) {
    data_generated = do.call(what = generate_data, args = data_param)
  }
  
  spks_time_mlist = data_generated$spks_time_mlist
  stim_onset_vec = data_generated$stim_onset_vec
  
  center_density_array_true = data_generated$center_density_array_true
  center_intensity_array_true = data_generated$center_intensity_array_true
  mem_true_vec = data_generated$mem_true_vec
  clus_true_list = data_generated$clus_true_list
  v_true_mat_list = data_generated$v_mat_list
  
  # Prepare data for kmeans_align() ######
  f_mat = c()
  time_vec = c()
  for (id_node in 1:N_node){
    res_smooth = density(spks_time_mlist[[id_node]], bw = bw, n = length(t_vec), from = min(t_vec), to = max(t_vec))
    f_mat = cbind(f_mat, res_smooth$y)
    time_vec = res_smooth$x
  }
  
  
  ### Apply kmeans_align() --------------------------------
  set.seed(SEED)
  time_start = Sys.time()
  fdakma_obj = fdasrvf::kmeans_align(f = f_mat,
                                     time = time_vec,
                                     K = N_clus,
                                     alignment = TRUE, 
                                    nonempty = 2, 
                                     showplot = FALSE)
  time_end = Sys.time()
  time_estimation = time_end - time_start
  time_estimation = as.numeric(time_estimation, units='secs')
  
  # Find the best permutation of cluster labels
  clusters_list_true = data_generated$clus_true_list
  clusters_list_est = mem2clus(as.numeric(fdakma_obj$labels))
  the_permn = find_permn_clus_label(clusters_list_true = clusters_list_true, 
                                    clusters_list_est = clusters_list_est)
  clusters_list_est_permn = clusters_list_est[the_permn]

  # Extract center densities 
  center_density_mat_est = t(fdakma_obj$templates)
  center_density_mat_est_permn = center_density_mat_est[the_permn, , drop = FALSE]
  center_density_array_est_permn = array(dim = c(N_clus, 1, length(t_vec)))
  center_density_array_est_permn[ , 1, ] = center_density_mat_est_permn
  
  # Rename estimates
  clusters_list_est = clusters_list_est
  clusters_list_est_permn = clusters_list_est_permn
  center_density_array_est_permn = center_density_array_est_permn
  v_mat_list_est = NA
  
  # Other estimates
  N_clus_est = N_clus
  center_intensity_array_est = NA
  center_Nspks_mat_est = NA
  loss_history = NA
  N_iteration = NA
  center_intensity_array_est_permn = NA
  center_Nspks_mat_est_permn = NA
  
  
  # Compute estimation error ------------------------------------------------
  if (TRUE) {
    # Compute errors of conn patts, i.e. F ---------
    if (N_component_true == N_component) {
      ### Calculate distance 
      dist_mse_mat = matrix(nrow=N_clus, ncol=N_component)
      for (id_clus in 1:N_clus) {
        for (id_component in 1:N_component) {
          dist_mse_mat[id_clus,id_component] = sum( (center_density_array_est_permn[id_clus,id_component,] - 
                                                       center_density_array_true[id_clus,id_component,])^2 * 
                                                      t_unit )
        }
      }
      F_mean_sq_err = mean(dist_mse_mat)
      F_mean_sq_err_vec = colMeans(dist_mse_mat)
      
      F_l2_squared_norm_mat = apply(center_density_array_true, 1:2, function(density){
        sum(density^2 * t_unit)
      })
      F_mse_squarel2_ratio_mat =  dist_mse_mat / F_l2_squared_norm_mat 
      
      weight_vec = sapply(clusters_list_est_permn, length) / length(unlist(clusters_list_est_permn))
      F_mse_squarel2_ratio = sum( rowMeans(F_mse_squarel2_ratio_mat) * weight_vec )
    } else {
      dist_mse_mat = NA
      F_mean_sq_err = NA
      F_mean_sq_err_vec = NA
      F_mse_squarel2_ratio_mat = NA
      F_mse_squarel2_ratio = NA
    }
    
    
    # Compute errors of clusters, i.e. Z ------------------------------------
    ARI = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est), 
                      memb_true_vec = mem_true_vec)
    
    
    # Compute error of time shifts, i.e. w, v --------------------------------------------
    v_mean_sq_err_vec = NA
    v_mean_sq_err = NA
  }
  
  
  # Output ------------------------------------------------------------------
  
  return(list(data_param=data_param, 
              # model selection result
              N_clus_est=N_clus_est, 
              correct_N_clus=I(N_clus_est==N_clus)*1, 
              ICL_vec=NA, 
              compl_log_lik_vec=NA, 
              log_lik_vec=NA,
              penalty_vec=NA,
              # parameter estimates of best cluster number
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
              # estimation error
              ARI=ARI, 
              F_mean_sq_err=F_mean_sq_err, 
              F_mean_sq_err_vec=F_mean_sq_err_vec, 
              F_mse_squarel2_ratio=F_mse_squarel2_ratio,
              F_mse_squarel2_ratio_mat =  F_mse_squarel2_ratio_mat,
              dist_mse_mat = dist_mse_mat,
              v_mean_sq_err=v_mean_sq_err,
              v_mean_sq_err_vec=v_mean_sq_err_vec,
              # other
              cand_N_clus_vec=NA,
              N_restart = NA,
              t_vec=t_vec, t_vec_extend=t_vec_extend,
              time_estimation = time_estimation,
              N_iteration=NA,
              loss_history=NA
  ))
  
}


