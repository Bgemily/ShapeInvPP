### Generate data, apply kcfc, and output measurements of errors.
# library(fdapace)

main_kcfc = function(### Parameters for generative model
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
  
  # Prepare data for FPCA ######
  yList = list()
  tList = list()
  for (id_node in 1:N_node){
    res_smooth = density(spks_time_mlist[[id_node]], bw = bw, n = 256, from = min(t_vec), to = max(t_vec))
    yList[[id_node]] = res_smooth$y
    tList[[id_node]] = res_smooth$x
  }
  
  
  ### Apply kCFC --------------------------------
  kcfcObj = kCFC(y = yList, t = tList, k = N_clus, 
                 kSeed = 123, maxIter = 20,
                 optnsSW = list(dataType='Dense', maxK=N_component, FVEthreshold = 1), 
                 optnsCS = list(dataType='Dense', maxK=N_component, FVEthreshold = 1))
  
  # Find the best permutation of cluster labels
  clusters_list_true = data_generated$clus_true_list
  clusters_list_est = mem2clus(as.numeric(kcfcObj$cluster))
  the_permn = find_permn_clus_label(clusters_list_true = clusters_list_true, 
                                    clusters_list_est = clusters_list_est)
  memb_est_vec_permn = clus2mem(clusters_list_est[the_permn])
  fpcaList_permn = kcfcObj$fpcaList[the_permn]    
  
  # Reconstruct density components and time shifts 
  calculate_deriv = function(x_vec,f_vec){
    n = length(f_vec)
    grad = rep(0, n)
    
    if (length(f_vec) > 1){
      grad[1] = (f_vec[2] - f_vec[1]) / (x_vec[2]-x_vec[1])
      grad[n] = (f_vec[n] - f_vec[n-1]) / (x_vec[n]-x_vec[n-1])
    }
    
    if (length(f_vec) > 2){
      grad[2:(n-1)] = (f_vec[3:n]-f_vec[1:(n-2)]) / (x_vec[3:n] - x_vec[1:(n-2)])
    }
    
    return(grad)
  }
  calculate_antideriv = function(x_vec, gradient_vec){
    dx_vec = c(0, diff(x_vec))
    f_vec = cumsum(gradient_vec*dx_vec)
    return(f_vec)
  }
  
  center_density_fpca_array_permn = array(dim = c(N_clus, N_component, length(t_vec_extend)))
  v_fpca_mat = matrix(nrow = N_node, ncol = N_component)
  for (id_clus in 1:N_clus){
    ### Extract FPCA estimates
    FPCAobj = fpcaList_permn[[id_clus]]    
    N_node_tmp = nrow(FPCAobj$xiEst)
    mean_density_vec = FPCAobj$mu
    eigenfuncs_mat = FPCAobj$phi[ , 1:N_component, drop=FALSE] # len(t_fpca_vec) x N_component
    fpc_scores_mat = FPCAobj$xiEst[ , 1:N_component, drop=FALSE] # N_node_tmp x N_component
    t_fpca_vec = FPCAobj$workGrid
    current_cluster = which(memb_est_vec_permn==id_clus)
    
    for (id_component in 1:N_component){
      derivative_mean_density_vec = calculate_deriv(x_vec = t_fpca_vec, f_vec = mean_density_vec)
      eigen_func_vec = eigenfuncs_mat[,id_component]
      inner_prod = sum( derivative_mean_density_vec * eigen_func_vec ) * (t_fpca_vec[2]-t_fpca_vec[1])
      derivative_density_vec = eigen_func_vec * inner_prod
      density_vec = calculate_antideriv(x_vec = t_fpca_vec, gradient_vec = derivative_density_vec)
      density_vec = approx(x = t_fpca_vec, xout = data_generated$t_vec_extend,
                           y = density_vec)$y
      center_density_fpca_array_permn[id_clus, id_component, ] = density_vec
      
      v_fpca_vec = -fpc_scores_mat[, id_component] * (inner_prod^(-1))
      v_fpca_vec = v_fpca_vec - min(v_fpca_vec)
      v_fpca_mat[current_cluster, id_component] = v_fpca_vec
    }
  }
  
  
  # Find the best permutation of density components
  v_fpca_mat_permn = 0 * v_fpca_mat
  for (id_clus in 1:N_clus){
    current_cluster = which(memb_est_vec_permn==id_clus)
    permn_list = combinat::permn(1:N_component)  
    mise_f_max = Inf
    the_permn = c()
    for (permn in permn_list) {
      mise_f_tmp_vec = c()
      for (id_component in 1:N_component) {
        t_unit = data_generated$t_vec[2] - data_generated$t_vec[1]
        mise_f_tmp_vec[id_component] = sum( t_unit * (center_density_fpca_array_permn[id_clus, permn[id_component], ] - 
                                                        data_generated$center_density_array_true[id_clus, id_component, ])^2 )
      }
      mise_f_tmp = mean(mise_f_tmp_vec)
      if (mise_f_tmp < mise_f_max){
        the_permn = permn
        mise_f_max = mise_f_tmp
      }
    }
    center_density_fpca_array_permn[id_clus, , ] = center_density_fpca_array_permn[id_clus, the_permn, ]
    v_fpca_mat_permn[current_cluster, ] = v_fpca_mat[current_cluster, the_permn]
  }

  # Rename estimates
  clusters_list_est = clusters_list_est
  clusters_list_est_permn = mem2clus(membership = memb_est_vec_permn, N_clus_min = N_clus)
  center_density_array_est_permn = center_density_fpca_array_permn
  v_mat_list_est = list()
  for (id_component in 1:N_component) {
    v_mat_list_est[[id_component]] = matrix(data = v_fpca_mat_permn[, id_component], nrow = N_node, ncol = N_replicate)
  }
  
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
    

    # Compute errors of clusters, i.e. Z ------------------------------------
    ARI = get_one_ARI(memb_est_vec = clus2mem(clusters_list_est), 
                      memb_true_vec = mem_true_vec)
    
    
    # Compute error of time shifts, i.e. w, v --------------------------------------------
    v_mean_sq_err_vec = sapply(1:N_component, function(id_component){
      mean((unlist(v_true_mat_list[[id_component]])-unlist(v_mat_list_est[[id_component]]))^2) /
        ( 2*timeshift_max_vec[id_component] )^2
    })
    v_mean_sq_err = mean(v_mean_sq_err_vec)
    
    
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
              time_estimation=NA,
              N_iteration=NA,
              loss_history=NA
  ))
  
}


