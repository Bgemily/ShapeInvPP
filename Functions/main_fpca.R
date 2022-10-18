### Generate data, apply fpca, and output measurements of errors.
# library(fdapace)

main_fpca = function(### Parameters for generative model
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
    res_smooth = density(spks_time_mlist[[id_node]], bw = bw, from = min(t_vec), to = max(t_vec))
    yList[[id_node]] = res_smooth$y
    tList[[id_node]] = res_smooth$x
  }
  
  
  # Fit FPCA ------------------------------------
  time_start = Sys.time()
  FPCAobj = FPCA(Ly=yList, Lt=tList, optns=list(dataType='Dense', maxK=N_component))
  time_end = Sys.time()
  time_estimation = time_end - time_start
  time_estimation = as.numeric(time_estimation, units='secs')
  
  mean_density_vec = FPCAobj$mu
  eigenfuncs_mat = FPCAobj$phi[ , 1:N_component, drop=FALSE] # len(t_vec) x N_component
  fpc_scores_mat = FPCAobj$xiEst[ , 1:N_component, drop=FALSE] # N_node x N_component
  t_fpca_vec = FPCAobj$workGrid
  
  
  # Reconstruct estimates in our model -----
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
  center_density_array_est = array(dim = c(N_clus, N_component, length(t_vec)))
  v_mat_list_est = list() 
  for (id_component in 1:N_component){
    derivative_mean_density_vec = calculate_deriv(x_vec = t_fpca_vec, f_vec = mean_density_vec)
    eigen_func_vec = eigenfuncs_mat[,id_component]
    inner_prod = sum( derivative_mean_density_vec * eigen_func_vec ) * (t_fpca_vec[2]-t_fpca_vec[1])
    derivative_density_shifted_vec = eigen_func_vec * inner_prod
    v_fpca_vec = -fpc_scores_mat[, id_component] * (inner_prod^(-1))
    v_mean_fpca_vec = -min(v_fpca_vec)
    if (v_mean_fpca_vec > (max(t_vec)-min(t_vec)) ) {
      v_mean_fpca_vec = 1/2 * data_param$timeshift_max_vec[id_component]
    }
    
    density_shifted_vec = calculate_antideriv(x_vec = t_fpca_vec, gradient_vec = derivative_density_shifted_vec)
    density_shifted_vec = approx(x = t_fpca_vec, y = density_shifted_vec, xout = t_vec)$y
    timeshift_density = -v_mean_fpca_vec
    if (timeshift_density < 0) {
      N_timegrid = round(abs(timeshift_density)/t_unit)
      density_vec = c(tail(density_shifted_vec, length(density_shifted_vec) - N_timegrid),
                      rep(0, N_timegrid))
    } else if (timeshift_density > 0) {
      N_timegrid = round(abs(timeshift_density)/t_unit)
      density_vec = c(rep(0, N_timegrid),
                      head(density_shifted_vec, length(density_shifted_vec) - N_timegrid))
    }
    center_density_array_est[1, id_component, ] = density_vec
    
    v_fpca_vec = v_fpca_vec + v_mean_fpca_vec
    v_mat_list_est[[id_component]] = matrix(data = v_fpca_vec, nrow = N_node, ncol = N_replicate)
  }
  center_density_array_est = round(center_density_array_est, digits = 4)
  
  
  # Find the best permutation of density components -----
  permn_list = combinat::permn(1:N_component)  
  mise_f_min = Inf
  the_permn = c()
  id_clus = 1
  for (permn in permn_list) {
    mise_f_tmp_vec = c()
    for (id_component in 1:N_component) {
      t_unit = data_generated$t_vec[2] - data_generated$t_vec[1]
      f_target = data_generated$center_density_array_true[id_clus, id_component, ]
      density_est = center_density_array_est[id_clus, permn[id_component], ]
      res_ccf = ccf(y = density_est, x = f_target, plot = FALSE)
      n0_init = res_ccf$lag[which.max(res_ccf$acf)]
      f_origin_mat = matrix(density_est, nrow = 1)
      n0 = align_multi_components(f_target = f_target,
                                  f_origin_mat = f_origin_mat,
                                  t_unit = t_unit, 
                                  n0_vec = c(n0_init),
                                  n0_min_vec = -length(f_target) %/% 2,
                                  n0_max_vec = length(f_target) %/% 2 )$n0_vec
      if (n0 > 0) {
        density_est_shift = c(rep(0, n0), head(density_est, length(density_est) - n0) )
      } else if (n0 < 0) {
        density_est_shift = c(tail(density_est, length(density_est) - abs(n0)), rep(0, abs(n0)) )
      } else if (n0 == 0) {
        density_est_shift = density_est
      }
      mise_f_tmp_vec[id_component] = sum( t_unit * (density_est_shift - 
                                                      data_generated$center_density_array_true[id_clus, id_component, ])^2 )
    }
    mise_f_tmp = mean(mise_f_tmp_vec)
    if (mise_f_tmp < mise_f_min){
      the_permn = permn
      mise_f_min = mise_f_tmp
    }
  }
  center_density_array_est[id_clus, , ] = center_density_array_est[id_clus, the_permn, ]
  v_mat_list_est = v_mat_list_est[the_permn]
  
  
  # Other estimates
  N_clus_est = 1
  clusters_list_est = list(1:N_node)
  center_intensity_array_est = NA
  center_Nspks_mat_est = NA
  loss_history = NA
  N_iteration = NA
  
  
  # Compute estimation error ------------------------------------------------
  if (TRUE) {
    # Compute errors of conn patts, i.e. F ---------
    
    permn = 1
    
    center_density_array_est_permn = center_density_array_est[permn, , ,drop=FALSE]
    center_intensity_array_est_permn = NA
    center_Nspks_mat_est_permn = NA
    clusters_list_est_permn = clusters_list_est[permn]
    
    
    ### Calculate distance 
    dist_mse_mat = matrix(nrow=N_clus, ncol=N_component)
    for (id_clus in 1:N_clus) {
      for (id_component in 1:N_component) {
        f_target = center_density_array_true[id_clus,id_component,]
        density_est = center_density_array_est_permn[id_clus,id_component,]
        res_ccf = ccf(y = density_est, x = f_target, plot = FALSE)
        n0_init = res_ccf$lag[which.max(res_ccf$acf)]
        f_origin_mat = matrix(density_est, nrow = 1)
        n0 = align_multi_components(f_target = f_target,
                                    f_origin_mat = f_origin_mat,
                                    t_unit = t_unit, 
                                    n0_vec = c(n0_init),
                                    n0_min_vec = -length(f_target) %/% 2,
                                    n0_max_vec = length(f_target) %/% 2 )$n0_vec
        if (n0 > 0) {
          density_est_shift = c(rep(0, n0), head(density_est, length(density_est) - n0) )
        } else if (n0 < 0) {
          density_est_shift = c(tail(density_est, length(density_est) - abs(n0)), rep(0, abs(n0)) )
        } else if (n0 == 0) {
          density_est_shift = density_est
        }
        dist_mse_mat[id_clus,id_component] = sum( (density_est_shift - 
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
              time_estimation=time_estimation,
              N_iteration=NA,
              loss_history=NA
  ))
  
}


