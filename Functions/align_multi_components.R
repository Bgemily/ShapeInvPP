align_multi_components = function(f_target_array,
                                  f_origin_mat,
                                  v_trialwise_vec_list,
                                  N_spks_mat,
                                  n0_init_mat,
                                  t_unit = 0.05, 
                                  n0_min_vec = 0,
                                  n0_max_vec = ncol(f_target_array), 
                                  pad = NULL,
                                  MaxIter=1000, 
                                  stopping_redu=0.01)
{
  
  N_subj = dim(f_target_array)[1]
  N_trial = dim(f_target_array)[2]
  N_component = nrow(f_origin_mat)
  n0_min_mat = matrix(n0_min_vec, byrow = TRUE, nrow = N_subj, ncol = N_component)
  n0_max_mat = matrix(n0_max_vec, byrow = TRUE, nrow = N_subj, ncol = N_component)
  
  ### Extend f_target_array and f_origin from [0,T] to [-T,2T] ----
  N_timetick = dim(f_target_array)[3]
  N_timetick_extend = (N_timetick-1) + N_timetick + 2*round(N_timetick/2)
  f_target_array_extend = array(dim = c(N_subj, N_trial, N_timetick_extend))
  if(!is.null(pad)){
    f_target_array_extend[ , , 1:(N_timetick-1)] = pad
    f_target_array_extend[ , , N_timetick:(2*N_timetick-1)] = f_target_array
    f_target_array_extend[ , , (2*N_timetick):(2*N_timetick+2*round(N_timetick/2)-1)] = pad
  } else{
    f_target_array_extend[ , , 1:(N_timetick-1)] = f_target_array[ , , 1]
    f_target_array_extend[ , , N_timetick:(2*N_timetick-1)] = f_target_array
    f_target_array_extend[ , , (2*N_timetick):(2*N_timetick+2*round(N_timetick/2)-1)] = f_target_array[ , , N_timetick]
  }
  f_target_array = f_target_array_extend
  
  f_origin_mat_extend = matrix(nrow = N_component, ncol = N_timetick_extend)
  if(!is.null(pad)){
    f_origin_mat_extend[ , 1:(N_timetick-1)] = pad
    f_origin_mat_extend[ , N_timetick:(2*N_timetick-1)] = f_origin_mat
    f_origin_mat_extend[ , (2*N_timetick):(2*N_timetick+2*round(N_timetick/2)-1)] = pad
  } else{
    f_origin_mat_extend[ , 1:(N_timetick-1)] = f_origin_mat[ , 1]
    f_origin_mat_extend[ , N_timetick:(2*N_timetick-1)] = f_origin_mat
    f_origin_mat_extend[ , (2*N_timetick):(2*N_timetick+2*round(N_timetick/2)-1)] = f_origin_mat[ , N_timetick]
  }
  f_origin_mat = f_origin_mat_extend

  ### Compute terms needed in gradients ----
  fft_f_target_array = 0 * f_target_array
  for (id_subj in 1:N_subj) {
    for (id_trial in 1:N_trial) {
      fft_f_target_array[id_subj, id_trial, ] = fft(f_target_array[id_subj, id_trial, ]) #/ length(f_target_array[id_trial, ])
    }
  }
  
  fft_f_origin_mat = 0 * f_origin_mat
  for (id_component in 1:N_component) {
    fft_f_origin_mat[id_component, ] = fft(f_origin_mat[id_component, ]) #/ length(f_origin_mat[id_component, ])
  }
  
  
  ### Gradient descent ----
  n0_mat = n0_init_mat
  iter_count = 0
  dist_redu_vec = rep(Inf, N_subj)
  dist_curr_vec = rep(Inf, N_subj)
  converge_vec = rep(FALSE, N_subj)
  while (!all(converge_vec) && iter_count<MaxIter) {
    iter_count = iter_count + 1
    ##### Update time shifts -----
    gd_mat = gradient_multi_component(fft_f_target_array = fft_f_target_array[!converge_vec, , , drop = FALSE], 
                                      fft_f_origin_mat = fft_f_origin_mat,
                                      n0_mat = n0_mat[!converge_vec, , drop = FALSE],
                                      v_trialwise_vec_list = v_trialwise_vec_list,
                                      N_spks_mat = N_spks_mat[!converge_vec, , drop = FALSE],
                                      t_unit = t_unit )
    hessian_array = hessian_multi_component(fft_f_target_array = fft_f_target_array[!converge_vec, , , drop = FALSE], 
                                          fft_f_origin_mat = fft_f_origin_mat,
                                          n0_mat = n0_mat[!converge_vec, , drop = FALSE],
                                          v_trialwise_vec_list = v_trialwise_vec_list,
                                          N_spks_mat = N_spks_mat[!converge_vec, , drop = FALSE],
                                          t_unit = t_unit )
    for (i in 1:sum(!converge_vec)) {
      id_subj = which(!converge_vec)[i]
      inv_hessian_mat_tmp = try(solve(hessian_array[i, , ]), silent = TRUE)
      if (identical(class(inv_hessian_mat_tmp), "try-error")) {
        delta_n0 = pmax(-round(N_timetick/10), pmin(round(N_timetick/10), (diag(as.matrix(hessian_array[i, , ])) + .Machine$double.eps)^(-1) * gd_mat[i, ] ) )
        n0_mat[id_subj, ] = n0_mat[id_subj, ] - delta_n0
      } else {
        delta_n0 = pmax(-round(N_timetick/10), pmin(round(N_timetick/10), inv_hessian_mat_tmp %*% gd_mat[i, ] ))
        n0_mat[id_subj, ] = n0_mat[id_subj, ] - delta_n0
      }
    }
    n0_mat[which(n0_mat<n0_min_mat)] = n0_min_mat[which(n0_mat<n0_min_mat)]
    n0_mat[which(n0_mat>n0_max_mat)] = n0_max_mat[which(n0_mat>n0_max_mat)]

    ### Calculate shifted center densities ----
    l_vec = 0:(N_timetick_extend-1)
    l_vec = c( head(l_vec, N_timetick_extend-(N_timetick_extend-1)%/%2),
               tail(l_vec, (N_timetick_extend-1)%/%2) - N_timetick_extend )
    
    fft_f_origin_shifted_array = array(data = 0, dim = c(sum(!converge_vec), N_trial, N_timetick_extend))
    for (id_component in 1:N_component) {
      n0_subjwise_tmp_vec = n0_mat[!converge_vec, id_component]
      n0_trialwise_tmp_vec = round(v_trialwise_vec_list[[id_component]] / t_unit)
      n0_subj_trial_tmp_mat = outer(n0_subjwise_tmp_vec, n0_trialwise_tmp_vec, FUN = "+")

      fft_f_origin_tmp = fft_f_origin_mat[id_component, ]
      fft_f_origin_tmp_array = outer(matrix(data = 1, nrow = sum(!converge_vec), ncol = N_trial), fft_f_origin_tmp)

      fft_curr_comp_shifted_array = exp(1i*2*pi*outer(-n0_subj_trial_tmp_mat, l_vec)/N_timetick_extend) * fft_f_origin_tmp_array
      fft_f_origin_shifted_array = fft_f_origin_shifted_array + fft_curr_comp_shifted_array
    }
    
    
    ### Calculate distance between observed point processes and shifted center densities ----
    diff_fft_squared_array = abs( fft_f_origin_shifted_array - fft_f_target_array[!converge_vec, , , drop = FALSE] )^2
    d_mat = apply(diff_fft_squared_array, MARGIN = c(1,2), FUN = sum)
    d_mat = sqrt(t_unit / N_timetick_extend * d_mat)
    dist_upd_vec = rowSums(d_mat * N_spks_mat[!converge_vec, , drop = FALSE])
    
    ### Evaluate convergence criterion ----
    dist_redu_vec = rep(0, N_subj)
    dist_redu_vec[!converge_vec] = (dist_curr_vec - dist_upd_vec) / (dist_upd_vec+.Machine$double.eps)
    dist_redu_vec[which(is.na(dist_redu_vec))] = 0
    converge_vec = I(dist_redu_vec < stopping_redu)
    
    dist_curr_vec = dist_upd_vec[!converge_vec]
  }

  if (iter_count == MaxIter) {
    warning("Reached max iteration number when estimating a time shift. Consider adjusting the step size.")
  }
  
  return(list(n0_mat = n0_mat))
  
  
}
