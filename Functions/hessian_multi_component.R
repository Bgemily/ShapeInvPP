hessian_multi_component = function(fft_f_target_array,
                                   fft_f_origin_mat,
                                   n0_mat,
                                   v_trialwise_vec_list,
                                   N_spks_mat,
                                   t_unit)
{
  N_subj = dim(fft_f_target_array)[1]
  N_trial = dim(fft_f_target_array)[2]
  N_timetick = dim(fft_f_target_array)[3]
  N_component = nrow(fft_f_origin_mat)
  if(N_timetick != ncol(fft_f_origin_mat)) 
    stop("Length of fft_f_target_array and fft_f_origin's do not match.")
  
  idx_timetick = 1:N_timetick
  idx_timetick_swap = c(tail(idx_timetick, (N_timetick-1)%/%2), 
                        head(idx_timetick, N_timetick-(N_timetick-1)%/%2) )
  fft_f_target_array = fft_f_target_array[ , , idx_timetick_swap, drop = FALSE]
  fft_f_origin_mat = fft_f_origin_mat[ , idx_timetick_swap, drop = FALSE]
  
  l_vec = 0:(N_timetick-1)
  l_vec = c(tail(l_vec, (N_timetick-1)%/%2) - N_timetick,
            head(l_vec, N_timetick-(N_timetick-1)%/%2) )
  l_array = outer(matrix(data = 1, nrow = N_subj, ncol = N_trial), l_vec)
  
  hessian_array = array(dim = c(N_subj, N_component, N_component))
  for (id_component in 1:N_component) {
    fft_f_origin_no_curr_comp_array = array(data = 0, dim = c(N_subj, N_trial, N_timetick))
    for (id_component_2 in setdiff(1:N_component, id_component)) {
      n0_subjwise_tmp_vec = n0_mat[, id_component_2]
      n0_trialwise_tmp_vec = round(v_trialwise_vec_list[[id_component_2]] / t_unit)
      n0_subj_trial_tmp_mat = outer(n0_subjwise_tmp_vec, n0_trialwise_tmp_vec, FUN = "+")
      fft_f_origin_tmp = fft_f_origin_mat[id_component_2, ]
      fft_f_origin_tmp_array = outer(matrix(data = 1, nrow = N_subj, ncol = N_trial), fft_f_origin_tmp)
      fft_shifted_tmp_array = exp(1i*2*pi*outer(-n0_subj_trial_tmp_mat, l_vec)/N_timetick) * fft_f_origin_tmp_array
      fft_f_origin_no_curr_comp_array = fft_f_origin_no_curr_comp_array + fft_shifted_tmp_array
    }
    n0_subjwise_tmp_vec = n0_mat[, id_component]
    n0_trialwise_tmp_vec = round(v_trialwise_vec_list[[id_component]] / t_unit)
    n0_subj_trial_tmp_mat = outer(n0_subjwise_tmp_vec, n0_trialwise_tmp_vec, FUN = "+")
    fft_f_origin_tmp = fft_f_origin_mat[id_component, ]
    fft_f_origin_tmp_array = outer(matrix(data = 1, nrow = N_subj, ncol = N_trial), fft_f_origin_tmp)
    
    gd_order2_mat = matrix(nrow = N_subj, ncol = N_trial)
    gd_order2_mat = 2 * apply( Re( (-1i*2*pi*(l_array/N_timetick))^2 * 
                                     exp(1i*2*pi*outer(-n0_subj_trial_tmp_mat, l_vec)/N_timetick) * fft_f_origin_tmp_array * 
                                     Conj(fft_f_origin_no_curr_comp_array - fft_f_target_array) ), 
                               MARGIN = c(1,2), FUN = sum ) 
    hessian_array[1:N_subj, id_component, id_component] = rowSums(gd_order2_mat * N_spks_mat)
    
    for (id_component_2 in setdiff(1:N_component, id_component)) {
      n0_subjwise_tmp_vec_2 = n0_mat[, id_component_2]
      n0_trialwise_tmp_vec_2 = round(v_trialwise_vec_list[[id_component_2]] / t_unit)
      n0_subj_trial_tmp_mat_2 = outer(n0_subjwise_tmp_vec_2, n0_trialwise_tmp_vec_2, FUN = "+")
      fft_f_origin_tmp_2 = fft_f_origin_mat[id_component_2, ]
      fft_f_origin_tmp_array_2 = outer(matrix(data = 1, nrow = N_subj, ncol = N_trial), fft_f_origin_tmp_2)
      
      gd_order2_mat = matrix(nrow = N_subj, ncol = N_trial)
      gd_order2_mat = 2 * apply( Re( (-1i*2*pi*(l_array/N_timetick) * 
                                        exp(1i*2*pi*outer(-n0_subj_trial_tmp_mat, l_vec)/N_timetick) * 
                                        fft_f_origin_tmp_array) * 
                                       Conj(-1i*2*pi*(l_array/N_timetick) * 
                                              exp(1i*2*pi*outer(-n0_subj_trial_tmp_mat_2, l_vec)/N_timetick) * 
                                              fft_f_origin_tmp_array_2) ), 
                                 MARGIN = c(1,2), FUN = sum ) 
      hessian_array[1:N_subj, id_component, id_component_2] = rowSums(gd_order2_mat * N_spks_mat)
    }
  }
  
  return(hessian_array)
}
