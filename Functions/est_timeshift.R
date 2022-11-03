### Estimate time shifts between each subject and each cluster

est_timeshift = function(subjtrial_density_smooth_array,
                         fft_subjtrial_density_unsmooth_array,
                         N_spks_mat,
                         v_trialwise_vec_list = NULL,
                         center_density_array,
                         v_vec=NULL,
                         v_mat_list=NULL,
                         N_component=1,
                         freq_trun=5, 
                         v0 = 0.15, v1 = 0.1,
                         t_vec=seq(0, v0, by=0.01),
                         key_times_vec = c(min(t_vec), 0, max(t_vec)),
                         step_size = 1e-4,
                         fix_timeshift=FALSE,
                         rand_init = FALSE,
                         fix_comp1_timeshift_only=FALSE,
                         bw=0)
{
  t_unit = t_vec[2]-t_vec[1]
  N_subj = dim(subjtrial_density_smooth_array)[1]
  N_trial = dim(subjtrial_density_smooth_array)[2]
  N_clus = dim(center_density_array)[1]
  
  v_mat = matrix(nrow = N_subj, ncol = N_clus)
  v_array_list = rep(list(array(dim = c(N_subj, N_trial, N_clus))), N_component)
  dist_mat = matrix(nrow = N_subj, ncol = N_clus)
  center_density_smooth_array = 0 * center_density_array
  for (id_clus in 1:N_clus) {
    # Smooth center densities -----
    if (freq_trun<Inf) {
      for (id_component in 1:N_component) {
        fft_tmp = fft(center_density_array[id_clus, id_component, ]) / length(t_vec)
        fft_tmp_trun = c(head(fft_tmp, freq_trun+1),
                         rep(0, length(t_vec)-2*freq_trun-1),
                         tail(fft_tmp, freq_trun))
        center_density_smooth_array[id_clus, id_component, ] = Re(fft(fft_tmp_trun, inverse = TRUE))
      }
    }
    f_origin_mat = matrix(nrow = N_component, ncol = length(t_vec))
    f_origin_mat[1:N_component, ] = center_density_smooth_array[id_clus, 1:N_component, ]
    
    # Smooth observed point process -------
    f_target_array = subjtrial_density_smooth_array
    
    # Align smoothed point process and smoothed cluster-wise density components ----
    n0_mat_current = matrix(0, nrow = N_subj, ncol = N_component)
    for (id_component in 1:N_component) {
      id_trial = 1
      v_subj_comp = v_mat_list[[id_component]][1:N_subj, id_trial] - v_trialwise_vec_list[[id_component]][id_trial]
      n0_subj_comp = round(v_subj_comp / t_unit)
      n0_mat_current[1:N_subj, id_component] = n0_subj_comp
    }
    if (fix_timeshift) {
      n0_mat_update = n0_mat_current
      for (id_component in 1:N_component) {
        v_array_list[[id_component]][1:N_subj, 1:N_trial, id_clus] = v_mat_list[[id_component]][1:N_subj, 1:N_trial]
      }
    } else {
      ### TODO: Set n0_max_vec according to identifiability assumptions
      if (FALSE) {
        n0_max_vec = c(round((v1/2)/t_unit), round((v0-v1/2)/t_unit) )
      } else {
        n0_max_vec = round( (key_times_vec[2:(N_component+1)] - key_times_vec[1:N_component]) / t_unit)
      }
      n0_min_vec = rep(0, N_component)
      if (rand_init) {
        n0_min_vec = -1 * n0_max_vec
      }
      n0_mat_update = align_multi_components(f_target_array = f_target_array,
                                             f_origin_mat = f_origin_mat,
                                             v_trialwise_vec_list = v_trialwise_vec_list,
                                             N_spks_mat = N_spks_mat,
                                             n0_init_mat = n0_mat_current,
                                             step_size = step_size,
                                             t_unit = t_unit, 
                                             n0_min_vec = n0_min_vec,
                                             n0_max_vec = n0_max_vec, 
                                             # pad = 0,
                                             periodic = TRUE)$n0_mat
      for (id_component in 1:N_component) {
        v_subjwise_vec = n0_mat_update[1:N_subj, id_component]*t_unit
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_array_list[[id_component]][1:N_subj,1:N_trial,id_clus] = outer(v_subjwise_vec, v_trialwise_vec, "+")
      }
      if (fix_comp1_timeshift_only) {
        v_array_list[[1]][1:N_subj,1:N_trial,id_clus] = v_mat_list[[1]][1:N_subj,1:N_trial]
      }
    }
    
    # Get un-smoothed center densities ----
    center_density_unsmooth_array = center_density_array
    fft_center_density_mat = matrix(nrow = N_component, ncol = dim(center_density_unsmooth_array)[3])
    for (id_component in 1:N_component) {
      center_density_vec_tmp = center_density_unsmooth_array[id_clus, id_component, ]   
      fft_center_density_mat[id_component, ] = fft(center_density_vec_tmp) / length(center_density_vec_tmp)
    }
    
    # Get un-smoothed subj density ----
    fft_subj_density_array = fft_subjtrial_density_unsmooth_array
    
    # Calculate distance between (id_subj, id_trial) and id_clus ----
    N_timetick = length(t_vec)
    l_vec = 0:(N_timetick-1)
    l_vec = c( head(l_vec, N_timetick-(N_timetick-1)%/%2),
               tail(l_vec, (N_timetick-1)%/%2) - N_timetick )
    fft_center_density_shifted_array = array(data = 0, dim = c(N_subj, N_trial, N_timetick))
    for (id_component in 1:N_component) {
      n0_subjwise_tmp_vec = n0_mat_update[1:N_subj, id_component]
      n0_trialwise_tmp_vec = round(v_trialwise_vec_list[[id_component]] / t_unit)
      n0_subj_trial_tmp_mat = outer(n0_subjwise_tmp_vec, n0_trialwise_tmp_vec, FUN = "+")
      fft_center_density_tmp = fft_center_density_mat[id_component, ]
      fft_center_density_tmp_array = outer(matrix(data = 1, nrow = N_subj, ncol = N_trial), fft_center_density_tmp)
      
      fft_curr_comp_shifted_array = exp(-1i*2*pi*outer(n0_subj_trial_tmp_mat, l_vec)/N_timetick) * fft_center_density_tmp_array
      fft_center_density_shifted_array = fft_center_density_shifted_array + fft_curr_comp_shifted_array
    }
    dist_tmp_mat = matrix(nrow = N_subj, ncol = N_trial)
    dist_tmp_mat = N_spks_mat * apply(abs(fft_center_density_shifted_array - fft_subj_density_array)^2, MARGIN = c(1,2), FUN = sum)
      
    dist_mat[1:N_subj, id_clus] = rowSums(dist_tmp_mat)
  }
  

  return(list(v_mat = v_mat,
              v_array_list = v_array_list,
              dist_mat = dist_mat))
}


