### Estimate time shifts between each subject and each cluster

est_timeshift_predecomp = function(subjtrial_density_smooth_array,
                         fft_subjtrial_density_unsmooth_array,
                         N_spks_mat,
                         v_trialwise_vec_list,
                         center_density_array,
                         v_mat_list=NULL,
                         N_component=1,
                         freq_trun=5, 
                         t_vec=seq(0, 1, by=0.01),
                         key_times_vec = c(min(t_vec), 0, max(t_vec)),
                         fix_timeshift=FALSE,
                         rand_init = FALSE,
                         fix_comp1_timeshift_only=FALSE,
                         bw=0)
{
  t_unit = t_vec[2]-t_vec[1]
  N_subj = dim(subjtrial_density_smooth_array)[1]
  N_trial = dim(subjtrial_density_smooth_array)[2]
  N_clus = dim(center_density_array)[1]
  
  # Decompose the counting processes into components (i.e. S^v_{i,m}f_{z_i,m})
  subj_component_array = array(0, dim=c(N_subj, N_component, length(t_vec)))
  for (id_subj in 1:N_subj) {
    subjtrial_density_smooth_array_tmp = subjtrial_density_smooth_array[id_subj, , ,drop=FALSE]
    fft_subjtrial_density_unsmooth_array_tmp = fft_subjtrial_density_unsmooth_array[id_subj, , ,drop=FALSE]
    N_spks_mat_tmp = N_spks_mat[id_subj, ,drop=FALSE]
    v_mat_list_tmp = lapply(v_trialwise_vec_list, function(vec)matrix(vec, nrow=1))
    tmp = get_center_intensity_array(subjtrial_density_unsmooth_array = subjtrial_density_smooth_array_tmp, 
                                     fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array_tmp, 
                                     N_spks_mat = N_spks_mat_tmp, 
                                     v_trialwise_vec_list = v_trialwise_vec_list, 
                                     clusters_list = list(1), 
                                     v_mat_list = v_mat_list_tmp, 
                                     N_component = N_component, 
                                     freq_trun = freq_trun, 
                                     t_vec = t_vec, 
                                     key_times_vec = key_times_vec, 
                                     bw = bw, 
                                     fix_timeshift = fix_timeshift)
    subj_component_array[id_subj, , ] = tmp$center_density_array
  }
  
  # 
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
    
    n0_mat_current = matrix(0, nrow = N_subj, ncol = N_component)
    n0_mat_update = matrix(0, nrow = N_subj, ncol = N_component)
    for (id_component in 1:N_component) {
      f_target_array = subj_component_array[ ,id_component, ,drop=FALSE]
      id_trial = 1
      v_subj_comp = v_mat_list[[id_component]][1:N_subj, id_trial] - v_trialwise_vec_list[[id_component]][id_trial]
      n0_subj_comp = round(v_subj_comp / t_unit)
      n0_mat_current[1:N_subj, id_component] = n0_subj_comp
      if (fix_timeshift) {
        n0_mat_update[1:N_subj, id_component] = n0_mat_current[1:N_subj, id_component]
        v_array_list[[id_component]][1:N_subj, 1:N_trial, id_clus] = v_mat_list[[id_component]][1:N_subj, 1:N_trial]
      } else{
        ### TODO: Set n0_max_vec according to identifiability assumptions
        n0_max_vec = round( (key_times_vec[1+id_component] - key_times_vec[id_component]) / t_unit)
        n0_min_vec = 0
        if (rand_init) {
          n0_min_vec = -1 * n0_max_vec
        }
        n0_mat_update[, id_component] = align_multi_components(f_target_array = f_target_array,
                                                               f_origin_mat = f_origin_mat[id_component, ,drop=FALSE],
                                                               v_trialwise_vec_list = list(0),
                                                               N_spks_mat = matrix(rowMeans(N_spks_mat), ncol=1),
                                                               n0_init_mat = n0_mat_current[1:N_subj, id_component,drop=FALSE],
                                                               n0_min_vec = n0_min_vec,
                                                               n0_max_vec = n0_max_vec, 
                                                               # pad = 0,
                                                               t_unit = t_unit)$n0_mat
        
        v_subjwise_vec = n0_mat_update[1:N_subj, id_component]*t_unit
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_array_list[[id_component]][1:N_subj,1:N_trial,id_clus] = outer(v_subjwise_vec, v_trialwise_vec, "+")
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


