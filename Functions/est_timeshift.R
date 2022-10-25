### Estimate time shifts between each subject and each cluster

est_timeshift = function(spks_time_mlist, 
                         stim_onset_vec, 
                         v_trialwise_vec_list = NULL,
                         center_density_array,
                         v_vec=NULL,
                         v_mat_list=NULL,
                         N_component=1,
                         freq_trun=5, 
                         v0 = 0.15, v1 = 0.1,
                         t_vec=seq(0, v0, by=0.01),
                         step_size = 1e-4,
                         fix_timeshift=FALSE,
                         rand_init = FALSE,
                         fix_comp1_timeshift_only=FALSE,
                         bw=0)
{
  t_unit = t_vec[2]-t_vec[1]
  N_subj = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  N_clus = dim(center_density_array)[1]
  
  
  v_mat = matrix(nrow = N_subj, ncol = N_clus)
  v_array_list = rep(list(array(dim = c(N_subj, N_replicate, N_clus))), N_component)
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
    
    for (id_subj in 1:N_subj) {
      # Smooth point process -------
      f_target_mat = matrix(nrow = N_replicate, ncol = length(t_vec))
      N_spks_trialwise_vec = rep(0, N_replicate)
      for (id_replicate in 1:N_replicate) {
        spks_time_nodetrial = unlist(spks_time_mlist[id_subj,id_replicate]) - stim_onset_vec[id_replicate]
        spks_time_vec = spks_time_nodetrial[which(spks_time_nodetrial >= min(t_vec) & spks_time_nodetrial <= max(t_vec))]
        tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                              freq_trun = freq_trun, 
                              t_vec = t_vec, 
                              bw = bw)
        subj_intensity_smooth = tmp$intens_vec
        subj_density_smooth = subj_intensity_smooth / length(spks_time_vec)
        f_target_mat[id_replicate, ] = subj_density_smooth
        N_spks_trialwise_vec[id_replicate] = length(spks_time_vec)
      }
      
      # Align smoothed point process and smoothed cluster-wise density components ----
      n0_vec_current = rep(0, N_component)
      for (id_component in 1:N_component) {
        id_replicate = 1
        v_subj_comp = v_mat_list[[id_component]][id_subj, id_replicate] - v_trialwise_vec_list[[id_component]][id_replicate]
        n0_subj_comp = round(v_subj_comp / t_unit)
        n0_vec_current[id_component] = n0_subj_comp
      }
      if (fix_timeshift) {
        n0_vec_update = n0_vec_current
        for (id_component in 1:N_component) {
          v_array_list[[id_component]][id_subj, 1:N_replicate, id_clus] = v_mat_list[[id_component]][id_subj, 1:N_replicate]
        }
      } else{
        # u_0 = v1, u_1 = v0
        ### TODO: Set n0_max_vec to a reasonable value
        if (N_component == 2) {
          n0_max_vec = c(round((v1/2)/t_unit), round((v0-v1/2)/t_unit) )
        } else {
          n0_max_vec = rep(round((v1/1)/t_unit), N_component)
        }
        n0_min_vec = rep(0, N_component)
        if (rand_init) {
          n0_min_vec = -1 * n0_max_vec
        }
        n0_vec_update = align_multi_components(f_target_mat = f_target_mat,
                                               f_origin_mat = f_origin_mat,
                                               v_trialwise_vec_list = v_trialwise_vec_list,
                                               N_spks_trialwise_vec = N_spks_trialwise_vec,
                                               step_size = step_size,
                                               t_unit = t_unit, 
                                               n0_vec = n0_vec_current,
                                               n0_min_vec = n0_min_vec,
                                               n0_max_vec = n0_max_vec, 
                                               # pad = 0,
                                               periodic = TRUE)$n0_vec
        for (id_component in 1:N_component) {
          v_array_list[[id_component]][id_subj,1:N_replicate,id_clus] = n0_vec_update[id_component]*t_unit + v_trialwise_vec_list[[id_component]]
        }
        if (fix_comp1_timeshift_only) {
          v_array_list[[1]][id_subj,1:N_replicate,id_clus] = v_mat_list[[1]][id_subj,1:N_replicate]
        }
      }
      
      # Get un-smoothed center densities ----
      center_density_unsmooth_array = center_density_array
      fft_center_density_mat = matrix(nrow = N_component, ncol = dim(center_density_unsmooth_array)[3])
      for (id_component in 1:N_component) {
        center_density_vec_tmp = center_density_unsmooth_array[id_clus, id_component, ]   
        fft_center_density_mat[id_component, ] = fft(center_density_vec_tmp) / length(center_density_vec_tmp)
      }
      
      dist_tmp_vec = rep(0, N_replicate)
      for (id_replicate in 1:N_replicate) {
        # Get un-smoothed node density ----
        spks_time_nodetrial = unlist(spks_time_mlist[id_subj,id_replicate]) - stim_onset_vec[id_replicate]
        spks_time_vec = spks_time_nodetrial[which(spks_time_nodetrial >= min(t_vec) & spks_time_nodetrial <= max(t_vec))]
        tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                              freq_trun = Inf, 
                              t_vec = t_vec, 
                              bw = 0)
        subj_intensity_unsmooth = tmp$intens_vec
        subj_density_unsmooth = subj_intensity_unsmooth / length(spks_time_vec)
        
        # Calculate distance between (id_subj, id_replicate) and id_clus ----
        fft_subj_density = fft(subj_density_unsmooth) / length(subj_density_unsmooth)
        N = length(subj_density_unsmooth)
        l_vec = 0:(N-1)
        l_vec = c( head(l_vec, N-(N-1)%/%2),
                   tail(l_vec, (N-1)%/%2) - N )
        fft_center_density_shifted = 0
        for (id_component in 1:N_component) {
          n0_trialwise = round(v_trialwise_vec_list[[id_component]][id_replicate] / t_unit)
          fft_curr_comp_shifted = exp(-1i*2*pi*l_vec*(n0_vec_update[id_component]+n0_trialwise)/N) * fft_center_density_mat[id_component, ]
          fft_center_density_shifted = fft_center_density_shifted + fft_curr_comp_shifted
        }
        dist_tmp_vec[id_replicate] = sum(abs( fft_center_density_shifted - fft_subj_density )^2) * length(spks_time_vec)
      }
      dist_mat[id_subj,id_clus] = sum(dist_tmp_vec)
    }
    
  }
  

  return(list(v_mat = v_mat,
              v_array_list = v_array_list,
              dist_mat = dist_mat))
}


