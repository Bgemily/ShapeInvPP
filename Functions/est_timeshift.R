### Estimate time shifts between each subject and each cluster

est_timeshift = function(spks_time_mlist, 
                         stim_onset_vec, 
                         center_density_array,
                         v_vec=NULL,
                         v_mat_list=NULL,
                         N_component=1,
                         freq_trun=5, 
                         v0 = 0.15, v1 = 0.1,
                         t_vec=seq(0, v0, by=0.01),
                         step_size = 1e-4,
                         fix_timeshift=FALSE,
                         fix_comp1_timeshift_only=FALSE,
                         bw=0)
{
  t_unit = t_vec[2]-t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  N_clus = dim(center_density_array)[1]
  
  
  v_mat = matrix(nrow = N_node, ncol = N_clus)
  v_array_list = list(v_mat)
  dist_mat = matrix(nrow = N_node, ncol = N_clus)
  
  
  v_array_list = rep(list(array(dim = c(N_node, N_replicate, N_clus))), N_component)
  
  dist_mat = matrix(nrow = N_node, ncol = N_clus)
  center_density_smooth_array = 0 * center_density_array
  for (id_clus in 1:N_clus) {
    # Get smoothed center densities ----------
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
    
    # Estimate time shifts, and calculate distance between each subject and each cluster ----
    for (id_node in 1:N_node) {
      dist_tmp_vec = c()
      for (id_replicate in 1:N_replicate) {
        ### Smooth point process ----
        spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_replicate]) - stim_onset_vec[id_replicate]
        spks_time_vec = spks_time_nodetrial[which(spks_time_nodetrial >= min(t_vec) & spks_time_nodetrial <= max(t_vec))]
        tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                              freq_trun = freq_trun, 
                              t_vec = t_vec, 
                              bw = bw)
        node_intensity_smooth = tmp$intens_vec
        node_density_smooth = node_intensity_smooth / length(spks_time_vec)
        f_target = node_density_smooth
        
        # Align smoothed point process and smoothed cluster-wise density components ----
        n0_vec_update = rep(0, N_component)
        n0_vec_current = sapply(v_mat_list, function(v_mat){round(v_mat[id_node, id_replicate]/t_unit)})
        if (fix_timeshift) {
          n0_vec_update = n0_vec_current
          for (id_component in 1:N_component) {
            v_array_list[[id_component]][id_node,id_replicate,id_clus] = v_mat_list[[id_component]][id_node,id_replicate]
          }
        } else{
          # u_0 = v1, u_1 = v0
          ### TODO: Set n0_max_vec to a reasonable value
          if (N_component == 2) {
            n0_max_vec = c(round((v1/2)/t_unit), round((v0-v1/2)/t_unit) )
          } else {
            n0_max_vec = rep(round((v1/2)/t_unit), N_component)
          }
          n0_min_vec = rep(0, N_component)
          n0_vec_update = align_multi_components(f_target = f_target,
                                                 f_origin_mat = f_origin_mat,
                                                 step_size = step_size,
                                                 t_unit = t_unit, 
                                                 n0_vec = n0_vec_current,
                                                 n0_min_vec = n0_min_vec,
                                                 n0_max_vec = n0_max_vec, 
                                                 # pad = 0,
                                                 periodic = TRUE)$n0_vec
          for (id_component in 1:N_component) {
            v_array_list[[id_component]][id_node,id_replicate,id_clus] = n0_vec_update[id_component]*t_unit
          }
          if (fix_comp1_timeshift_only) {
            v_array_list[[1]][id_node,id_replicate,id_clus] = v_mat_list[[1]][id_node,id_replicate]
          }
        }
        
        # Get un-smoothed node density and center densities ----
        tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                              freq_trun = Inf, 
                              t_vec = t_vec, 
                              bw = 0)
        node_intensity_unsmooth = tmp$intens_vec
        node_density_unsmooth = node_intensity_unsmooth / length(spks_time_vec)
        center_density_unsmooth_array = center_density_array
        
        # Calculate distance between (id_node, id_replicate) and id_clus ----
        ### TODO: Move calculation of fft_center_density_mat out of this for loop
        fft_node_density = fft(node_density_smooth) / length(node_density_smooth)
        fft_center_density_mat = matrix(nrow = N_component, ncol = dim(center_density_unsmooth_array)[3])
        for (id_component in 1:N_component) {
          center_density_vec_tmp = center_density_unsmooth_array[id_clus, id_component, ]   
          fft_center_density_mat[id_component, ] = fft(center_density_vec_tmp) / length(center_density_vec_tmp)
        }
        N = length(node_density_smooth)
        l_vec = 0:(N-1)
        l_vec = c( head(l_vec, N-(N-1)%/%2),
                   tail(l_vec, (N-1)%/%2) - N )
        fft_center_density_shifted = 0
        for (id_component in 1:N_component) {
          fft_curr_comp_shifted = exp(-1i*2*pi*l_vec*n0_vec_update[id_component]/N) * fft_center_density_mat[id_component, ]
          fft_center_density_shifted = fft_center_density_shifted + fft_curr_comp_shifted
        }
        dist_tmp_vec[id_replicate] = sum(abs( fft_center_density_shifted - fft_node_density )^2) * length(spks_time_vec)
      }
      
      
      dist_mat[id_node,id_clus] = sum(dist_tmp_vec)
    }
    
  }
  

  return(list(v_mat = v_mat,
              v_array_list = v_array_list,
              dist_mat = dist_mat))
}


