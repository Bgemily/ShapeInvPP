### Centering step: Update time shifts and connecting patterns alternatively

est_timeshift_density = function(spks_time_mlist, 
                                 stim_onset_vec, 
                                 reaction_time_vec, 
                                 clusters_list, 
                                 v_vec,
                                 N_component=1,
                                 freq_trun=5, 
                                 v0 = 0.15, v1 = 0.1,
                                 t_vec=seq(0, v0, by=0.01),
                                 step_size=0.0001,
                                 max_iter_centering=5, 
                                 epsilon=0.001,
                                 fix_timeshift=FALSE,
                                 bw_nodedsty=0.02)
{
  t_unit = t_vec[2]-t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  N_clus = length(clusters_list)
  
  
  # Update time shifts and connecting patterns alternatively ----------------
  n_iter = 1
  converge = FALSE
  v_vec_update = v_vec_current = v_vec
  tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                   stim_onset_vec = stim_onset_vec, 
                                   reaction_time_vec = reaction_time_vec, 
                                   clusters_list = clusters_list, 
                                   v_vec = v_vec_current,
                                   N_component = N_component,
                                   freq_trun = freq_trun, 
                                   t_vec = t_vec,
                                   v0 = v0, v1 = v1)
  center_density_array = tmp$center_density_array
  center_Nspks_mat = tmp$center_Nspks_mat
  center_intensity_array = tmp$center_intensity_array
  
  center_density_array_update = center_density_array_current = center_density_array
  center_Nspks_mat_update = center_Nspks_mat_current = center_Nspks_mat
  while(!converge && n_iter<= max_iter_centering)
  {
    ## Update connecting patterns ----
    tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                                        stim_onset_vec = stim_onset_vec, 
                                                        reaction_time_vec = reaction_time_vec, 
                                                        clusters_list = clusters_list, 
                                                        v_vec = v_vec_current,
                                                        N_component = N_component,
                                                        freq_trun = freq_trun, 
                                                        t_vec = t_vec,
                                                        v0 = v0, v1 = v1,
                                                        rmv_conn_prob = TRUE)
    center_density_array_update = tmp$center_density_array
    center_Nspks_mat_update = tmp$center_Nspks_mat
    center_intensity_array = tmp$center_intensity_array
    
    # Update time shifts ----
    if (fix_timeshift) {
      v_vec_update = v_vec_current
    } else{
      node_density_array = get_center_intensity_array(spks_time_mlist = spks_time_mlist,
                                                        stim_onset_vec = stim_onset_vec,
                                                        reaction_time_vec = reaction_time_vec,
                                                        clusters_list = mem2clus(1:N_node),
                                                        v_vec = v_vec_current*0,
                                                        N_component = N_component,
                                                        freq_trun = freq_trun,
                                                        v0 = v0, v1 = v1,
                                                        t_vec = t_vec, 
                                                        bw=bw_nodedsty)$center_density_array
      for (q in 1:N_clus) {
        n0_vec_tmp = rep(0,N_node)
        shifted_spks_time_q = c()
        for (i in clusters_list[[q]]) {
          f_target_list = lapply(1:N_component, function(l) node_density_array[i,l, ])
          f_origin_list = lapply(1:N_component, function(l) center_density_array_update[q,l, ])
          
          n0_current = round(v_vec[i]/t_unit)
          n0_max = which(cumsum(f_target_list[[1]])*t_unit > 0.5)[1]
          n0_min = which(cumsum(f_target_list[[1]])*t_unit > 0.5)[1] - length(t_vec)
          n0_vec_tmp[i] = align_multi_curves_gd_v2(f_origin_list = f_origin_list, 
                                                   f_target_list = f_target_list,
                                                   step_size = step_size,
                                                   t_unit = t_unit, 
                                                   n0 = n0_current,
                                                   n0_min = n0_min,
                                                   n0_max = n0_max, 
                                                   pad = 0,
                                                   periodic = TRUE)$n0
          # browser(expr = (q==2 ))
          # plot(f_origin_list[[1]],type='l'); lines(f_target_list[[1]],col=2)
          for (id_trial in 1:N_trial) {
            spks_time_vec_tmp = unlist(spks_time_mlist[i, id_trial]) - stim_onset_vec[id_trial]
            spks_time_vec_tmp = spks_time_vec_tmp[which(spks_time_vec_tmp>=min(t_vec)&
                                                          spks_time_vec_tmp<=max(t_vec))]
            shifted_spks_time_vec_tmp = spks_time_vec_tmp - n0_vec_tmp[i]*t_unit
            shifted_spks_time_vec_tmp = shifted_spks_time_vec_tmp[which(shifted_spks_time_vec_tmp>=min(t_vec)&
                                                                          shifted_spks_time_vec_tmp<=max(t_vec))]
            shifted_spks_time_q = c(shifted_spks_time_q, shifted_spks_time_vec_tmp)
          }
        }
        n0_vec_q = n0_vec_tmp[clusters_list[[q]]]
        v_vec_q = n0_vec_q*t_unit
        ### Shift time shifts in cluster-q s.t. median(shifted_spks_time_q) = median(t_vec)
        v_vec_q = v_vec_q + (median(shifted_spks_time_q)-median(t_vec))
      }
      v_vec_update = n0_vec_tmp*t_unit
      # v_vec_update = v_vec_update-min(v_vec_update)
    }
    
    ### Evaluate stopping criterion
    delta_center_density_vec = sapply(1:N_component, function(id_component){
      sqrt( sum((abs(center_density_array_update[,id_component,]-center_density_array_current[,id_component,]))^2) / 
              (sum(abs(center_density_array_current[,id_component,])^2) + .Machine$double.eps) )
    } )
    
    delta_v_vec = sqrt( sum((v_vec_update-v_vec_current)^2) / 
                          (sum((v_vec_current)^2) + .Machine$double.eps) )
    
    converge = mean(c(delta_center_density_vec,delta_v_vec)) < epsilon
    
    ### *update -> *current
    n_iter = n_iter + 1
    v_vec_current = v_vec_update
    center_density_array_current = center_density_array_update
    center_Nspks_mat_current = center_Nspks_mat_update
  }
  
  ### Get time shifts and connecting patterns
  v_vec = v_vec_current
  center_density_array = center_density_array_current
  center_Nspks_mat = center_Nspks_mat_current
  center_intensity_array = center_intensity_array
  
  return(list(center_density_array=center_density_array, 
              center_Nspks_mat = center_Nspks_mat,
              center_intensity_array = center_intensity_array,
              v_vec=v_vec))
}


