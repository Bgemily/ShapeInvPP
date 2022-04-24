### Estimate time shifts between each node and each cluster

est_timeshift = function(spks_time_mlist, 
                         stim_onset_vec, 
                         center_density_array,
                         center_Nspks_mat,
                         v_vec,
                         N_component=1,
                         freq_trun=5, 
                         v0 = 0.15, v1 = 0.1,
                         t_vec=seq(0, v0, by=0.01),
                         step_size=0.0001,
                         fix_timeshift=FALSE,
                         bw_nodedsty=0.02)
{
  t_unit = t_vec[2]-t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  N_clus = dim(center_density_array)[1]
  
  
  v_mat = matrix(nrow = N_node, ncol = N_clus)
  if(fix_timeshift){
    v_mat = matrix(v_vec, nrow=N_node, ncol=N_clus)
  } else{
    tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist,
                                     stim_onset_vec = stim_onset_vec,
                                     clusters_list = mem2clus(1:N_node),
                                     v_vec = rep(0,N_node),
                                     N_component = N_component,
                                     freq_trun = freq_trun,
                                     v0 = v0, v1 = v1,
                                     t_vec = t_vec, 
                                     bw=bw_nodedsty)
    node_density_array = tmp$center_density_array

    for (q in 1:N_clus) {
      for (i in 1:N_node) {
        f_target_list = lapply(1:N_component, function(l) node_density_array[i,l, ])
        f_origin_list = lapply(1:N_component, function(l) center_density_array[q,l, ])
        
        n0_current = round(v_vec[i]/t_unit)
        n0_max = which(cumsum(f_target_list[[1]])*t_unit > 0.3)[1]
        n0_min = which(cumsum(f_target_list[[1]])*t_unit > 0.7)[1] - length(t_vec)
        n0_tmp = align_multi_curves_gd_v2(f_origin_list = f_origin_list, 
                                          f_target_list = f_target_list,
                                          step_size = step_size,
                                          t_unit = t_unit, 
                                          n0 = n0_current,
                                          n0_min = n0_min,
                                          n0_max = n0_max, 
                                          pad = 0,
                                          periodic = TRUE)$n0
        v_mat[i,q] = n0_tmp*t_unit
        # browser(expr = (q==2 ))
        # plot(f_origin_list[[1]],type='l'); lines(f_target_list[[1]],col=2)
        
      }
      
    }
    
  }
 
  
  
  return(list(v_mat=v_mat))
}


