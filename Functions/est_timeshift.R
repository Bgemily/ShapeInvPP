### Estimate time shifts between each node and each cluster

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
  
  if(N_component==1){
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
                                       bw=bw)
      node_density_array = tmp$center_density_array
      
      for (id_clus in 1:N_clus) {
        for (id_node in 1:N_node) {
          f_target_list = lapply(1:N_component, function(l) node_density_array[id_node,l, ])
          f_origin_list = lapply(1:N_component, function(l) center_density_array[id_clus,l, ])
          
          n0_current = round(v_vec[id_node]/t_unit)
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
          v_mat[id_node,id_clus] = n0_tmp*t_unit
          
        }
        
      }
    }
    
    
  } else if(N_component == 2) {
    v_array_list = list(array(dim = c(N_node, N_replicate, N_clus)),
                        array(dim = c(N_node, N_replicate, N_clus)) )
    dist_mat = matrix(nrow = N_node, ncol = N_clus)
    for (id_clus in 1:N_clus) {
      ### Get smoothed center densities
      center_density_1 = center_density_array[id_clus,1, ]   
      center_density_2 = center_density_array[id_clus,2, ]
      if(freq_trun<Inf){
        fft_tmp = fft(center_density_1) / length(t_vec)
        fft_tmp_trun = c(head(fft_tmp, freq_trun+1), 
                         rep(0, length(t_vec)-2*freq_trun-1),
                         tail(fft_tmp, freq_trun))
        center_density_1 = Re(fft(fft_tmp_trun, inverse = TRUE))
        
        fft_tmp = fft(center_density_2) / length(t_vec)
        fft_tmp_trun = c(head(fft_tmp, freq_trun+1), 
                         rep(0, length(t_vec)-2*freq_trun-1),
                         tail(fft_tmp, freq_trun))
        center_density_2 = Re(fft(fft_tmp_trun, inverse = TRUE))
      }
      
      for (id_node in 1:N_node) {
        dist_tmp_vec = c()
        for (id_replicate in 1:N_replicate) {
          spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_replicate]) - stim_onset_vec[id_replicate]
          spks_time_vec = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) & 
                                                            spks_time_nodetrial<=max(t_vec))]
          
          ### Smooth point process 
          tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                                freq_trun = freq_trun, 
                                t_vec = t_vec, 
                                bw=bw)
          node_intensity = tmp$intens_vec
          node_density = node_intensity/length(spks_time_vec)
          
          
          ### Update time shifts of two components
          f_target = node_density
          f_origin_1 = center_density_1
          f_origin_2 = center_density_2
          
          n0_vec_current = sapply(v_mat_list, function(v_mat){round(v_mat[id_node, id_replicate]/t_unit)})
          if(fix_timeshift){
            n0_tmp_vec = n0_vec_current
            v_array_list[[1]][id_node,id_replicate,id_clus] = v_mat_list[[1]][id_node,id_replicate]
            v_array_list[[2]][id_node,id_replicate,id_clus] = v_mat_list[[2]][id_node,id_replicate]
            
          } else{
            # u_0 = v1, u_1 = v0
            n0_max_vec = c(round((v1/2)/t_unit),
                           round((v0-v1/2)/t_unit) )
            n0_min_vec = c(-round((v1/2)/t_unit),
                           -round((v0-v1/2)/t_unit) )
            n0_tmp_vec = align_two_components(f_target = f_target,
                                              f_origin_1 = f_origin_1,
                                              f_origin_2 = f_origin_2,
                                              step_size = step_size,
                                              t_unit = t_unit, 
                                              n0_vec = n0_vec_current,
                                              n0_min_vec = n0_min_vec,
                                              n0_max_vec = n0_max_vec, 
                                              # pad = 0,
                                              periodic = TRUE)$n0_vec
            
            
            if( fix_comp1_timeshift_only ){
              v_array_list[[1]][id_node,id_replicate,id_clus] = v_mat_list[[1]][id_node,id_replicate]
              v_array_list[[2]][id_node,id_replicate,id_clus] = n0_tmp_vec[2]*t_unit
            } else{
              v_array_list[[1]][id_node,id_replicate,id_clus] = n0_tmp_vec[1]*t_unit
              v_array_list[[2]][id_node,id_replicate,id_clus] = n0_tmp_vec[2]*t_unit
            }
          }
          
          ### Get un-smoothed node density and center densities
          tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                                freq_trun = Inf, 
                                t_vec = t_vec, 
                                bw=0)
          node_intensity = tmp$intens_vec
          node_density = node_intensity/length(spks_time_vec)
          center_density_1 = center_density_array[id_clus,1, ]   
          center_density_2 = center_density_array[id_clus,2, ]

          ### Calculate distance between (id_node, id_replicate) and id_clus
          fft_node_density = fft(node_density)/length(node_density)
          fft_center_density_1 = fft(center_density_1)/length(center_density_1)
          fft_center_density_2 = fft(center_density_2)/length(center_density_2)
          N = length(node_density)
          l_vec = 0:(N-1)
          l_vec = c( head(l_vec, N-(N-1)%/%2),
                     tail(l_vec, (N-1)%/%2) - N )
          dist_tmp_vec[id_replicate] = sum(abs( exp(-1i*2*pi*l_vec*n0_tmp_vec[1]/N)*fft_center_density_1 + 
                                                  exp(-1i*2*pi*l_vec*n0_tmp_vec[2]/N)*fft_center_density_2 -
                                                  fft_node_density )^2) * length(spks_time_vec)
          
        }
        
        
        dist_mat[id_node,id_clus] = sum(dist_tmp_vec)
      }
      
    }
    
    
    
  }
  
  
  
  return(list(v_mat=v_mat,
              v_array_list=v_array_list,
              dist_mat=dist_mat))
}


