### Input: spks_time_mlist: N_node * N_replicate with each element being a list of spike times
### Perform algorithm based on intensities 
do_cluster_pdf = function(spks_time_mlist, stim_onset_vec, 
                          reaction_time_vec=NULL,
                          # Initial values
                          clusters_list_init,
                          v_vec_init=NULL,
                          v_mat_list_init=NULL,
                          # Tuning parameter
                          N_clus=length(clusters_list_init), 
                          N_component=1,
                          freq_trun=5, 
                          bw = 0,
                          v0 = 0.15, v1 = 0.1,
                          t_vec=seq(0, v0, length.out=200),
                          t_vec_extend=t_vec,
                          MaxIter=10, conv_thres=5e-3, 
                          fix_timeshift=FALSE,
                          fix_comp1_timeshift_only=FALSE,
                          fix_membership=FALSE,
                          gamma=0.06,
                          # Unused arguments
                          n0_vec_list_init=NULL,
                          opt_radius=max(t_vec)/2,
                          ...)
{
  if(is.matrix(stim_onset_vec)){
    stim_onset_vec = stim_onset_vec[,1]
  }
  
  
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  
  v_vec_init->v_vec
  v_mat_list_init->v_mat_list
  if(N_component==1){
    clusters_history = list()
    loss_history = c()
    
    # Update cluster-specific intensities and clusters ---------------------
    ### Get init values
    clusters_list = clusters_list_init
    v_vec = v_vec_init
    
    ### Save estimation
    clusters_history = c(clusters_history, list(clusters_list))
    
    ### Estimate parameters 
    clusters_list_update = clusters_list_current = clusters_list
    v_vec_update = v_vec_current = v_vec
    res = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                     stim_onset_vec = stim_onset_vec, 
                                     reaction_time_vec = reaction_time_vec, 
                                     clusters_list = clusters_list, 
                                     v_vec = v_vec,
                                     N_component = N_component,
                                     freq_trun = freq_trun,
                                     v0 = v0, v1 = v1,
                                     t_vec = t_vec,
                                     rmv_conn_prob = TRUE)
    center_density_array = res$center_density_array
    center_Nspks_mat = res$center_Nspks_mat
    center_intensity_array = res$center_intensity_array
    # browser()
    center_density_array_update = center_density_array_current = center_density_array
    center_Nspks_mat_update = center_Nspks_mat_current = center_Nspks_mat
    
    n_iter = 1
    stopping = FALSE
    while (!stopping & n_iter<=MaxIter){
      ### Update cluster-specific intensities and clusters once
      res = cluster_kmeans_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec, 
                               reaction_time_vec = reaction_time_vec, 
                               clusters_list = clusters_list_current, 
                               v_vec = v_vec_current,
                               N_component = N_component, 
                               freq_trun = freq_trun, 
                               v0 = v0, v1 = v1,
                               t_vec = t_vec,
                               fix_timeshift=fix_timeshift,
                               gamma=gamma,
                               ...
      )
      clusters_list_update = res$clusters_list
      center_density_array_update = res$center_density_array
      center_Nspks_mat_update = res$center_Nspks_mat
      center_intensity_array = res$center_intensity_array
      v_vec_update = res$v_vec
      
      l2_loss = res$l2_loss
      loss_history = c(loss_history, l2_loss)
      
      ### Evaluate stopping criterion
      clusters_update = clusters_list_update
      clusters_current = clusters_list_current
      delta_clusters = 1 - get_one_ARI(memb_est_vec = clus2mem(clusters_update), 
                                       memb_true_vec = clus2mem(clusters_current))
      
      delta_center_density_vec = sapply(1:N_component, function(id_component){
        sqrt( sum((abs(center_density_array_update[,id_component,]-center_density_array_current[,id_component,]))^2) / 
                (sum(abs(center_density_array_current[,id_component,])^2) + .Machine$double.eps) )
      } )
      
      delta_center_Nspks_vec = sapply(1:N_component, function(id_component){
        sqrt( sum((abs(center_Nspks_mat_update[,id_component]-center_Nspks_mat_current[,id_component]))^2) / 
                (sum(abs(center_Nspks_mat_current[,id_component])^2) + .Machine$double.eps) )
      } )
      
      delta_v_vec = sqrt( sum((v_vec_update-v_vec_current)^2) / 
                            (sum((v_vec_current)^2) + .Machine$double.eps) )
      
      stopping = mean(c(delta_clusters,delta_center_density_vec,delta_center_Nspks_vec,delta_v_vec)) < conv_thres
      
      
      ### *update -> *current
      n_iter = n_iter+1
      clusters_list_update -> clusters_list_current
      center_density_array_update -> center_density_array_current 
      center_Nspks_mat_update -> center_Nspks_mat_current
      v_vec_update -> v_vec_current
      
      clusters_history = c(clusters_history, list(clusters_list_current))
    }
    
    if (n_iter>MaxIter) {
      message("[do_cluster_pdf]: Reached maximum iteration number.")
    }
    N_iteration = n_iter
    
    
    
    # Get final result --------------------------------------------------------
    
    clusters_list_current -> clusters_list
    center_density_array_current -> center_density_array
    center_Nspks_mat_current -> center_Nspks_mat
    center_intensity_array = center_intensity_array
    v_vec_current -> v_vec
    
  } else if(N_component==2){
    clusters_history = list()
    loss_history = c()
    
    # Update cluster-specific intensities and clusters ---------------------
    ### Get init values
    clusters_list = clusters_list_init
    v_mat_list = v_mat_list_init
    
    ### Save estimation
    clusters_history = c(clusters_history, list(clusters_list))
    
    ### Estimate parameters 
    clusters_list_update = clusters_list_current = clusters_list
    res = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                     stim_onset_vec = stim_onset_vec, 
                                     clusters_list = clusters_list, 
                                     v_mat_list = v_mat_list,
                                     N_component = N_component,
                                     freq_trun = Inf,
                                     bw = bw,
                                     v0 = v0, v1 = v1,
                                     t_vec = t_vec,
                                     rmv_conn_prob = TRUE)
    center_density_array = res$center_density_array
    center_Nspks_mat = res$center_Nspks_mat
    center_intensity_array = res$center_intensity_array
    v_mat_list = res$v_mat_list
    
    center_density_array_update = center_density_array_current = center_density_array
    center_Nspks_mat_update = center_Nspks_mat_current = center_Nspks_mat
    v_mat_list_update = v_mat_list_current = v_mat_list
    
    n_iter = 1
    stopping = FALSE
    while (!stopping & n_iter<=MaxIter){
      ### Update cluster-specific intensities and clusters once
      res = cluster_kmeans_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec, 
                               clusters_list = clusters_list_current, 
                               v_mat_list = v_mat_list_current,
                               N_component = N_component, 
                               freq_trun = freq_trun, 
                               bw = bw,
                               v0 = v0, v1 = v1,
                               t_vec = t_vec,
                               fix_timeshift=fix_timeshift,
                               fix_comp1_timeshift_only=fix_comp1_timeshift_only,
                               fix_membership=fix_membership,
                               gamma=gamma,
                               ...
      )
      clusters_list_update = res$clusters_list
      center_density_array_update = res$center_density_array
      center_Nspks_mat_update = res$center_Nspks_mat
      center_intensity_array = res$center_intensity_array
      v_mat_list_update = res$v_mat_list
      
      l2_loss = res$l2_loss
      loss_history = c(loss_history, l2_loss)
      
      ### Evaluate stopping criterion
      clusters_update = clusters_list_update
      clusters_current = clusters_list_current
      delta_clusters = 1 - get_one_ARI(memb_est_vec = clus2mem(clusters_update), 
                                       memb_true_vec = clus2mem(clusters_current))
      
      delta_center_density_vec = sapply(1:N_component, function(id_component){
        sqrt( sum((abs(center_density_array_update[,id_component,]-center_density_array_current[,id_component,]))^2) / 
                (sum(abs(center_density_array_current[,id_component,])^2) + .Machine$double.eps) )
      } )
      
      delta_center_Nspks_vec = sapply(1:N_component, function(id_component){
        sqrt( sum((abs(center_Nspks_mat_update[,id_component]-center_Nspks_mat_current[,id_component]))^2) / 
                (sum(abs(center_Nspks_mat_current[,id_component])^2) + .Machine$double.eps) )
      } )
      
      delta_v_vec = sapply(1:N_component, function(id_component){
        sqrt( sum((v_mat_list_update[[id_component]]-v_mat_list_current[[id_component]])^2) / 
                (sum((v_mat_list_current[[id_component]])^2) + .Machine$double.eps) )
      })
        
      stopping = mean(c(delta_clusters,delta_center_density_vec,delta_center_Nspks_vec,delta_v_vec)) < conv_thres
      
      
      ### *update -> *current
      n_iter = n_iter+1
      clusters_list_update -> clusters_list_current
      center_density_array_update -> center_density_array_current 
      center_Nspks_mat_update -> center_Nspks_mat_current
      v_mat_list_update -> v_mat_list_current
      
      clusters_history = c(clusters_history, list(clusters_list_current))
    }
    
    if (n_iter>MaxIter) {
      message("[do_cluster_pdf]: Reached maximum iteration number.")
    }
    N_iteration = n_iter
    
    
    
    # Get final result --------------------------------------------------------
    
    clusters_list_current -> clusters_list
    center_density_array_current -> center_density_array
    center_Nspks_mat_current -> center_Nspks_mat
    center_intensity_array = center_intensity_array
    v_mat_list_current -> v_mat_list
    
    
    ### Extend estimated densities and intensities to t_vec_extend
    if (length(t_vec)<length(t_vec_extend)){
      center_intensity_array_extend = array(dim=c(N_clus, N_component, length(t_vec_extend)))
      center_density_array_extend = array(dim=c(N_clus, N_component, length(t_vec_extend)))
      center_Nspks_mat_extend = matrix(nrow=N_clus, ncol=N_component)
      for (id_clus in 1:N_clus) {
        intensity_q_1 = center_intensity_array[id_clus, 1, ]
        u_0 = v1; u_1 = v0;
        intensity_q_1 = c(rep(0,length(t_vec_extend)-length(t_vec)),
                          intensity_q_1 )
        N_spks_q_1 = sum(intensity_q_1*t_unit)
        
        intensity_q_2 = center_intensity_array[id_clus, 2, ]
        intensity_q_2 = c(rep(0,length(t_vec_extend)-length(t_vec)),
                          intensity_q_2 )
        N_spks_q_2 = sum(intensity_q_2*t_unit)
        
        
        center_intensity_array_extend[id_clus, 1, ] = intensity_q_1
        center_intensity_array_extend[id_clus, 2, ] = intensity_q_2
        center_Nspks_mat_extend[id_clus, 1] = N_spks_q_1
        center_Nspks_mat_extend[id_clus, 2] = N_spks_q_2
        center_density_array_extend[id_clus, 1, ] = intensity_q_1/(N_spks_q_1+N_spks_q_2+.Machine$double.eps)
        center_density_array_extend[id_clus, 2, ] = intensity_q_2/(N_spks_q_1+N_spks_q_2+.Machine$double.eps)
      }
      center_intensity_array = center_intensity_array_extend
      center_Nspks_mat = center_Nspks_mat_extend
      center_density_array = center_density_array_extend
    }
    
    
  }
  
  
  return(list(clusters_list=clusters_list, 
              loss_history=loss_history,
              clusters_history=clusters_history, 
              center_density_array=center_density_array,
              center_Nspks_mat=center_Nspks_mat,
              center_intensity_array = center_intensity_array,
              v_vec=v_vec,
              v_mat_list=v_mat_list,
              t_vec=t_vec,
              t_vec_extend=t_vec_extend,
              N_iteration=N_iteration))
  
}


