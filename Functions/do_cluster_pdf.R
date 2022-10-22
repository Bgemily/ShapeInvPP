### Input: spks_time_mlist: N_node * N_replicate with each element being a list of spike times
do_cluster_pdf = function(spks_time_mlist, stim_onset_vec, 
                          reaction_time_vec=NULL,
                          # Initial values
                          center_density_array_init,
                          center_Nspks_mat_init,
                          clusters_list_init,
                          v_vec_init=NULL,
                          v_mat_list_init=NULL,
                          # Tuning parameter
                          N_clus,
                          N_component=1,
                          freq_trun=5, 
                          bw = 0,
                          v0 = 0.15, v1 = 0.1,
                          t_vec=seq(0, v0, length.out=200),
                          t_vec_extend=t_vec,
                          key_times_vec = c(min(t_vec), 0, max(t_vec)),
                          MaxIter=10, conv_thres=5e-3, 
                          fix_timeshift=FALSE,
                          rand_init = FALSE,
                          fix_comp1_timeshift_only=FALSE,
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
  
  clusters_history = list()
  loss_history = c()
  
  ### Save init estimation
  clusters_history = c(clusters_history, list(clusters_list_init))
  
  ### Estimate parameters 
  clusters_list_update = clusters_list_current = clusters_list_init
  center_density_array_update = center_density_array_current = center_density_array_init
  center_Nspks_mat_update = center_Nspks_mat_current = center_Nspks_mat_init
  v_mat_list_update = v_mat_list_current = v_mat_list_init
  l2_loss_update = l2_loss_current = Inf
  n_iter = 1
  stopping = FALSE
  while (!stopping & n_iter<=MaxIter){
    ### Update time shifts and clusters 
    tmp = get_timeshift_and_clusters(spks_time_mlist = spks_time_mlist,
                                     stim_onset_vec = stim_onset_vec,
                                     center_density_array = center_density_array_current,
                                     center_Nspks_mat = center_Nspks_mat_current,
                                     v_mat_list = v_mat_list_current,
                                     freq_trun = freq_trun,
                                     bw = bw,
                                     v0 = v0, v1 = v1,
                                     t_vec = t_vec,
                                     fix_timeshift = fix_timeshift,
                                     rand_init = rand_init,
                                     fix_comp1_timeshift_only = fix_comp1_timeshift_only,
                                     gamma = gamma)
    clusters_list_update = tmp$clusters_list
    v_mat_list_update = tmp$v_mat_list
    l2_loss = tmp$l2_loss
    loss_history = c(loss_history, l2_loss)
    
    ### Update intensities 
    tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                     stim_onset_vec = stim_onset_vec, 
                                     reaction_time_vec = reaction_time_vec, 
                                     clusters_list = clusters_list_update, 
                                     v_mat_list = v_mat_list_update,
                                     N_component = N_component,
                                     freq_trun = Inf, 
                                     bw = bw,
                                     v0 = v0, v1 = v1,
                                     t_vec = t_vec,
                                     key_times_vec = key_times_vec,
                                     fix_timeshift = fix_timeshift,
                                     align_density = FALSE)
    center_density_array_update = tmp$center_density_array
    center_Nspks_mat_update = tmp$center_Nspks_mat
    center_intensity_array = tmp$center_intensity_array
    v_mat_list_update = tmp$v_mat_list
    
    
    ### Evaluate stopping criterion
    l2_loss_update = l2_loss
    delta_loss = (l2_loss_current - l2_loss_update) / (l2_loss_update + .Machine$double.eps)
    stopping = delta_loss < conv_thres
    
    
    ### *update -> *current
    n_iter = n_iter+1
    clusters_list_update -> clusters_list_current
    center_density_array_update -> center_density_array_current 
    center_Nspks_mat_update -> center_Nspks_mat_current
    v_mat_list_update -> v_mat_list_current
    l2_loss_update -> l2_loss_current
    
    clusters_history = c(clusters_history, list(clusters_list_current))
  }
  
  if (n_iter>MaxIter) {
    message("[do_cluster_pdf]: Reached maximum iteration number.")
  }
  N_iteration = n_iter
  
  
  
  # Get final result --------------------------------------------------------
  
  clusters_list = clusters_list_current
  center_density_array = center_density_array_current
  center_Nspks_mat = center_Nspks_mat_current
  center_intensity_array = center_intensity_array
  v_mat_list = v_mat_list_current
  
  
  ### Extend estimated densities and intensities to t_vec_extend
  if (length(t_vec)<length(t_vec_extend)){
    center_intensity_array_extend = array(dim=c(N_clus, N_component, length(t_vec_extend)))
    center_density_array_extend = array(dim=c(N_clus, N_component, length(t_vec_extend)))
    center_Nspks_mat_extend = matrix(nrow=N_clus, ncol=N_component)
    for (id_clus in 1:N_clus) {
      for (id_component in 1:N_component) {
        intensity_q_curr_comp = center_intensity_array[id_clus, id_component, ]
        intensity_q_curr_comp = c(rep(0, length(t_vec_extend) - length(t_vec)),
                                  intensity_q_curr_comp )
        N_spks_q_curr_comp = sum(intensity_q_curr_comp * t_unit)
        center_intensity_array_extend[id_clus, id_component, ] = intensity_q_curr_comp
        center_Nspks_mat_extend[id_clus, id_component] = N_spks_q_curr_comp
      }
      for (id_component in 1:N_component) {
        intensity_q_curr_comp = center_intensity_array_extend[id_clus, id_component, ]
        N_spks_q_all_comp = sum(center_Nspks_mat_extend[id_clus, ])
        center_density_array_extend[id_clus, id_component, ] = intensity_q_curr_comp / (N_spks_q_all_comp + .Machine$double.eps)
      }
    }
    center_intensity_array = center_intensity_array_extend
    center_Nspks_mat = center_Nspks_mat_extend
    center_density_array = center_density_array_extend
  }
  
  
  return(list(clusters_list=clusters_list, 
              loss_history=loss_history,
              clusters_history=clusters_history, 
              center_density_array=center_density_array,
              center_Nspks_mat=center_Nspks_mat,
              center_intensity_array = center_intensity_array,
              v_mat_list = v_mat_list,
              t_vec = t_vec,
              t_vec_extend = t_vec_extend,
              N_iteration = N_iteration))
  
}


