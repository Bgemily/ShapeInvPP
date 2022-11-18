### Input: spks_time_mlist: N_subj * N_trial with each element being a list of spike times
### Perform algorithm based on intensities 
do_cluster_pdf = function(spks_time_mlist, 
                          v_trialwise_vec_list = NULL,
                          # Initial values
                          clusters_list_init,
                          v_vec_init=NULL,
                          v_mat_list_init=NULL,
                          # Tuning parameter
                          N_clus=length(clusters_list_init), 
                          N_component=1,
                          freq_trun=5, 
                          bw = 0,
                          t_vec=seq(0, 1, length.out=200),
                          t_vec_extend=t_vec,
                          key_times_vec = c(min(t_vec), 0, max(t_vec)),
                          MaxIter=10, conv_thres=5e-3, 
                          fix_timeshift=FALSE,
                          rand_init = FALSE,
                          fix_comp1_timeshift_only=FALSE,
                          gamma=0.06,
                          eta = 0,
                          ...)
{
  
  t_unit = t_vec[2] - t_vec[1]
  N_subj = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  clusters_history = list()
  loss_history = c()
  
  ### Save init estimation
  clusters_history = c(clusters_history, list(clusters_list_init))
  
  ### Calculate subj-trial-wise densities
  res = get_indiv_intensity_array(spks_time_mlist = spks_time_mlist, 
                                  freq_trun = freq_trun,
                                  bw = bw,
                                  t_vec = t_vec )
  subjtrial_density_smooth_array = res$subjtrial_density_array
  fft_subjtrial_density_smooth_array = res$fft_subjtrial_density_array
  
  res = get_indiv_intensity_array(spks_time_mlist = spks_time_mlist, 
                                  freq_trun = Inf,
                                  bw = 0,
                                  t_vec = t_vec )
  subjtrial_density_unsmooth_array = res$subjtrial_density_array
  fft_subjtrial_density_unsmooth_array = res$fft_subjtrial_density_array
  N_spks_mat = res$N_spks_mat 
  
  
  ### Estimate parameters 
  res = get_center_intensity_array(subjtrial_density_unsmooth_array = subjtrial_density_unsmooth_array,
                                   fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                                   N_spks_mat = N_spks_mat,
                                   v_trialwise_vec_list = v_trialwise_vec_list,
                                   clusters_list = clusters_list_init, 
                                   v_mat_list = v_mat_list_init,
                                   N_component = N_component,
                                   key_times_vec = key_times_vec,
                                   fix_timeshift = fix_timeshift,
                                   freq_trun = Inf,
                                   bw = bw,
                                   eta = eta,
                                   t_vec = t_vec)
  center_density_array = res$center_density_array
  center_Nspks_mat = res$center_Nspks_mat
  center_intensity_array = res$center_intensity_array
  v_mat_list = res$v_mat_list
  
  clusters_list_update = clusters_list_current = clusters_list_init
  center_density_array_update = center_density_array_current = center_density_array
  center_Nspks_mat_update = center_Nspks_mat_current = center_Nspks_mat
  v_mat_list_update = v_mat_list_current = v_mat_list
  l2_loss_update = l2_loss_current = Inf
  
  n_iter = 1
  stopping = FALSE
  while (!stopping & n_iter<=MaxIter){
    ### Update intensities 
    tmp = get_center_intensity_array(subjtrial_density_unsmooth_array = subjtrial_density_unsmooth_array,
                                     fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                                     N_spks_mat = N_spks_mat,
                                     v_trialwise_vec_list = v_trialwise_vec_list,
                                     clusters_list = clusters_list_current, 
                                     v_mat_list = v_mat_list_current,
                                     N_component = N_component,
                                     freq_trun = Inf, 
                                     bw = bw,
                                     eta = eta,
                                     t_vec = t_vec,
                                     key_times_vec = key_times_vec,
                                     fix_timeshift = fix_timeshift )
    center_density_array_update = tmp$center_density_array
    center_Nspks_mat_update = tmp$center_Nspks_mat
    center_intensity_array = tmp$center_intensity_array
    v_mat_list_tmp = tmp$v_mat_list
    
    ### Update time shifts and clusters 
    tmp = get_timeshift_and_clusters(subjtrial_density_smooth_array = subjtrial_density_smooth_array,
                                     fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                                     N_spks_mat = N_spks_mat,
                                     v_trialwise_vec_list = v_trialwise_vec_list,
                                     center_density_array = center_density_array_update,
                                     center_Nspks_mat = center_Nspks_mat_update,
                                     v_mat_list = v_mat_list_tmp,
                                     freq_trun = freq_trun,
                                     bw = bw,
                                     t_vec = t_vec,
                                     key_times_vec = key_times_vec,
                                     fix_timeshift = fix_timeshift,
                                     fix_comp1_timeshift_only = fix_comp1_timeshift_only,
                                     rand_init = FALSE,
                                     gamma = gamma)
    clusters_list_update = tmp$clusters_list
    v_mat_list_update = tmp$v_mat_list
    l2_loss = tmp$l2_loss
    dist_to_centr_vec = tmp$dist_to_centr_vec
    loss_history = c(loss_history, l2_loss)
    
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
  
  
  return(list(clusters_list = clusters_list, 
              loss_history = loss_history,
              dist_to_centr_vec = dist_to_centr_vec,
              clusters_history = clusters_history, 
              center_density_array = center_density_array,
              center_Nspks_mat = center_Nspks_mat,
              center_intensity_array = center_intensity_array,
              v_mat_list = v_mat_list,
              t_vec = t_vec,
              t_vec_extend = t_vec_extend,
              N_iteration = N_iteration))
  
}


