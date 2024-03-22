### Implement algorithm 
do_cluster_pdf = function(# Observables
                          spks_time_mlist, 
                          v_trialwise_vec_list = NULL,
                          key_times_vec = c(min(t_vec), 0, max(t_vec)),
                          N_component=1,
                          # Initial values
                          clusters_list_init,
                          v_mat_list_init=NULL,
                          # Tuning parameter
                          N_clus=length(clusters_list_init), 
                          freq_trun=5, 
                          bw = 0,
                          t_vec=seq(0, 1, length.out=200),
                          v_subjwise_max=NULL,
                          MaxIter=10, conv_thres=5e-3, 
                          fix_timeshift=FALSE,
                          rand_init = FALSE,
                          gamma=0.06,
                          alpha=0 )
{
  
  t_unit = t_vec[2] - t_vec[1]
  N_subj = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  
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
                                   freq_trun = freq_trun,
                                   bw = bw,
                                   t_vec = t_vec, 
                                   alpha = alpha)
  center_density_array = res$center_density_array
  center_Nspks_mat = res$center_Nspks_mat
  center_intensity_array = res$center_intensity_array
  center_intensity_baseline_vec = res$center_intensity_baseline_vec
  v_mat_list = res$v_mat_list
  
  clusters_list_update = clusters_list_current = clusters_list_init
  center_density_array_update = center_density_array_current = center_density_array
  center_Nspks_mat_update = center_Nspks_mat_current = center_Nspks_mat
  center_intensity_array_update = center_intensity_array_current = center_intensity_array
  center_intensity_baseline_vec_update = center_intensity_baseline_vec_current = center_intensity_baseline_vec
  v_mat_list_update = v_mat_list_current = v_mat_list
  l2_loss_update = l2_loss_current = Inf
  
  clusters_history = list()
  center_density_array_history = list()
  clusters_history = c(clusters_history, list(clusters_list_init))
  center_density_array_history = c(center_density_array_history, list(center_density_array))
  
  v_subjwise_vec_list_current = list()
  for (id_component in 1:N_component) {
    id_trial_tmp = 1
    v_subjwise_vec_list_current[[id_component]] = v_mat_list[[id_component]][ ,id_trial_tmp] - v_trialwise_vec_list[[id_component]][id_trial_tmp]
  }
  v_subjwise_vec_list_history = list()
  v_subjwise_vec_list_history = c(v_subjwise_vec_list_history, list(v_subjwise_vec_list_current))
  

  n_iter = 1
  stopping = FALSE
  loss_history = c()
  while (!stopping & n_iter<=MaxIter){
    ### *update -> *current
    n_iter = n_iter+1
    clusters_list_update -> clusters_list_current
    center_density_array_update -> center_density_array_current 
    center_Nspks_mat_update -> center_Nspks_mat_current
    center_intensity_array_update -> center_intensity_array_current 
    center_intensity_baseline_vec_update -> center_intensity_baseline_vec_current
    v_mat_list_update -> v_mat_list_current
    l2_loss_update -> l2_loss_current
    
    ### Update intensities 
    tmp = get_center_intensity_array(subjtrial_density_unsmooth_array = subjtrial_density_unsmooth_array,
                                     fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                                     N_spks_mat = N_spks_mat,
                                     v_trialwise_vec_list = v_trialwise_vec_list,
                                     clusters_list = clusters_list_current, 
                                     v_mat_list = v_mat_list_current,
                                     N_component = N_component,
                                     freq_trun = freq_trun, 
                                     bw = bw,
                                     t_vec = t_vec,
                                     key_times_vec = key_times_vec,
                                     fix_timeshift = fix_timeshift, 
                                     alpha = alpha )
    center_density_array_update = tmp$center_density_array
    center_Nspks_mat_update = tmp$center_Nspks_mat
    center_intensity_array_update = tmp$center_intensity_array
    v_mat_list_tmp = tmp$v_mat_list
    center_density_baseline_vec_update = tmp$center_density_baseline_vec
    center_intensity_baseline_vec_update = tmp$center_intensity_baseline_vec
    
    v_subjwise_vec_list_current = list()
    for (id_component in 1:N_component) {
      id_trial_tmp = 1
      v_subjwise_vec_list_current[[id_component]] = v_mat_list_tmp[[id_component]][ ,id_trial_tmp] - v_trialwise_vec_list[[id_component]][id_trial_tmp]
    }
    v_subjwise_vec_list_history = c(v_subjwise_vec_list_history, list(v_subjwise_vec_list_current))
    
    clusters_history = c(clusters_history, list(clusters_list_current))
    center_density_array_history = c(center_density_array_history, list(center_density_array_update))
    
    ### Add baseline density to the first density component
    center_density_array_update_add_baseline = center_density_array_update
    center_Nspks_mat_update_add_baseline = center_Nspks_mat_update
    for (id_clus in 1:N_clus) {
      id_component_tmp = 1
      baseline_density = center_density_baseline_vec_update[id_clus]
      center_density_array_update_add_baseline[id_clus, id_component_tmp, ] = baseline_density + center_density_array_update[id_clus, id_component_tmp, ]
      baseline_intensity = center_intensity_baseline_vec_update[id_clus]
      center_Nspks_mat_update_add_baseline[id_clus, id_component_tmp] = baseline_intensity*(max(t_vec)-min(t_vec)) + center_Nspks_mat_update[id_clus, id_component_tmp]
    }
    
    ### Update time shifts and clusters 
    tmp = get_timeshift_and_clusters(subjtrial_density_smooth_array = subjtrial_density_smooth_array,
                                     fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                                     N_spks_mat = N_spks_mat,
                                     v_trialwise_vec_list = v_trialwise_vec_list,
                                     center_density_array = center_density_array_update_add_baseline,
                                     center_Nspks_mat = center_Nspks_mat_update_add_baseline,
                                     v_mat_list = v_mat_list_tmp,
                                     freq_trun = freq_trun,
                                     bw = bw,
                                     t_vec = t_vec,
                                     v_subjwise_max = v_subjwise_max,
                                     key_times_vec = key_times_vec,
                                     fix_timeshift = fix_timeshift,
                                     rand_init = FALSE,
                                     gamma = gamma)
    clusters_list_update = tmp$clusters_list
    v_mat_list_update = tmp$v_mat_list
    l2_loss = tmp$l2_loss
    for (id_clus in 1:N_clus) {
      for (id_component in 1:N_component) {
        deviance_from_constant = mean((center_density_array_update[id_clus,id_component,])^2)
        regularization = alpha * deviance_from_constant
        l2_loss = l2_loss + regularization
      }
    }
    
    dist_to_centr_vec = tmp$dist_to_centr_vec
    loss_history = c(loss_history, l2_loss)
    
    ### Evaluate stopping criterion
    l2_loss_update = l2_loss
    delta_loss = (l2_loss_current - l2_loss_update) / (l2_loss_update + .Machine$double.eps)
    stopping = delta_loss < conv_thres
    ### Save the two l2_loss terms concerning distribution and N_spks
    if (!stopping & n_iter<=MaxIter) {
      l2_loss_part_1 = tmp$l2_loss_part_1
      l2_loss_part_2 = tmp$l2_loss_part_2
    }
    
  }
  
  if (n_iter>MaxIter) {
    message("[do_cluster_pdf]: Reached maximum iteration number.")
  }
  N_iteration = n_iter
  
  ### Update intensities 
  if (FALSE) {
    tmp = get_center_intensity_array(subjtrial_density_unsmooth_array = subjtrial_density_unsmooth_array,
                                     fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                                     N_spks_mat = N_spks_mat,
                                     v_trialwise_vec_list = v_trialwise_vec_list,
                                     clusters_list = clusters_list_current, 
                                     v_mat_list = v_mat_list_current,
                                     N_component = N_component,
                                     freq_trun = freq_trun, 
                                     bw = bw,
                                     t_vec = t_vec,
                                     key_times_vec = key_times_vec,
                                     fix_timeshift = fix_timeshift,
                                     alpha = alpha)
    center_density_array_current = tmp$center_density_array
    center_Nspks_mat_current = tmp$center_Nspks_mat
    center_intensity_array_current = tmp$center_intensity_array
    v_mat_list_current = tmp$v_mat_list
  }
  
  
  
  # Get final result --------------------------------------------------------
  
  clusters_list = clusters_list_current
  center_density_array = center_density_array_current
  center_Nspks_mat = center_Nspks_mat_current
  center_intensity_array = center_intensity_array_current
  center_intensity_baseline_vec = center_intensity_baseline_vec_current
  v_mat_list = v_mat_list_current
  v_subjwise_vec_list = list()
  for (id_component in 1:N_component) {
    id_trial_tmp = 1
    v_subjwise_vec_list[[id_component]] = v_mat_list[[id_component]][ ,id_trial_tmp] - v_trialwise_vec_list[[id_component]][id_trial_tmp]
  }
  
  
  
  return(list(clusters_list = clusters_list, 
              loss_history = loss_history,
              l2_loss_part_1 = l2_loss_part_1,
              l2_loss_part_2 = l2_loss_part_2,
              clusters_history = clusters_history, 
              center_density_array_history = center_density_array_history,
              center_density_array = center_density_array,
              center_Nspks_mat = center_Nspks_mat,
              center_intensity_array = center_intensity_array,
              center_intensity_baseline_vec = center_intensity_baseline_vec,
              v_mat_list = v_mat_list,
              v_subjwise_vec_list = v_subjwise_vec_list,
              v_subjwise_vec_list_history = v_subjwise_vec_list_history,
              t_vec = t_vec,
              N_iteration = N_iteration))
  
}


