### Input: spks_time_mlist: N_node * N_trial with each element being a list of spike times
### Perform algorithm based on cumulative intensities 
do_cluster_pdf = function(spks_time_mlist, stim_onset_vec, reaction_time_vec,
                          # Initial values
                          clusters_list_init,
                          # Tuning parameter
                          N_component=2,
                          freq_trun=5, 
                          v0 = 0.2, v1 = 0.1,
                          t_vec=seq(0, max(reaction_time_vec-stim_onset_vec+v0)+0.01, by=0.01),
                          MaxIter=10, conv_thres=5e-3, 
                          # Unused arguments
                          N_clus=NULL, 
                          n0_vec_list_init=NULL,
                          step_size=200,
                          opt_radius=max(t_vec)/2,
                          fix_timeshift=FALSE,
                          gamma=0.001,
                          ...)
{
  if(is.matrix(stim_onset_vec)){
    stim_onset_vec = stim_onset_vec[,1]
  }
  
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  clusters_history = list()
  loss_history = c()
 
  # Update cluster-specific intensities and clusters ---------------------
  ### Get init clusters
  clusters_list = clusters_list_init
  
  ### Save estimation
  clusters_history = c(clusters_history, list(clusters_list))

  ### Estimate parameters 
  clusters_list_update = clusters_list_current = clusters_list
  center_intensity_array = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                                      stim_onset_vec = stim_onset_vec, 
                                                      reaction_time_vec = reaction_time_vec, 
                                                      clusters_list = clusters_list, 
                                                      N_component = N_component,
                                                      freq_trun = freq_trun,
                                                      v0 = v0, v1 = v1,
                                                      rmv_conn_prob = TRUE)
  center_intensity_array_update = center_intensity_array_current = center_intensity_array
  
  n_iter = 1
  stopping = FALSE
  while (!stopping & n_iter<=MaxIter){
    ### Update cluster-specific intensities and clusters once
    res = cluster_kmeans_pdf(spks_time_mlist = spks_time_mlist,
                             stim_onset_vec = stim_onset_vec, 
                             reaction_time_vec = reaction_time_vec, 
                             clusters_list = clusters_list_current, 
                             N_component = N_component, 
                             freq_trun = freq_trun, 
                             v0 = v0, v1 = v1)
    clusters_list_update = res$clusters_list
    center_intensity_array_update = res$center_intensity_array
    l2_loss = res$l2_loss
    
    loss_history = c(loss_history, l2_loss)
    
    ### Evaluate stopping criterion
    clusters_update = clusters_list_update
    clusters_current = clusters_list_current
    delta_clusters = 1 - get_one_ARI(memb_est_vec = clus2mem(clusters_update), 
                                     memb_true_vec = clus2mem(clusters_current))
    
    
    N_component = dim(center_intensity_array_update)[2]
    delta_center_intensity_vec = sapply(1:N_component, function(id_component){
      sqrt( sum((abs(center_intensity_array_update[,id_component,]-center_intensity_array_current[,id_component,]))^2) / 
        (sum(abs(center_intensity_array_current[,id_component,])^2) + .Machine$double.eps) )
    }
      )
    
    stopping = mean(c(delta_center_intensity_vec,delta_clusters)) < conv_thres
    
    
    ### *update -> *current
    n_iter = n_iter+1
    clusters_list_update -> clusters_list_current
    center_intensity_array_update -> center_intensity_array_current 
    
    clusters_history = c(clusters_history, list(clusters_list_current))
  }
  
  if (n_iter>MaxIter) {
    message("[do_cluster_pdf]: Reached maximum iteration number.")
  }
  N_iteration = n_iter
  
  
  
  # Get final result --------------------------------------------------------
  
  clusters_list_current -> clusters_list
  center_intensity_array_current -> center_intensity_array
  

  return(list(clusters_list=clusters_list, 
              loss_history=loss_history,
              clusters_history=clusters_history, 
              center_intensity_array=center_intensity_array, 
              t_vec=t_vec,
              N_iteration=N_iteration))
  
}


