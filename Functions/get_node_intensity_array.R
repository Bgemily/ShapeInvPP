# TODO: Add v_vec in input
### Obtain truncated fourier series (smoothed shifted point process) for each node and each trial 

get_node_intensity_array = function(spks_time_mlist, stim_onset_vec, reaction_time_vec,
                                    clusters_list, 
                                    v_vec,
                                    N_component=2,
                                    freq_trun=5,
                                    v0 = 0.2, v1 = 0.1,
                                    t_vec=seq(0, max(reaction_time_vec-stim_onset_vec+v0)+0.01, by=0.01),
                                    rmv_conn_prob=FALSE)
{  
  time_unit = t_vec[2]-t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  node_intensity_array = array(dim=c(N_node, N_trial, length(t_vec)))
  for (id_node in 1:N_node) {
    for (id_trial in 1:N_trial) {
      t_vec_tmp = t_vec + stim_onset_vec[id_trial]
      spks_time = unlist(spks_time_mlist[id_node, id_trial])
      spks_time = spks_time[which(spks_time<=max(t_vec_tmp) & spks_time>=min(t_vec_tmp))]
      if (length(spks_time)>0) {
        fft_res = get_adaptive_fft(event_time_vec = spks_time, 
                                   freq_trun_max = freq_trun, 
                                   t_vec = t_vec_tmp)
        fft_tmp = fft_res$fft_vec_best
      }
      else{
        fft_tmp = rep(0,2*freq_trun+1)
      }
      # TODO: Check whether tail() and head() are flipped
      intensity_tmp = Re(fft(c(tail(fft_tmp, freq_trun+1), 
                                   rep(0, length(t_vec)-2*freq_trun-1),
                                   head(fft_tmp, freq_trun)), inverse = TRUE))
      node_intensity_array[id_node, id_trial, ] = intensity_tmp
    }
  }
  
  return(node_intensity_array)
}

