
### Obtain truncated fourier series (smoothed point process) for each cluster and intensity component
get_center_intensity_array = function(spks_time_mlist, stim_onset_vec, reaction_time_vec,
                                      clusters_list, 
                                      v_vec,
                                      N_component=1,
                                      freq_trun=5,
                                      v0 = 0.15, v1 = 0.1,
                                      t_vec=seq(0, v0, by=0.01),
                                      n0_mat_list=NULL, 
                                      rmv_conn_prob=FALSE)
{  
  t_unit = t_vec[2]-t_vec[1]
  N_clus = length(clusters_list)
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  
  center_intensity_array = array(dim=c(N_clus, 1, length(t_vec)))
  N_spks_vec = rep(0,N_clus)
  for (q in 1:N_clus) {
    intensity_q = 0
    N_spks_q = 0
    if(length(clusters_list[[q]])>0){
      spks_time_q = c()
      for (id_node in clusters_list[[q]]) {
        for (id_trial in 1:N_trial) {
          spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_trial]) - stim_onset_vec[id_trial] - v_vec[id_node]
          spks_time_nodetrial = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) & 
                                                            spks_time_nodetrial<=max(t_vec))]
          spks_time_q = c(spks_time_q, spks_time_nodetrial)
          
        }
      }
      if (length(spks_time_q)>0) {
        fft_q_res = get_adaptive_fft(event_time_vec = spks_time_q, 
                                       freq_trun_max = freq_trun, 
                                       t_vec = t_vec)
        fft_q = fft_q_res$fft_vec_best
      } else if(freq_trun<Inf) {
        fft_q = rep(0,2*freq_trun+1)
      } else{
        fft_q = rep(0,length(t_vec))
      }
      intensity_q = fft2f(fft = fft_q, freq_trun = freq_trun, t_vec = t_vec)
      intensity_q = intensity_q / (N_trial*length(clusters_list[[q]]))

      N_spks_q = length(spks_time_q)
    }
    
    if(rmv_conn_prob){
      if(N_spks_q>0){
        intensity_q = intensity_q*(N_trial*length(clusters_list[[q]]))/N_spks_q
      }
    }
    
    center_intensity_array[q,1,] = intensity_q
    N_spks_vec[q] = N_spks_q
    
  }
  
  return(list(center_intensity_array=center_intensity_array,
              N_spks_vec=N_spks_vec))
}


