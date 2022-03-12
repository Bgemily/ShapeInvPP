

### Obtain truncated fourier series (smoothed point process) for each cluster and intensity component
get_center_intensity_array = function(spks_time_mlist, stim_onset_vec, reaction_time_vec,
                                      clusters_list, 
                                      N_component=2,
                                      freq_trun=5,
                                      v0 = 0.2, v1 = 0.1,
                                      t_vec=seq(0, max(reaction_time_vec-stim_onset_vec+v0)+0.01, by=0.01),
                                      n0_mat_list=NULL, 
                                      rmv_conn_prob=FALSE)
{  
  time_unit = t_vec[2]-t_vec[1]
  N_clus = length(clusters_list)
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  
  center_intensity_array = array(dim=c(N_clus, N_component, length(t_vec)))

  for (q in 1:N_clus) {
    intensity_vis = intensity_act = 0
    if(length(clusters_list[[q]])>0){
      for (id_trial in 1:N_trial) {
        spks_time_rltvto_stim_onset = unlist(spks_time_mlist[clusters_list[[q]], id_trial]) - stim_onset_vec[id_trial]
        critical_time = reaction_time_vec[id_trial]-v1-stim_onset_vec[id_trial]
        id_spks_time_vis = which(spks_time_rltvto_stim_onset >= 0 & 
                                   spks_time_rltvto_stim_onset <=
                                   critical_time)
        spks_time_rltvto_stim_onset_vis = spks_time_rltvto_stim_onset[id_spks_time_vis]
        if (length(spks_time_rltvto_stim_onset_vis)>0) {
          fft_vis_res = get_adaptive_fft(event_time_vec = spks_time_rltvto_stim_onset_vis, 
                                         freq_trun_max = freq_trun, 
                                         t_vec = t_vec)
          fft_vis = fft_vis_res$fft_vec_best
        }
        else{
          fft_vis = rep(0,2*freq_trun+1)
        }
        intensity_vis_tmp = Re(fft(c(tail(fft_vis, freq_trun+1), 
                                     rep(0, length(t_vec)-2*freq_trun-1),
                                     head(fft_vis, freq_trun)), inverse = TRUE))
        intensity_vis_tmp[which(t_vec>critical_time)] = 0
        intensity_vis = intensity_vis + intensity_vis_tmp
        
        
        spks_time_rltvto_reaction_time = unlist(spks_time_mlist[clusters_list[[q]], id_trial]) - reaction_time_vec[id_trial]
        t_vec_2 = t_vec+max(stim_onset_vec)-max(reaction_time_vec)
        critical_time_2 = reaction_time_vec[id_trial]+v0-stim_onset_vec[id_trial]
        id_spks_time_act = which(spks_time_rltvto_stim_onset > critical_time &
                                   spks_time_rltvto_stim_onset <= critical_time_2)
        spks_time_rltvto_reaction_time = spks_time_rltvto_reaction_time[id_spks_time_act]
        if (length(spks_time_rltvto_reaction_time)>0) {
          fft_act_res = get_adaptive_fft(event_time_vec = spks_time_rltvto_reaction_time, 
                                         freq_trun_max = freq_trun, 
                                         t_vec = t_vec_2)
          fft_act = fft_act_res$fft_vec_best
        } else{
          fft_act = rep(0, 2*freq_trun+1)
        }
        intensity_act_tmp = Re(fft(c(tail(fft_act, freq_trun+1), 
                                     rep(0, length(t_vec)-2*freq_trun-1),
                                     head(fft_act, freq_trun)), inverse = TRUE))
        intensity_act_tmp[which(t_vec_2<=-v1)] = 0
        intensity_act = intensity_act + intensity_act_tmp
      }
      intensity_vis = intensity_vis / (N_trial*length(clusters_list[[q]]))
      intensity_act = intensity_act / (N_trial*length(clusters_list[[q]]))
      
    }
    
    center_intensity_array[q,1,] = intensity_vis
    center_intensity_array[q,2,] = intensity_act

  }
  
  return(center_intensity_array)
}


