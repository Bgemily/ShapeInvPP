get_indiv_intensity_array = function(spks_time_mlist, 
                                     freq_trun,
                                     t_vec,
                                     bw)
{  
  t_unit = t_vec[2]-t_vec[1]
  N_subj = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  
  subj_density_array = array(dim = c(N_subj, N_trial, length(t_vec)))
  subj_intensity_array = array(dim = c(N_subj, N_trial, length(t_vec)))
  fft_subj_density_array = array(dim = c(N_subj, N_trial, length(t_vec)))
  N_spks_mat = matrix(0, nrow = N_subj, ncol = N_trial)
  for (id_subj in 1:N_subj) {
    for (id_trial in 1:N_trial) {
      spks_time_subjtrial = unlist(spks_time_mlist[id_subj,id_trial]) 
      spks_time_vec = spks_time_subjtrial[which(spks_time_subjtrial >= min(t_vec) & spks_time_subjtrial <= max(t_vec))]
      tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                            freq_trun = freq_trun, 
                            t_vec = t_vec, 
                            bw = bw)
      subj_intensity = tmp$intens_vec
      subj_intensity_array[id_subj, id_trial, ] = subj_intensity
      
      subj_density = subj_intensity / (length(spks_time_vec)+.Machine$double.eps)
      subj_density_array[id_subj, id_trial, ] = subj_density
      
      fft_subj_density_tmp = fft(subj_density) / length(subj_density)
      fft_subj_density_array[id_subj, id_trial, ] = fft_subj_density_tmp
      
      N_spks_mat[id_subj, id_trial] = length(spks_time_vec)
    }
  }
  
  
  return(list(subj_density_array = subj_density_array,
              subj_intensity_array = subj_intensity_array,
              fft_subj_density_array = fft_subj_density_array,
              N_spks_mat = N_spks_mat ))
}


