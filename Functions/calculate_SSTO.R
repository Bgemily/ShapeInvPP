calculate_SSTO = function(spks_time_mlist, 
                          t_vec)
{
  N_subj = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  t_unit = t_vec[2]-t_vec[1]
  
  # Calculate empirical distribution of spike times, and empirical N_spks
  intensity_empirical_mat = matrix(data = NA, nrow = N_subj*N_trial, ncol = length(t_vec))
  distribution_empirical_mat = matrix(data = NA, nrow = N_subj*N_trial, ncol = length(t_vec))
  N_spks_empirical_vec = rep(NA, N_subj*N_trial)
  id_row = 0
  for (id_subj in 1:N_subj){
    for (id_trial in 1:N_trial) {
      # empirical distribution of spike times for current (subj, trial)
      event_time_vec_tmp = unlist(spks_time_mlist[id_subj, id_trial])
      if(length(event_time_vec_tmp)==0){
        next
      }
      breaks = c(t_vec[1]-t_unit, t_vec)
      counts_tmp = hist(event_time_vec_tmp, breaks=breaks, plot=FALSE, right=FALSE)$counts
      intensity_empirical = counts_tmp / t_unit
      distribution_empirical = intensity_empirical/sum(intensity_empirical*t_unit)
      
      # empirical N_spks
      N_spks_empirical = length(event_time_vec_tmp)
      
      # Save empirical distribution and empirical N_spks
      id_row = id_row + 1
      intensity_empirical_mat[id_row, ] = intensity_empirical
      distribution_empirical_mat[id_row, ] = distribution_empirical
      N_spks_empirical_vec[id_row] = N_spks_empirical
      
    }
  }
  intensity_empirical_mat = intensity_empirical_mat[1:id_row, , drop=FALSE]
  distribution_empirical_mat = distribution_empirical_mat[1:id_row, , drop=FALSE]
  N_spks_empirical_vec = N_spks_empirical_vec[1:id_row]
  
  # Calculate SSTO for empirical distribution and empirical N_spks
  disrtibution_empirical_avg_vec = colMeans(distribution_empirical_mat)
  disrtibution_empirical_avg_mat = matrix(disrtibution_empirical_avg_vec, byrow = TRUE, nrow = nrow(distribution_empirical_mat), ncol = ncol(distribution_empirical_mat))
  SSTO_distribution = sum(N_spks_empirical_vec*max(t_vec)^(-1)*rowSums((distribution_empirical_mat - disrtibution_empirical_avg_mat)^2*t_unit))
  
  SSTO_N_spks = sum((mean(N_spks_empirical_vec))^(-1) * (N_spks_empirical_vec - mean(N_spks_empirical_vec))^2)
  
  # Output ------------------------------------------------------------------
  
  return(list(SSTO_distribution = SSTO_distribution,
              SSTO_N_spks = SSTO_N_spks))
}
