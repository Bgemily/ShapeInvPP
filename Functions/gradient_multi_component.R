gradient_multi_component = function(fft_f_target_mat, 
                                    fft_f_origin_mat,
                                    n0_vec,
                                    v_trialwise_vec_list,
                                    N_spks_trialwise_vec,
                                    t_unit)
{
  N_replicate = nrow(fft_f_target_mat)
  N_timegrid = ncol(fft_f_target_mat)
  N_component = nrow(fft_f_origin_mat)
  if(N_timegrid != ncol(fft_f_origin_mat)) 
    stop("Length of fft_f_target_mat and fft_f_origin's do not match.")
  
  for (id_replicate in 1:N_replicate) {
    fft_curr_comp = fft_f_target_mat[id_replicate, ]
    fft_f_target_mat[id_replicate, ] = c(tail(fft_curr_comp, (N_timegrid-1)%/%2), 
                                         head(fft_curr_comp, N_timegrid-(N_timegrid-1)%/%2) )
  }
  for (id_component in 1:N_component) {
    fft_curr_comp = fft_f_origin_mat[id_component, ]
    fft_f_origin_mat[id_component, ] = c(tail(fft_curr_comp, (N_timegrid-1)%/%2), 
                                         head(fft_curr_comp, N_timegrid-(N_timegrid-1)%/%2) )
  }
  
  l_vec = 0:(N_timegrid-1)
  l_vec = c(tail(l_vec, (N_timegrid-1)%/%2)-N_timegrid,
            head(l_vec, N_timegrid-(N_timegrid-1)%/%2) )
  
  gd_vec = c()
  for (id_component in 1:N_component) {
    gd_curr_comp_vec = c()
    for (id_replicate in 1:N_replicate) {
      fft_f_target = fft_f_target_mat[id_replicate, ]
      fft_f_origin_no_curr_comp = 0
      for (id_component_2 in 1:N_component) {
        if (id_component_2 != id_component) {
          n0_trialwise = round(v_trialwise_vec_list[[id_component_2]][id_replicate] / t_unit)
          fft_tmp = exp(1i*2*pi*(-(n0_vec[[id_component_2]]+n0_trialwise))*l_vec/N_timegrid) * fft_f_origin_mat[id_component_2, ]
          fft_f_origin_no_curr_comp = fft_f_origin_no_curr_comp + fft_tmp
        }
      }
      n0_trialwise = round(v_trialwise_vec_list[[id_component]][id_replicate] / t_unit)
      gd_curr_comp = 2 * sum( Re(1i*2*pi*(l_vec/N_timegrid) * 
                                   exp(1i*2*pi*(-(n0_vec[[id_component]]+n0_trialwise))*l_vec/N_timegrid)*fft_f_origin_mat[id_component, ] * 
                                   Conj(fft_f_target - fft_f_origin_no_curr_comp) ) ) 
      gd_curr_comp_vec[id_replicate] = gd_curr_comp
    }
    gd_vec[id_component] = sum(gd_curr_comp_vec * N_spks_trialwise_vec)
  }
  
  return(gd_vec)
}
