gradient_multi_component = function(fft_f_target, 
                                    fft_f_origin_mat,
                                    n0_vec)
{
  N_timegrid = length(fft_f_target)
  N_component = nrow(fft_f_origin_mat)
  if(N_timegrid != ncol(fft_f_origin_mat)) 
    stop("Length of fft_f_target and fft_f_origin's do not match.")
  
  fft_f_target = c(tail(fft_f_target, (N_timegrid-1)%/%2), 
                   head(fft_f_target, N_timegrid-(N_timegrid-1)%/%2) )
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
    fft_f_origin_no_curr_comp = 0
    for (id_component_2 in 1:N_component) {
      if (id_component_2 != id_component) {
        fft_tmp = exp(1i*2*pi*(-n0_vec[[id_component_2]])*l_vec/N_timegrid) * fft_f_origin_mat[id_component_2, ]
        fft_f_origin_no_curr_comp = fft_f_origin_no_curr_comp + fft_tmp
      }
    }
    gd_curr_comp = 2 * sum( Re(1i*2*pi*(l_vec/N_timegrid) * 
                                 exp(1i*2*pi*(-n0_vec[[id_component]])*l_vec/N_timegrid)*fft_f_origin_mat[id_component, ] * 
                                 Conj(fft_f_target - fft_f_origin_no_curr_comp) ) ) 
    gd_vec[id_component] = gd_curr_comp
  }
  
  return(gd_vec)
}
