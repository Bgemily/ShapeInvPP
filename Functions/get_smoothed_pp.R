### Smooth point process using frequency truncation or kernel smoothing
### If freq_trun == Inf, no frequncy truncation is applied
### If bw == 0, no kernel smoothing is applied
get_smoothed_pp = function(event_time_vec, 
                           t_vec, 
                           freq_trun=Inf, 
                           bw=0){
  if(length(event_time_vec)==0){
    fft_vec = 0 * t_vec
    intens_vec = 0 * t_vec
  } else if(length(event_time_vec)>=1) {
    tmp = get_adaptive_fft(event_time_vec = event_time_vec, 
                           freq_trun_max = freq_trun, 
                           t_vec = t_vec, 
                           bw=bw,
                           freq_trun_min = freq_trun)
    fft_tmp = tmp$fft_vec_best
    if (freq_trun < Inf){
      fft_vec = c(head(fft_tmp, freq_trun+1), 
                   rep(0, length(t_vec)-2*freq_trun-1),
                   tail(fft_tmp, freq_trun))
    } else if (freq_trun == Inf){
      fft_vec = fft_tmp
    }
    intens_vec = Re(fft(fft_vec, inverse = TRUE))
  }
  
  return(list(freq_trun = freq_trun, 
              fft_vec = fft_vec,
              intens_vec = intens_vec))
}


