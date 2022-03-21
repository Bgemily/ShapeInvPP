fft2f = function(fft, freq_trun, t_vec) {
  if (freq_trun<Inf) {
    f = Re(fft(c(head(fft, freq_trun+1), 
                 rep(0, length(t_vec)-2*freq_trun-1),
                 tail(fft, freq_trun)), inverse = TRUE))
  } else{
    f = Re(fft(fft, inverse = TRUE))
  }
  
  return(f)
}