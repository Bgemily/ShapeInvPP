### Choose adaptive frequency truncation parameter to smooth a vector of event times
get_adaptive_fft = function(event_time_vec, 
                            freq_trun_max, 
                            t_vec, 
                            bw=0.03,
                            freq_trun_min=freq_trun_max)
{
  t_unit = t_vec[2]-t_vec[1]
  
  loss_min = Inf
  ### Get empirical intensity of event times
  ## V1
  # breaks = c(t_vec[1]-t_unit,t_vec)+t_unit/2
  # emp_intens_vec = hist(event_time_vec, breaks=breaks, plot=FALSE)$counts
  # emp_intens_vec = emp_intens_vec/t_unit
  ## V2
  if(length(event_time_vec)>1){
    density = density(event_time_vec,
                      bw=bw,
                      from=min(t_vec), to=max(t_vec),
                      n=length(t_vec))$y
    emp_intens_vec = density*length(event_time_vec)
  } else{
    event_time_vec = rep(event_time_vec[1],2)
    density = density(event_time_vec,
                      bw=bw,
                      from=min(t_vec), to=max(t_vec),
                      n=length(t_vec))$y
    emp_intens_vec = density
  }
  
  
  ### Get normalized fourier series
  emp_intens_fft = fft(emp_intens_vec) / length(t_vec)
  ### Get truncated fourier series
  if (freq_trun_max<Inf) {
    fft_vec_max = c(head(emp_intens_fft, freq_trun_max+1),
                    tail(emp_intens_fft, freq_trun_max))
  } else{
    fft_vec_max = emp_intens_fft
  }
  fft_vec_best = fft_vec_max
  freq_trun_best = freq_trun_max
  
  return(list(freq_trun_best = freq_trun_best, 
         fft_vec_best = fft_vec_best))
}