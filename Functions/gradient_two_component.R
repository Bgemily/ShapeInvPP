gradient_two_component = function(fft_f_target, 
                                  fft_f_origin_1, 
                                  fft_f_origin_2,
                                  n0_vec)
  {
    N = length(fft_f_target)
    if(N != length(fft_f_origin_1) | N != length(fft_f_origin_2)) 
      stop("Length of fft_f_target and fft_f_origin's do not match.")
    
    fft_f_target = c(tail(fft_f_target, (N-1)%/%2), 
                     head(fft_f_target, N-(N-1)%/%2) )
    fft_f_origin_1 = c(tail(fft_f_origin_1, (N-1)%/%2), 
                       head(fft_f_origin_1, N-(N-1)%/%2) )
    fft_f_origin_2 = c(tail(fft_f_origin_2, (N-1)%/%2), 
                       head(fft_f_origin_2, N-(N-1)%/%2) )
    
    l_vec = 0:(N-1)
    l_vec = c(tail(l_vec, (N-1)%/%2)-N,
              head(l_vec, N-(N-1)%/%2) )
    
    
    gd_1 = 2 * sum( Re(1i*2*pi*(l_vec/N)*
                         exp(1i*2*pi*(-n0_vec[[1]])*l_vec/N)*fft_f_origin_1*
                         Conj(fft_f_target-exp(1i*2*pi*(-n0_vec[[2]])*l_vec/N)*fft_f_origin_2) ) ) 
    gd_2 = 2 * sum( Re(1i*2*pi*(l_vec/N)*
                         exp(1i*2*pi*(-n0_vec[[2]])*l_vec/N)*fft_f_origin_2*
                         Conj(fft_f_target-exp(1i*2*pi*(-n0_vec[[1]])*l_vec/N)*fft_f_origin_1) ) ) 
    gd_vec = c(gd_1, gd_2)
    
    return (gd_vec)
  }
