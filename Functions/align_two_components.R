align_two_components = function(f_target,
                                f_origin_1,
                                f_origin_2,
                                step_size = 0.02,
                                t_unit = 0.05, 
                                n0_vec = c(0,0),
                                n0_min_vec = 0,
                                n0_max_vec = length(f_target), 
                                pad = NULL,
                                periodic = FALSE,
                                MaxIter=1000, 
                                stopping_redu=0.0001, 
                                weights=NULL)
{
  
  ### Extend f1 and f2 from [0,T] to [-T,2T]
  if(!is.null(pad)){
    extend = function(f) {return(c( rep(pad,length(f)-1), f, rep(pad,2*round(length(f)/2)) ))} # make N:=length_of_func odd
  } else{
    extend = function(f) {return(c( rep(head(f,1),length(f)-1), f, rep(tail(f,1),2*round(length(f)/2)) ))} # make N:=length_of_func odd
  }
  f_target = extend(f_target)
  f_origin_1 = extend(f_origin_1)
  f_origin_2 = extend(f_origin_2)
  
  ### Compute terms needed in gradients
  fft_f_target = fft(f_target) #/ length(f_target)
  fft_f_origin_1 = fft(f_origin_1) #/ length(f_origin_1)
  fft_f_origin_2 = fft(f_origin_2) #/ length(f_origin_2)
  
  
  ### Gradient descent
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  converge = FALSE
  while (!converge && iter_count<MaxIter) {
    iter_count = iter_count+1
    
    gd_vec = gradient_two_component(fft_f_target = fft_f_target, 
                                     fft_f_origin_1 = fft_f_origin_1, 
                                     fft_f_origin_2 = fft_f_origin_2,
                                     n0_vec = n0_vec)
    n0_vec = n0_vec - step_size*gd_vec
    # print((step_size)*gd_vec)
    n0_vec = round(n0_vec) 
    
    n0_vec[which(n0_vec<n0_min_vec)] = n0_min_vec[which(n0_vec<n0_min_vec)]
    n0_vec[which(n0_vec>n0_max_vec)] = n0_max_vec[which(n0_vec>n0_max_vec)]

    
    N = length(fft_f_target)
    l_vec = 0:(N-1)
    l_vec = c( head(l_vec, N-(N-1)%/%2),
               tail(l_vec, (N-1)%/%2) - N )
    d = sum(abs( exp(1i*2*pi*l_vec*(-n0_vec[1])/N)*fft_f_origin_1 +
                   exp(1i*2*pi*l_vec*(-n0_vec[2])/N)*fft_f_origin_2 -
                   fft_f_target )^2)  
    d = sqrt(t_unit / N * d)
    dist_upd = d
    dist_redu = (dist_curr-dist_upd)/dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    
    ### DEBUG
    # if(dist_redu<0){
      # print(step_size*gd_vec)
    # }
    
    dist_curr = dist_upd
    converge = dist_redu < stopping_redu
  
  }
  
  
  if (iter_count==MaxIter) {
    warning("Reached max iteration number when estimating a time shift. Consider adjusting the step size.")
  }
  
  return(list(n0_vec=n0_vec))
  
  
}