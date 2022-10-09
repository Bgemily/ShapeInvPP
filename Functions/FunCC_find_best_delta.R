FunCC_find_best_delta = function(fun_mat, delta_min, delta_max, num_delta = 10, template.type = "mean", 
                                 theta = 1.5, number = 100, alpha = 0, beta = 0, const_alpha = FALSE, 
                                 const_beta = FALSE, shift.alignement = FALSE, shift.max = 0.1, 
                                 max.iter.align = 100) 
{
  if (length(dim(fun_mat)) != 3) {
    stop("Error: fun_mat should be an array of three dimensions")
  }
  if (template.type == "medoid" & (alpha != 0 | beta != 0)) {
    stop("Error: Medoid template is defined only for alpha and beta equal to 0")
  }
  if (shift.max > 1 | shift.max < 0) {
    stop("Error: shift.max must be in [0,1]")
  }
  if (!(alpha %in% c(0, 1)) | !(beta %in% c(0, 1))) {
    stop("Error: alpha and beta must be 0 or 1")
  }
  if (!(template.type %in% c("mean", "medoid"))) {
    stop(paste0("Error: template.type ", template.type, 
                " is not defined"))
  }
  if (number <= 0 | number%%1 != 0) {
    stop("Error: number must be an integer greater than 0")
  }
  if (!(shift.alignement %in% c(TRUE, FALSE))) {
    stop(paste0("Error: shift.alignement should be a logicol variable"))
  }
  if (max.iter.align <= 0 | max.iter.align%%1 != 0) {
    stop("Error: max.iter.align must be an integer greater than 0")
  }
  if (base::length(dim(fun_mat)) != 3) {
    stop("Error: fun_mat must an array of three dimensions")
  }
  if (!(const_alpha %in% c(FALSE, TRUE)) | !(const_beta %in% 
                                             c(FALSE, TRUE))) {
    stop("Error: const_alpha and const_beta must TRUE or FALSE")
  }
  if (delta_min < 0 | !is.numeric(delta_min) | delta_max < 
      0 | !is.numeric(delta_max)) {
    stop("Error: delta_min and delta_max must be a number greater than 1")
  }
  if (delta_min > delta_max) {
    stop("Error: delta_max must be a number greater than delta_min")
  }
  cl <- delta <- NULL
  delta_check <- seq(delta_min, delta_max, (delta_max - delta_min)/num_delta)
  Htot_best <- delta_max
  best_d <- delta_max
  Htot_sum <- numeric()
  Htot_all_mean <- numeric()
  num_clust <- numeric()
  not_assigned <- numeric()
  for (d in delta_check) {
    res_fun_list <- FunCC:::funcc_biclust(fun_mat, delta = d, template.type = template.type, 
                                  theta = theta, number = number, alpha = alpha, beta = beta, 
                                  const_alpha, const_beta, shift.alignement, shift.max, 
                                  max.iter.align)
    res_fun <- res_fun_list[[1]]
    if (res_fun@Number == 0) {
      Htot_all_mean <- c(Htot_all_mean, NA)
      Htot_sum <- c(Htot_sum, NA)
      num_clust <- c(num_clust, 0)
      not_assigned <- c(not_assigned, nrow(fun_mat) * 
                          ncol(fun_mat))
    }
    if (res_fun@Number == 1) {
      fun_mat_cl <- array(fun_mat[c(res_fun@RowxNumber), 
                                  c(res_fun@NumberxCol), ], dim = c(sum(c(res_fun@RowxNumber)), 
                                                                    sum(c(res_fun@NumberxCol)), dim(fun_mat)[3]))
      dist_mat <- FunCC:::evaluate_mat_dist(fun_mat_cl, template.type, 
                                    alpha, beta, const_alpha, const_beta, shift.alignement, 
                                    shift.max, max.iter.align)
      H_cl <- FunCC:::ccscore_fun(dist_mat)
      elements <- nrow(fun_mat) * ncol(fun_mat)
      elements <- elements - nrow(fun_mat_cl) * ncol(fun_mat_cl)
      not_assigned <- c(not_assigned, elements)
      Htot_d <- mean(H_cl)
      Htot_all_mean <- c(Htot_all_mean, Htot_d)
      Htot_sum <- c(Htot_sum, sum(H_cl))
      num_clust <- c(num_clust, 1)
    }
    if (res_fun@Number > 1) {
      num_clust <- c(num_clust, res_fun@Number)
      H_cl <- numeric()
      elements <- nrow(fun_mat) * ncol(fun_mat)
      for (cl in 1:res_fun@Number) {
        dist_mat <- FunCC:::evaluate_mat_dist(array(fun_mat[c(res_fun@RowxNumber[, 
                                                                         cl]), c(res_fun@NumberxCol[cl, ]), ], dim = c(sum(c(res_fun@RowxNumber[, 
                                                                                                                                                cl])), sum(c(res_fun@NumberxCol[cl, ])), dim(fun_mat)[1])), 
                                      template.type, alpha, beta, const_alpha, const_beta, 
                                      shift.alignement, shift.max, max.iter.align)
        H_cl_temp <- FunCC:::ccscore_fun(dist_mat)
        H_cl <- c(H_cl, H_cl_temp)
        fun_mat_cl <- array(fun_mat[c(res_fun@RowxNumber[, 
                                                         cl]), c(res_fun@NumberxCol[cl, ]), ], dim = c(sum(c(res_fun@RowxNumber[, 
                                                                                                                                cl])), sum(c(res_fun@NumberxCol[cl, ])), dim(fun_mat)[3]))
        elements <- elements - nrow(fun_mat_cl) * ncol(fun_mat_cl)
      }
      not_assigned <- c(not_assigned, elements)
      Htot_d <- mean(H_cl)
      Htot_all_mean <- c(Htot_all_mean, Htot_d)
      Htot_sum <- c(Htot_sum, sum(H_cl))
      if (Htot_d < Htot_best) {
        Htot_best <- Htot_d
        best_d <- d
      }
    }
  }
  h <- data.frame(Htot_sum = Htot_sum, Htot_all_mean = Htot_all_mean, 
                  num_clust = num_clust, delta = delta_check, not_assigned = not_assigned)
  
  return(h)
}
