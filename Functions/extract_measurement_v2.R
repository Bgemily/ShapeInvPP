
# Extract certain measurement from results of multiple replicates
# Output: dataframe, cols: param_value | measurement(s) | 

extract_measurement_v2 = function(folder_path, param_name=NULL, measurement=c("ARI_mean", "F_mean_sq_err", "v_mean_sq_err")){
  
  measurement_df = data.frame()
  
  param_value_vec = list.files(folder_path)
  for (param_value in param_value_vec) {
    file_name_vec = list.files(path = paste0(folder_path,"/",param_value), full.names = T, recursive = TRUE)
    for (file in file_name_vec) {
      load(file)
      # Size of meas_value_mat: N_meas*N_replicate
      if (FALSE){
        ind = which(measurement=='time_estimation')
        time_est_value_vec = sapply(results, function(one_replicate) tryCatch(as.numeric(one_replicate[["time_estimation"]],
                                                                                      units = 'secs'), 
                                                                          error=function(x)NA))  
        meas_value_mat = sapply(results[sapply(results,is.list)], function(one_replicate) tryCatch(unlist(one_replicate[measurement[-ind]]), 
                                                                      error=function(x)NA)) 
        meas_value_mat = rbind(meas_value_mat, time_est_value_vec)
        
      } else{
        meas_value_mat = sapply(results[sapply(results,is.list)], function(one_replicate) tryCatch(unlist(one_replicate[measurement]), 
                                                                      error=function(x)NA))  
      }
      if (is.vector(meas_value_mat)) 
        meas_value_mat = matrix(meas_value_mat)
      else
        meas_value_mat = t(meas_value_mat)

      # Add column names
      if (ncol(meas_value_mat) == length(measurement)) 
        colnames(meas_value_mat) = measurement
      else if (length(measurement) == 1)
        colnames(meas_value_mat) = rep(measurement,ncol(meas_value_mat))
      
      # Append parameter values as a new column 
      if (!is.null(param_name)){
        param_value = results[sapply(results,is.list)][[1]]$data_param[[param_name]]
        meas_value_df = as.data.frame(cbind(matrix(param_value, byrow = TRUE, 
                                                   nrow = nrow(meas_value_mat), 
                                                   ncol = length(param_value), 
                                                   dimnames = list(NULL,
                                                                   paste('param_value', 1:length(param_value), sep='_'))),
                                            meas_value_mat))
      } else{
        meas_value_df = as.data.frame(cbind("param_value"=as.numeric(param_value), meas_value_mat))
      }
      
      # Append SEEDs
      SEED_vec = rep(0, length(results))
      for (idx_res in 1:length(results)){
        if (is.list(results[[idx_res]])) {
          SEED_vec[idx_res] = results[[idx_res]]$data_param$SEED
        } else{
          SEED_vec[idx_res] = NA
        }
      }
      SEED_vec = SEED_vec[!is.na(SEED_vec)]
      meas_value_df = dplyr::bind_cols(SEED=SEED_vec, meas_value_df)
      
      # Append result to large data frame
      if(nrow(measurement_df)==0 | ncol(measurement_df) == ncol(meas_value_df)){
        measurement_df = dplyr::bind_rows(measurement_df, meas_value_df)
      }
    }
  }
  
  
  return(measurement_df)
}
