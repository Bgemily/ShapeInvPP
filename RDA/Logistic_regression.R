rm(list=ls())
file_path = "../Code/Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)
library(Hmisc)
library(glmnet)
library(broom)
library(randomForest)

load('../Data/mouse_name.rdata')

method_vec = c('Model4_multi_mouse_initNclus5_split1_gamma0',
               'Model4_multi_mouse_initNclus4_split1_gamma0',
               'Model4_multi_mouse_initNclus6_split1_gamma0')
id_session_considered = c(1:39)
# id_session_considered = 1


### Split one cluster into two
for (method in method_vec){
  ## Response type v.s. percent of neurons in clusters for all (observed) brain areas
  ### Logistic regression, weights = 1
  coef_list_1 = list()
  confusion_matrix_list_1 = list()
  cv_fit_list_1 = list()
  fit_multinom_list_1 = list()
  data_list_1 = list()
  cv_error_vec_1 = c()
  
  
  for (id_session_tmp in id_session_considered){
    mouse_name = mouse_name_vec[id_session_tmp]
    file_path = paste0("../Results/Rdata/RDA",
                       "/", method,
                       "/mouse_",mouse_name,
                       "/pre_stim/Nclus2.Rdata")
    
    load(file_path)
    
    
    res$clusters_allneuron_split_list -> clusters_list
    mem = clus2mem(clusters_list)
    
    df_tmp3 = tibble(mem=mem, 
                     id_node = data_res$id_node_vec, 
                     id_trial = data_res$id_trial_vec,
                     id_session = data_res$id_session_vec,
                     N_spks = data_res$N_spks_vec,
                     response_type = data_res$response_type_vec,
                     brain_area = data_res$brain_area_vec) %>%
      mutate(id_node = as_factor(id_node))
    df_tmp3$response_type[is.na(df_tmp3$response_type)] = 'passive'
    
    df_tmp3.1 = df_tmp3
    df_tmp3.1$response_type = recode(df_tmp3.1$response_type, `-1`='-1', `0`='0/1', `1`='0/1',)
    
    df_tmp3.1.1 = df_tmp3.1
    
    ### Get proportion of cluster sizes for each trial and brain area
    df_tmp3.1 = df_tmp3.1 %>% 
      mutate(response_type=as_factor(response_type),
             #            brain_area=as_factor(brain_area,),
             mem=as_factor(mem)) %>%
      group_by(id_trial, id_session, response_type, brain_area, mem, .drop=F) %>%
      summarise(N_node_pertrialmem=n()) %>%
      group_by(id_trial, id_session, response_type, brain_area, .drop=F) %>%
      mutate(N_node_pertrial=sum(N_node_pertrialmem), 
             prop_node_pertrialmem=N_node_pertrialmem/(N_node_pertrial+.Machine$double.eps)) %>%
      filter(N_node_pertrial>0) %>%
      group_by(response_type, brain_area, mem) %>%
      mutate(N_node_perresptype_mem = sum(N_node_pertrialmem)) %>%
      group_by(response_type, brain_area) %>%
      mutate(prop_node_perresptype_mem = N_node_perresptype_mem/sum(N_node_pertrialmem)) %>%
      group_by(response_type, brain_area, id_session, mem) %>%
      mutate(N_node_perresptype_session_mem = sum(N_node_pertrialmem)) %>%
      group_by(response_type, brain_area, id_session) %>%
      mutate(prop_node_perresptype_session_mem = N_node_perresptype_session_mem/sum(N_node_pertrialmem)) %>%
      mutate(brain_area=factor(brain_area)) %>%
      ungroup()
    df_long = df_tmp3.1
    
    
    
    ### Fit multinomial logistic regression
    df_multinom = df_long %>% 
      select(id_trial, id_session, response_type, brain_area, mem, prop_node_pertrialmem) %>%
      filter(mem != max(levels(df_long$mem)), id_session==id_session_tmp) %>%
      pivot_wider(id_cols = c(id_trial, id_session, response_type), 
                  names_from = c(brain_area, mem), 
                  values_from = prop_node_pertrialmem) %>%
      select(!c(id_trial, id_session))
    
    
    data_list_1[[id_session_tmp]] = df_multinom
    
    ### Standardize explanatory variables
    df_multinom[,-1] = scale(df_multinom[,-1])
    df_multinom = df_multinom[, colSums(is.na(df_multinom)) != nrow(df_multinom)]
    
    
    class_size_vec = ( df_multinom %>% group_by(response_type) %>% mutate(class_size = n()) )$class_size 
    
    
    cv.lasso <- cv.glmnet(x=as.matrix(df_multinom[,-1]), 
                          y=df_multinom$response_type, 
                          alpha = 1, nfolds = nrow(df_multinom[,-1]),
                          family = "multinomial", 
                          weights = rep(1,nrow(df_multinom)),
                          maxit = 5000,
                          type.measure='class')   
    cv_fit_list_1[[id_session_tmp]] = cv.lasso
    cv_error_vec_1[id_session_tmp] = min(cv.lasso$cvm)
    
    model = cv.lasso$glmnet.fit
    fit_multinom_list_1[[id_session_tmp]] = model
    
    
    coef_mat = sapply(1:length(coef(model)), function(i)coef(model)[[i]][,which.min(cv.lasso$cvm)])
    colnames(coef_mat) = names(model$beta)
    rownames(coef_mat)[1] = '(Intercept)'
    coef_list_1[[id_session_tmp]] = coef_mat
    
    
    X = cbind(as.matrix(df_multinom[,-1]))
    prob_mat = exp(cbind(1,X) %*% coef_mat)
    prob_mat = prob_mat / rowSums(prob_mat)
    est_resp = names(model$beta)[apply(prob_mat, 1, which.max)]
    confusion_matrix_tmp = table(true = factor(df_multinom$response_type, levels=c('-1','0/1','passive')), 
                                 est = factor(est_resp, levels=c('-1','0/1','passive')) )
    confusion_matrix_tmp = cbind(confusion_matrix_tmp, 'class.error' = 1-diag(confusion_matrix_tmp)/rowSums(confusion_matrix_tmp))
    confusion_matrix_list_1[[id_session_tmp]] = confusion_matrix_tmp
  }
  
  
  
  
  ## Pre-feedback type v.s. percent of neurons in clusters for all (observed) brain areas
  ### Logistic regression, weights = 1
  coef_list_2 = list()
  confusion_matrix_list_2 = list()
  cv_fit_list_2 = list()
  fit_multinom_list_2 = list()
  data_list_2 = list()
  cv_error_vec_2 = c()
  
  
  for (id_session_tmp in id_session_considered){
    mouse_name = mouse_name_vec[id_session_tmp]
    file_path = paste0("../Results/Rdata/RDA",
                       "/", method,
                       "/mouse_",mouse_name,
                       "/pre_stim/Nclus2.Rdata")
    
    load(file_path)
    
    
    res$clusters_allneuron_split_list -> clusters_list
    mem = clus2mem(clusters_list)
    
    df_tmp3 = tibble(mem=mem, 
                     id_node = data_res$id_node_vec, 
                     id_trial = data_res$id_trial_vec,
                     id_session = data_res$id_session_vec,
                     N_spks = data_res$N_spks_vec,
                     pre_feedback_type = data_res$pre_feedback_type_vec,
                     brain_area = data_res$brain_area_vec) %>%
      mutate(id_node = as_factor(id_node))
    
    df_tmp3.1 = df_tmp3
    
    
    ### Get proportion of cluster sizes for each trial and brain area
    df_tmp3.1 = df_tmp3.1 %>% 
      mutate(pre_feedback_type=as_factor(pre_feedback_type),
             #            brain_area=as_factor(brain_area,),
             mem=as_factor(mem)) %>%
      group_by(id_trial, id_session, pre_feedback_type, brain_area, mem, .drop=F) %>%
      summarise(N_node_pertrialmem=n()) %>%
      group_by(id_trial, id_session, pre_feedback_type, brain_area, .drop=F) %>%
      mutate(N_node_pertrial=sum(N_node_pertrialmem), 
             prop_node_pertrialmem=N_node_pertrialmem/(N_node_pertrial+.Machine$double.eps)) %>%
      filter(N_node_pertrial>0) %>%
      group_by(pre_feedback_type, brain_area, mem) %>%
      mutate(N_node_perresptype_mem = sum(N_node_pertrialmem)) %>%
      group_by(pre_feedback_type, brain_area) %>%
      mutate(prop_node_perresptype_mem = N_node_perresptype_mem/sum(N_node_pertrialmem)) %>%
      group_by(pre_feedback_type, brain_area, id_session, mem) %>%
      mutate(N_node_perresptype_session_mem = sum(N_node_pertrialmem)) %>%
      group_by(pre_feedback_type, brain_area, id_session) %>%
      mutate(prop_node_perresptype_session_mem = N_node_perresptype_session_mem/sum(N_node_pertrialmem)) %>%
      mutate(brain_area=factor(brain_area)) %>%
      ungroup()
    df_long = df_tmp3.1
    
    
    
    ### Fit multinomial logistic regression
    df_multinom = df_long %>% 
      filter( !is.na(pre_feedback_type) ) %>%
      select(id_trial, id_session, pre_feedback_type, brain_area, mem, prop_node_pertrialmem) %>%
      filter(mem != max(levels(df_long$mem)), id_session==id_session_tmp) %>%
      pivot_wider(id_cols = c(id_trial, id_session, pre_feedback_type), 
                  names_from = c(brain_area, mem), 
                  values_from = prop_node_pertrialmem) %>%
      select(!c(id_trial, id_session))
    
    
    data_list_2[[id_session_tmp]] = df_multinom
    
    ### Standardize explanatory variables
    df_multinom[,-1] = scale(df_multinom[,-1])
    df_multinom = df_multinom[, colSums(is.na(df_multinom)) != nrow(df_multinom)]
    
    
    class_size_vec = ( df_multinom %>% group_by(pre_feedback_type) %>% mutate(class_size = n()) )$class_size 
    
    
    cv.lasso <- cv.glmnet(x=as.matrix(df_multinom[,-1]), 
                          y=df_multinom$pre_feedback_type, 
                          alpha = 1, nfolds = nrow(df_multinom[,-1]),
                          family = "multinomial", 
                          weights = rep(1,nrow(df_multinom)),
                          maxit = 5000,
                          type.measure='class')   
    cv_fit_list_2[[id_session_tmp]] = cv.lasso
    cv_error_vec_2[id_session_tmp] = min(cv.lasso$cvm)
    
    model = cv.lasso$glmnet.fit
    fit_multinom_list_2[[id_session_tmp]] = model
    
    
    coef_mat = sapply(1:length(coef(model)), function(i)coef(model)[[i]][,which.min(cv.lasso$cvm)])
    colnames(coef_mat) = names(model$beta)
    rownames(coef_mat)[1] = '(Intercept)'
    coef_list_2[[id_session_tmp]] = coef_mat
    
    
    X = cbind(as.matrix(df_multinom[,-1]))
    prob_mat = exp(cbind(1,X) %*% coef_mat)
    prob_mat = prob_mat / rowSums(prob_mat)
    est_resp = names(model$beta)[apply(prob_mat, 1, which.max)]
    confusion_matrix_tmp = table(true = factor(df_multinom$pre_feedback_type, levels=c('1','-1')), 
                                 est = factor(est_resp, levels=c('1','-1')) )
    confusion_matrix_tmp = cbind(confusion_matrix_tmp, 'class.error' = 1-diag(confusion_matrix_tmp)/rowSums(confusion_matrix_tmp))
    confusion_matrix_list_2[[id_session_tmp]] = confusion_matrix_tmp
  }
  
  
  
  res_list = list(coef_list_1 = coef_list_1,
                  confusion_matrix_list_1 = confusion_matrix_list_1,
                  cv_fit_list_1 = cv_fit_list_1,
                  fit_multinom_list_1 = fit_multinom_list_1,
                  data_list_1 = data_list_1,
                  cv_error_vec_1 = cv_error_vec_1,
                  # 
                  coef_list_2 = coef_list_2,
                  confusion_matrix_list_2 = confusion_matrix_list_2,
                  cv_fit_list_2 = cv_fit_list_2,
                  fit_multinom_list_2 = fit_multinom_list_2,
                  data_list_2 = data_list_2,
                  cv_error_vec_2 = cv_error_vec_2)
  
  
  dir.create(path = paste0("../Results/Rdata/RDA",
                           "/", method,
                           "/Logistic_regression"), recursive = TRUE)
  saveRDS(res_list, file = paste0("../Results/Rdata/RDA",
                                  "/", method,
                                  "/Logistic_regression",
                                  "/Split.RDS")) 
  
  
}


# ### NO split
# for (method in method_vec){
#   ## Response type v.s. percent of neurons in clusters for all (observed) brain areas
#   ### Logistic regression, weights = 1
#   coef_list_1 = list()
#   confusion_matrix_list_1 = list()
#   cv_fit_list_1 = list()
#   fit_multinom_list_1 = list()
#   data_list_1 = list()
#   cv_error_vec_1 = c()
#   
#   
#   for (id_session_tmp in id_session_considered){
#     mouse_name = mouse_name_vec[id_session_tmp]
#     file_path = paste0("../Results/Rdata/RDA",
#                        "/", method,
#                        "/mouse_",mouse_name,
#                        "/pre_stim/Nclus2.Rdata")
#     
#     load(file_path)
#     
#     
#     res$clusters_allneuron_nosplit_list -> clusters_list
#     mem = clus2mem(clusters_list)
#     
#     df_tmp3 = tibble(mem=mem, 
#                      id_node = data_res$id_node_vec, 
#                      id_trial = data_res$id_trial_vec,
#                      id_session = data_res$id_session_vec,
#                      N_spks = data_res$N_spks_vec,
#                      response_type = data_res$response_type_vec,
#                      brain_area = data_res$brain_area_vec) %>%
#       mutate(id_node = as_factor(id_node))
#     df_tmp3$response_type[is.na(df_tmp3$response_type)] = 'passive'
#     
#     df_tmp3.1 = df_tmp3
#     df_tmp3.1$response_type = recode(df_tmp3.1$response_type, `-1`='-1', `0`='0/1', `1`='0/1',)
#     
#     df_tmp3.1.1 = df_tmp3.1
#     
#     ### Get proportion of cluster sizes for each trial and brain area
#     df_tmp3.1 = df_tmp3.1 %>% 
#       mutate(response_type=as_factor(response_type),
#              #            brain_area=as_factor(brain_area,),
#              mem=as_factor(mem)) %>%
#       group_by(id_trial, id_session, response_type, brain_area, mem, .drop=F) %>%
#       summarise(N_node_pertrialmem=n()) %>%
#       group_by(id_trial, id_session, response_type, brain_area, .drop=F) %>%
#       mutate(N_node_pertrial=sum(N_node_pertrialmem), 
#              prop_node_pertrialmem=N_node_pertrialmem/(N_node_pertrial+.Machine$double.eps)) %>%
#       filter(N_node_pertrial>0) %>%
#       group_by(response_type, brain_area, mem) %>%
#       mutate(N_node_perresptype_mem = sum(N_node_pertrialmem)) %>%
#       group_by(response_type, brain_area) %>%
#       mutate(prop_node_perresptype_mem = N_node_perresptype_mem/sum(N_node_pertrialmem)) %>%
#       group_by(response_type, brain_area, id_session, mem) %>%
#       mutate(N_node_perresptype_session_mem = sum(N_node_pertrialmem)) %>%
#       group_by(response_type, brain_area, id_session) %>%
#       mutate(prop_node_perresptype_session_mem = N_node_perresptype_session_mem/sum(N_node_pertrialmem)) %>%
#       mutate(brain_area=factor(brain_area)) %>%
#       ungroup()
#     df_long = df_tmp3.1
#     
#     
#     
#     ### Fit multinomial logistic regression
#     df_multinom = df_long %>% 
#       select(id_trial, id_session, response_type, brain_area, mem, prop_node_pertrialmem) %>%
#       filter(mem != max(levels(df_long$mem)), id_session==id_session_tmp) %>%
#       pivot_wider(id_cols = c(id_trial, id_session, response_type), 
#                   names_from = c(brain_area, mem), 
#                   values_from = prop_node_pertrialmem) %>%
#       select(!c(id_trial, id_session))
#     
#     
#     data_list_1[[id_session_tmp]] = df_multinom
#     
#     ### Standardize explanatory variables
#     df_multinom[,-1] = scale(df_multinom[,-1])
#     df_multinom = df_multinom[, colSums(is.na(df_multinom)) != nrow(df_multinom)]
#     
#     
#     class_size_vec = ( df_multinom %>% group_by(response_type) %>% mutate(class_size = n()) )$class_size 
#     
#     
#     cv.lasso <- cv.glmnet(x=as.matrix(df_multinom[,-1]), 
#                           y=df_multinom$response_type, 
#                           alpha = 1, nfolds = nrow(df_multinom[,-1]),
#                           family = "multinomial", 
#                           weights = rep(1,nrow(df_multinom)),
#                           maxit = 5000,
#                           type.measure='class')   
#     cv_fit_list_1[[id_session_tmp]] = cv.lasso
#     cv_error_vec_1[id_session_tmp] = min(cv.lasso$cvm)
#     
#     model = cv.lasso$glmnet.fit
#     fit_multinom_list_1[[id_session_tmp]] = model
#     
#     
#     coef_mat = sapply(1:length(coef(model)), function(i)coef(model)[[i]][,which.min(cv.lasso$cvm)])
#     colnames(coef_mat) = names(model$beta)
#     rownames(coef_mat)[1] = '(Intercept)'
#     coef_list_1[[id_session_tmp]] = coef_mat
#     
#     
#     X = cbind(as.matrix(df_multinom[,-1]))
#     prob_mat = exp(cbind(1,X) %*% coef_mat)
#     prob_mat = prob_mat / rowSums(prob_mat)
#     est_resp = names(model$beta)[apply(prob_mat, 1, which.max)]
#     confusion_matrix_tmp = table(true = factor(df_multinom$response_type, levels=c('-1','0/1','passive')), 
#                                  est = factor(est_resp, levels=c('-1','0/1','passive')) )
#     confusion_matrix_tmp = cbind(confusion_matrix_tmp, 'class.error' = 1-diag(confusion_matrix_tmp)/rowSums(confusion_matrix_tmp))
#     confusion_matrix_list_1[[id_session_tmp]] = confusion_matrix_tmp
#   }
#   
#   
#   
#   
#   ## Pre-feedback type v.s. percent of neurons in clusters for all (observed) brain areas
#   ### Logistic regression, weights = 1
#   coef_list_2 = list()
#   confusion_matrix_list_2 = list()
#   cv_fit_list_2 = list()
#   fit_multinom_list_2 = list()
#   data_list_2 = list()
#   cv_error_vec_2 = c()
#   
#   
#   for (id_session_tmp in id_session_considered){
#     mouse_name = mouse_name_vec[id_session_tmp]
#     file_path = paste0("../Results/Rdata/RDA",
#                        "/", method,
#                        "/mouse_",mouse_name,
#                        "/pre_stim/Nclus2.Rdata")
#     
#     load(file_path)
#     
#     
#     res$clusters_allneuron_nosplit_list -> clusters_list
#     mem = clus2mem(clusters_list)
#     
#     df_tmp3 = tibble(mem=mem, 
#                      id_node = data_res$id_node_vec, 
#                      id_trial = data_res$id_trial_vec,
#                      id_session = data_res$id_session_vec,
#                      N_spks = data_res$N_spks_vec,
#                      pre_feedback_type = data_res$pre_feedback_type_vec,
#                      brain_area = data_res$brain_area_vec) %>%
#       mutate(id_node = as_factor(id_node))
#     
#     df_tmp3.1 = df_tmp3
#     
#     
#     ### Get proportion of cluster sizes for each trial and brain area
#     df_tmp3.1 = df_tmp3.1 %>% 
#       mutate(pre_feedback_type=as_factor(pre_feedback_type),
#              #            brain_area=as_factor(brain_area,),
#              mem=as_factor(mem)) %>%
#       group_by(id_trial, id_session, pre_feedback_type, brain_area, mem, .drop=F) %>%
#       summarise(N_node_pertrialmem=n()) %>%
#       group_by(id_trial, id_session, pre_feedback_type, brain_area, .drop=F) %>%
#       mutate(N_node_pertrial=sum(N_node_pertrialmem), 
#              prop_node_pertrialmem=N_node_pertrialmem/(N_node_pertrial+.Machine$double.eps)) %>%
#       filter(N_node_pertrial>0) %>%
#       group_by(pre_feedback_type, brain_area, mem) %>%
#       mutate(N_node_perresptype_mem = sum(N_node_pertrialmem)) %>%
#       group_by(pre_feedback_type, brain_area) %>%
#       mutate(prop_node_perresptype_mem = N_node_perresptype_mem/sum(N_node_pertrialmem)) %>%
#       group_by(pre_feedback_type, brain_area, id_session, mem) %>%
#       mutate(N_node_perresptype_session_mem = sum(N_node_pertrialmem)) %>%
#       group_by(pre_feedback_type, brain_area, id_session) %>%
#       mutate(prop_node_perresptype_session_mem = N_node_perresptype_session_mem/sum(N_node_pertrialmem)) %>%
#       mutate(brain_area=factor(brain_area)) %>%
#       ungroup()
#     df_long = df_tmp3.1
#     
#     
#     
#     ### Fit multinomial logistic regression
#     df_multinom = df_long %>% 
#       filter( !is.na(pre_feedback_type) ) %>%
#       select(id_trial, id_session, pre_feedback_type, brain_area, mem, prop_node_pertrialmem) %>%
#       filter(mem != max(levels(df_long$mem)), id_session==id_session_tmp) %>%
#       pivot_wider(id_cols = c(id_trial, id_session, pre_feedback_type), 
#                   names_from = c(brain_area, mem), 
#                   values_from = prop_node_pertrialmem) %>%
#       select(!c(id_trial, id_session))
#     
#     
#     data_list_2[[id_session_tmp]] = df_multinom
#     
#     ### Standardize explanatory variables
#     df_multinom[,-1] = scale(df_multinom[,-1])
#     df_multinom = df_multinom[, colSums(is.na(df_multinom)) != nrow(df_multinom)]
#     
#     
#     class_size_vec = ( df_multinom %>% group_by(pre_feedback_type) %>% mutate(class_size = n()) )$class_size 
#     
#     
#     cv.lasso <- cv.glmnet(x=as.matrix(df_multinom[,-1]), 
#                           y=df_multinom$pre_feedback_type, 
#                           alpha = 1, nfolds = nrow(df_multinom[,-1]),
#                           family = "multinomial", 
#                           weights = rep(1,nrow(df_multinom)),
#                           maxit = 5000,
#                           type.measure='class')   
#     cv_fit_list_2[[id_session_tmp]] = cv.lasso
#     cv_error_vec_2[id_session_tmp] = min(cv.lasso$cvm)
#     
#     model = cv.lasso$glmnet.fit
#     fit_multinom_list_2[[id_session_tmp]] = model
#     
#     
#     coef_mat = sapply(1:length(coef(model)), function(i)coef(model)[[i]][,which.min(cv.lasso$cvm)])
#     colnames(coef_mat) = names(model$beta)
#     rownames(coef_mat)[1] = '(Intercept)'
#     coef_list_2[[id_session_tmp]] = coef_mat
#     
#     
#     X = cbind(as.matrix(df_multinom[,-1]))
#     prob_mat = exp(cbind(1,X) %*% coef_mat)
#     prob_mat = prob_mat / rowSums(prob_mat)
#     est_resp = names(model$beta)[apply(prob_mat, 1, which.max)]
#     confusion_matrix_tmp = table(true = factor(df_multinom$pre_feedback_type, levels=c('1','-1')), 
#                                  est = factor(est_resp, levels=c('1','-1')) )
#     confusion_matrix_tmp = cbind(confusion_matrix_tmp, 'class.error' = 1-diag(confusion_matrix_tmp)/rowSums(confusion_matrix_tmp))
#     confusion_matrix_list_2[[id_session_tmp]] = confusion_matrix_tmp
#   }
#   
#   
#   
#   res_list = list(coef_list_1 = coef_list_1,
#                   confusion_matrix_list_1 = confusion_matrix_list_1,
#                   cv_fit_list_1 = cv_fit_list_1,
#                   fit_multinom_list_1 = fit_multinom_list_1,
#                   data_list_1 = data_list_1,
#                   cv_error_vec_1 = cv_error_vec_1,
#                   # 
#                   coef_list_2 = coef_list_2,
#                   confusion_matrix_list_2 = confusion_matrix_list_2,
#                   cv_fit_list_2 = cv_fit_list_2,
#                   fit_multinom_list_2 = fit_multinom_list_2,
#                   data_list_2 = data_list_2,
#                   cv_error_vec_2 = cv_error_vec_2)
#   
#   dir.create(path = paste0("../Results/Rdata/RDA",
#                            "/", method,
#                            "/Logistic_regression"), recursive = TRUE)
#   saveRDS(res_list, file = paste0("../Results/Rdata/RDA",
#                                   "/", method,
#                                   "/Logistic_regression",
#                                   "/No_split.RDS")) 
#   
#   
# }



