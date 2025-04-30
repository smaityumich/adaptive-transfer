library(dplyr)
library(pROC)
library(glmnet)
library(doParallel)
library(PRROC)

get.weights <- function(y, alpha = NULL){
  mean.y = mean(y)
  if(is.null(alpha)){
    alpha = 1
  }
  beta = alpha / (1 + alpha)
  weights = (beta / mean.y) * y + ((1 - beta) / (1 - mean.y)) * (1 - y)
  # weights = rep(1, length(y))
  return(weights)
}

get.metrics <- function(y, pred){
  roc_obj <- suppressMessages(roc(y, as.numeric(pred)))
  auc <- as.numeric(roc_obj$auc)
  pr <- pr.curve(scores.class0 = pred[y == 1], 
                 scores.class1 = pred[y == 0], 
                 curve = TRUE)
  auprc <- as.numeric(pr$auc.integral)
  pred.y <- as.numeric(pred > 0.5)
  cm <- table(as.factor(y), as.factor(pred.y))
  TP <- cm["1", "1"]
  FN <- cm["1", "0"]
  Sensitivity <- TP / (TP + FN)
  
  # Calculate TNR (Specificity) = TN / (TN + FP)
  TN <- cm["0", "0"]
  FP <- cm["0", "1"]
  Specificity <- TN / (TN + FP)
  return(c(auc, auprc, Sensitivity, Specificity))
}

train_test_split <- function(df, train_proportion = 0.75, seed = 67665){
  set.seed(seed)
  train_indices <- sample(1:nrow(df), size = round(train_proportion * nrow(df)))
  train_set <- df[train_indices, ]
  test_set <- df[-train_indices, ]
  return(list(train_set, test_set, train_indices))
}

fit.source.model <- function(eICU, ethnicity = "Caucasian", nfold = 5, seed = 12121)
{
  set.seed(seed)
  if(is.null(ethnicity)){
    eICU.source = eICU %>% select(-ethnicity)
  }
  else{
    eICU.source = eICU %>% filter(ethnicity == "Caucasian") %>% select(-ethnicity)
  }
  
  
  # eICU.source.split = train_test_split(eICU.source, train_proportion = 0.8)
  # eICU.source.train = eICU.source.split[[1]]
  # eICU.source.test = eICU.source.split[[2]]
  
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.source)[, -1]
  y = as.numeric(eICU.source$actualhospitalmortality == "EXPIRED")
  weights = get.weights(y)
  measure = "auc"
  
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  num_cores <- parallel::detectCores()
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  glm.source = cv.glmnet(x, y, weights = weights, foldid = foldid,
                         trace.it = T, nfolds = nfold,
                         type.measure = measure,
                         family = 'binomial',
                         parallel = TRUE, keep = TRUE)
  stopImplicitCluster()
  stopCluster(cl)
  saveRDS(glm.source, file = "~/projects/eICU_TL/glm_source.rds")
}


only.target <- function(eICU.min.train, eICU.min.test, 
                        nfold = 5, verbose = F, seed = 12121, 
                        measure = "auc"){
  set.seed(seed)
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.train)[, -1]
  y = as.numeric(eICU.min.train$actualhospitalmortality == "EXPIRED")
  weights = get.weights(y)
  
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  num_cores <- parallel::detectCores()
  cl <- makeCluster(num_cores)      
  registerDoParallel(cl)      
  glm.min = cv.glmnet(x, y, weights = weights, foldid = foldid,
                      trace.it = F, nfolds = nfold,
                      type.measure = measure,
                      family = 'binomial',
                      parallel = TRUE, keep = TRUE)
  stopImplicitCluster()
  stopCluster(cl)
  
  newx = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.test)[, -1]
  newy = as.numeric(eICU.min.test$actualhospitalmortality == "EXPIRED")
  pred = predict(glm.min, newx = newx, type = "response", s = "lambda.min")
  roc_obj <- suppressMessages(roc(newy, as.numeric(pred)))
  auc <- as.numeric(roc_obj$auc)
  scores = pred
  labels = newy
  pr <- pr.curve(scores.class0 = scores[labels == 1], 
                 scores.class1 = scores[labels == 0], 
                 curve = TRUE)
  
  auprc <- as.numeric(pr$auc.integral)
  
  if(verbose){
    cat("Only target ...")
    cat("AUC", auc, "AUPRC", auprc, "\n\n")
  }
  
  beta = coef(glm.min, s = "lambda.min")
  return(list(c(get.metrics(newy, pred), 0), as.numeric(beta)))
}


only.main.effect <- function(eICU.min.train, eICU.min.test, 
                             glm.source, nfold = 5, verbose = F, seed = 12121, measure = "auc"){
  set.seed(seed)
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.train)[, -1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  x = model.matrix(actualhospitalmortality ~ ., data = eICU.min.train)[, -1]
  y = as.numeric(eICU.min.train$actualhospitalmortality == "EXPIRED")
  weights = get.weights(y)
  
  
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  num_cores <- parallel::detectCores()
  cl <- makeCluster(num_cores)      
  registerDoParallel(cl)      
  glm.min = cv.glmnet(x, y, weights = weights, foldid = foldid,
                      trace.it = F, nfolds = nfold, offset = offset,
                      type.measure = measure,
                      family = 'binomial',
                      parallel = TRUE, keep = TRUE)
  stopImplicitCluster()
  stopCluster(cl)
  
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.test)[, -1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  x = model.matrix(actualhospitalmortality ~ ., data = eICU.min.test)[, -1]
  newy = as.numeric(eICU.min.test$actualhospitalmortality == "EXPIRED")
  
  pred = predict(glm.min, newx = x, newoffset = offset, type = "response", s = "lambda.min")
  roc_obj <- suppressMessages(roc(newy, as.numeric(pred)))
  auc <- as.numeric(roc_obj$auc)
  scores = pred
  labels = newy
  pr <- pr.curve(scores.class0 = scores[labels == 1], 
                 scores.class1 = scores[labels == 0], 
                 curve = TRUE)
  
  auprc <- as.numeric(pr$auc.integral)
  
  if(verbose){
    cat("Two step only main ...")
    cat("AUC", auc, "AUPRC", auprc, "\n\n")
  }
  
  beta = coef(glm.source, s = "lambda.min") 
  beta_D = coef(glm.min, s = "lambda.min")
  beta[1:length(beta_D)]  = beta[1:length(beta_D)] + beta_D
  
  return(list(c(get.metrics(newy, pred), 1), as.numeric(beta)))
}



all.effect <- function(eICU.min.train, eICU.min.test, glm.source,
                       nfold = 5, verbose = F, seed = 12121, measure = "auc"){
  set.seed(seed)
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.train)[, -1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  y = as.numeric(eICU.min.train$actualhospitalmortality == "EXPIRED")
  weights = get.weights(y)
  
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  num_cores <- parallel::detectCores()
  cl <- makeCluster(num_cores)      
  registerDoParallel(cl)      
  glm.min = cv.glmnet(x, y, weights = weights, foldid = foldid,
                      trace.it = F, nfolds = nfold, offset = offset,
                      type.measure = measure,
                      family = 'binomial',
                      parallel = TRUE, keep = TRUE)
  stopImplicitCluster()
  stopCluster(cl)
  
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.test)[, -1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  newy = as.numeric(eICU.min.test$actualhospitalmortality == "EXPIRED")
  
  pred = predict(glm.min, newx = x, newoffset = offset, type = "response", s = "lambda.min")
  roc_obj <- suppressMessages(roc(newy, as.numeric(pred)))
  auc <- as.numeric(roc_obj$auc)
  scores = pred
  labels = newy
  pr <- pr.curve(scores.class0 = scores[labels == 1], 
                 scores.class1 = scores[labels == 0], 
                 curve = TRUE)
  
  auprc <- as.numeric(pr$auc.integral)
  if(verbose){
    cat("Two step only main ...")
    cat("AUC", auc, "AUPRC", auprc, "\n\n")
  }
  beta = coef(glm.source, s = "lambda.min") 
  beta_D = coef(glm.min, s = "lambda.min")
  beta  = beta + beta_D
  
  return(list(c(get.metrics(newy, pred), 1), as.numeric(beta)))
}


adaptive <- function(eICU.min.train, eICU.min.test, glm.source, 
                     nfold = 5, verbose = F, seed = 12121, measure = "auc"){
  set.seed(seed)
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.train)[,-1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  x = cbind(offset, x)
  y = as.numeric(eICU.min.train$actualhospitalmortality == "EXPIRED")
  weights = get.weights(y)
  
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  num_cores <- parallel::detectCores()
  cl <- makeCluster(num_cores)      
  registerDoParallel(cl)   
  penalty_factor <- c(0, rep(1, ncol(x) - 1))
  glm.min = cv.glmnet(x, y, weights = weights, foldid = foldid,
                      trace.it = F, nfolds = nfold,
                      type.measure = measure,
                      family = 'binomial',
                      parallel = TRUE, keep = TRUE, penalty.factor = penalty_factor)
  stopImplicitCluster()
  stopCluster(cl)
  gamma = coef(glm.min, s = "lambda.min")[2]
  
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.test)[, -1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  x = cbind(offset, x)
  newy = as.numeric(eICU.min.test$actualhospitalmortality == "EXPIRED")
  
  pred = predict(glm.min, newx = x, type = "response", s = "lambda.min")
  roc_obj <- suppressMessages(roc(newy, as.numeric(pred)))
  auc <- as.numeric(roc_obj$auc)
  scores = pred
  labels = newy
  pr <- pr.curve(scores.class0 = scores[labels == 1], 
                 scores.class1 = scores[labels == 0], 
                 curve = TRUE)
  
  auprc <- as.numeric(pr$auc.integral)
  if(verbose){
    cat("Two step only main ...")
    cat("AUC", auc, "AUPRC", auprc, "\n\n")
  }
  
  beta = coef(glm.source, s = "lambda.min") 
  beta_D = coef(glm.min, s = "lambda.min")
  beta  = beta + beta_D[-2] * beta_D[2]
  
  return(list(c(get.metrics(newy, pred), gamma), as.numeric(beta)))
}



adaptive.main <- function(eICU.min.train, eICU.min.test, glm.source, 
                     nfold = 5, verbose = F, seed = 12121, measure = "auc"){
  set.seed(seed)
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.train)[,-1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  x = model.matrix(actualhospitalmortality ~ ., data = eICU.min.train)[, -1]
  x = cbind(offset, x)
  y = as.numeric(eICU.min.train$actualhospitalmortality == "EXPIRED")
  weights = get.weights(y)
  
  fold0 <- sample.int(sum(y==0)) %% nfold
  fold1 <- sample.int(sum(y==1)) %% nfold
  foldid <- numeric(length(y))
  foldid[y==0] <- fold0
  foldid[y==1] <- fold1
  foldid <- foldid + 1
  
  num_cores <- parallel::detectCores()
  cl <- makeCluster(num_cores)      
  registerDoParallel(cl)    
  penalty_factor <- c(0, rep(1, ncol(x) - 1))
  glm.min = cv.glmnet(x, y, weights = weights, foldid = foldid,
                      trace.it = F, nfolds = nfold,
                      type.measure = measure,
                      family = 'binomial',
                      parallel = TRUE, keep = TRUE, penalty.factor = penalty_factor)
  stopImplicitCluster()
  stopCluster(cl)
  gamma = coef(glm.min, s = "lambda.min")[2]
  
  x = model.matrix(actualhospitalmortality ~ .^2, data = eICU.min.test)[, -1]
  offset = predict(glm.source, newx = x, response = "link", s = "lambda.min")
  x = model.matrix(actualhospitalmortality ~ ., data = eICU.min.test)[, -1]
  x = cbind(offset, x)
  newy = as.numeric(eICU.min.test$actualhospitalmortality == "EXPIRED")
  
  pred = predict(glm.min, newx = x, type = "response", s = "lambda.min")
  roc_obj <- suppressMessages(roc(newy, as.numeric(pred)))
  auc <- as.numeric(roc_obj$auc)
  scores = pred
  labels = newy
  pr <- pr.curve(scores.class0 = scores[labels == 1], 
                 scores.class1 = scores[labels == 0], 
                 curve = TRUE)
  
  auprc <- as.numeric(pr$auc.integral)
  if(verbose){
    cat("Two step only main ...")
    cat("AUC", auc, "AUPRC", auprc, "\n\n")
  }
  beta = coef(glm.source, s = "lambda.min") 
  beta_D = coef(glm.min, s = "lambda.min")
  beta[1:(length(beta_D)-1)]  = beta[1:(length(beta_D)-1)] + beta_D[-2] * beta_D[2]
  return(list(c(get.metrics(newy, pred), gamma), as.numeric(beta)))
}

