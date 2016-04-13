# elastic net optimization
# generate fold ids (1 thorugh 10)
# run cv.glmnet (lapply) over alpha (0 to 1 by 0.1)
# do this with a few sets of foldids

## CHANGE EVERYTHING TO RELATIVE PATHS!!
# make a function out of this (in FDAlibrary)

alpha_range <- seq(0, 1, 0.1)
set.seed(1)
fold_sets_1 <- sample(1:10, size = length(train), replace = TRUE)
set.seed(2)
fold_sets_2 <- sample(1:10, size = length(train), replace = TRUE)
set.seed(3)
fold_sets_3 <- sample(1:10, size = length(train), replace = TRUE)
set.seed(1)
fold_sets <- list(fold_sets_1, fold_sets_2, fold_sets_3)

chosen_folds <- fold_sets[[1]]
cv.out <- cv.glmnet(x[train,], y[train], family = "binomial", alpha = 0, foldid = chosen_folds)
lambdas_0 <- cv.out$lambda
cv.out <- cv.glmnet(x[train,], y[train], family = "binomial", alpha = 0.1, foldid = chosen_folds)
lambdas_0.1 <- cv.out$lambda
test_seq <- c(lambdas_0, lambdas_0.1)

#
cv.glmnet.a <- function(a, chosen_folds){
  if(a != 0){
    cv.out <- cv.glmnet(x[train,], y[train], family = "binomial", alpha = a, foldid = chosen_folds)
  }
  else {
    cv.out <- cv.glmnet(x[train,], y[train], family = "binomial", alpha = a, lambda = test_seq, foldid = chosen_folds)

  }
  results_frame <- data.frame(alpha = rep(a, length(cv.out$cvm)),
                              lambda = cv.out$lambda,
                              error = cv.out$cvm)
  return(results_frame)
}


for (i in 1:3){
  # chosen_folds <- fold_sets[[i]]
  error_list <- lapply(alpha_range, cv.glmnet.a, chosen_folds = fold_sets[[i]])
  error_frame <- do.call(rbind, error_list)

  # plot x = lambda, y = error, factor = alpha

  p <- ggplot(error_frame, aes(x = log(lambda), y= error, group=alpha))
  print(p + geom_line(aes(colour = alpha)))
  lowest_error <- error_frame[which.min(error_frame$error), ]
  print(lowest_error)
}

# next ROC curves, sparse grouped lasso
library(ROCR)
rocplot <- function (pred , truth , ...){
  predob <- prediction (pred , truth)
  perf <-  performance (predob , "tpr", "fpr")
  plot(perf ,...)}

pred_probs <- predict(lasso.mod, newx = x[-train,], type = "response", s = lam)
truth <- y[-train]

library(AUC)
roc_data <- roc(pred_probs, truth)
auc(roc_data)

df <- data.frame(False_Positive_Rate = roc_data$fpr, True_Positive_Rate = roc_data$tpr)
ggplot(df, aes(x = False_Positive_Rate, y = True_Positive_Rate)) + geom_line()

# ggplot object for ROC
df <- data.frame(False_Positive_Rate = perf@x.values[[1]], True_Positive_Rate = perf@y.values[[1]])
ggplot(df, aes(x = False_Positive_Rate, y = True_Positive_Rate)) + geom_line()

library(SGL)
# data: list that has $x data matrix $y results
# index: p-vector indicating group membership
# then also have those assignments to delve into the immunology,and then color them as needed as blocks
# (Drew's stuff)

