library(dplyr)
library(reshape2)
library(glmnet)
library(ggplot2)

dataDirectory <- "~/Documents/FDAlibrary/saved_data_structures/"
demographics <- read.csv(paste(dataDirectory, "demographics.csv", sep = ""))
compiled_df <- read.csv(paste(dataDirectory, "human_signaling_fold_asinh_0.2.csv", sep = ""), row.names = 1)

demographics <- demographics[(demographics$Donor %in% compiled_df$Donor),]

df <- dplyr::select(compiled_df, Donor, Gender, Condition_Feature, value)
wide_df <- dcast(df, Donor + Gender ~ Condition_Feature)

#####--------Age Analysis-------###########

dems <- dplyr::select(demographics, Donor, Age)
data <- merge(wide_df, dems)

# Setting up training and test sets
set.seed(1)
train <- sample(1:86, 50, replace = FALSE)
x <- dplyr::select(data, -Age, -Donor, -Gender)
x <- as.matrix(x)
y <- data$Age

# training lasso model
grid <- 10^seq(10,-2, length =100)
lasso.mod <- glmnet(x[train,], y[train], alpha=1, lambda =grid)

# use cross validation to choose lambda
cv.out <- cv.glmnet(x[train,], y[train], alpha = 1)
plot(cv.out)
lam <- cv.out$lambda.1se

# prediction
lasso.pred <- predict(lasso.mod, s = lam, newx = x[-train,])
lasso.pred.train <- predict(lasso.mod, s = lam, newx = x[train,])

# training error
plot(lasso.pred.train, y[train])
lasso.cor.train <- cor.test(lasso.pred.train[, 1], y[train])
# training R value is 0.831
lasso.df.train <- data.frame(predicted_training = lasso.pred.train[, 1], actual = y[train])
p <- ggplot(lasso.df.train, aes(predicted_training, actual))
p + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Age prediction (training set)") + coord_fixed(ratio = 1, xlim = c(20,65), ylim = c(15, 65))


# plotting
plot(lasso.pred, y[-train])
lasso.cor <- cor.test(lasso.pred[, 1], y[-train])
# test R value is 0.64
lasso.df <- data.frame(predicted = lasso.pred[, 1], actual = y[-train])
p <- ggplot(lasso.df, aes(predicted, actual))
p + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Age prediction based on immune features") + coord_fixed(ratio = 1, xlim = c(20,65), ylim = c(15, 65))

# extract coefficients
lasso.coefs <- predict(lasso.mod, type = "coefficients", s = lam)
lasso.coefs <- matrix(lasso.coefs, dimnames = list(rownames(lasso.coefs), "coef"))
lasso.coefs[lasso.coefs != 0,]

# training ridge model
grid <- 10^seq(10,-2, length =100)
ridge.mod <- glmnet(x[train,], y[train], alpha=0, lambda =grid)

# use cross validation to choose lambda
cv.ridge <- cv.glmnet(x[train,], y[train], alpha = 0, lambda = grid)
plot(cv.ridge)
lam <- cv.ridge$lambda.1se

# prediction
ridge.pred <- predict(ridge.mod, s = lam, newx = x[-train,])
ridge.pred.train <- predict(ridge.mod, s = lam, newx = x[train,])

# plotting
plot(ridge.pred, y[-train])
ridge.cor <- cor.test(lasso.pred[, 1], y[-train])
#
ridge.df <- data.frame(predicted = ridge.pred[, 1], actual = y[-train])
p <- ggplot(ridge.df, aes(predicted, actual))
p + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Age prediction based on immune features") + coord_fixed(ratio = 1, xlim = c(20,65), ylim = c(15, 65))

#----age with grouping info:

df <- df_moduled %>%
  dplyr::select(Condition_Feature, Module) %>%
  unique()
data_names <- convert_names(colnames(x))
module_names <- as.character(df$Condition_Feature)
identical(data_names, module_names)

index <- df$Module

index <- df$Module
cv.out.SGL <- cvSGL(data = list(x = x[train, ], y = y[train]), index = index, type = "linear")
plot(cv.out.SGL)
which.min(cv.out.SGL$lldiff)
SGLfit <- SGL(data = list(x = x[train, ], y = y[train]), index = index, type = "linear")
predictions.train <- predictSGL(SGLfit, newX = x[train, ])

predictions <- predictSGL(SGLfit, newX = x[-train, ])
cor(y[-train][-34], predictions[-34, 14])
# get 64% accuracy with grouping info!

SGL_coeffs <- data.frame(df, SGL_coefficients = SGLfit$beta[,14])
# down to 6 modules, 19 coefficients
# stat signaling (could site Mark's stuff) and inflammatory signaling (CREB, P38, ERK)


#----Age binning----#

bins_3 <- rep(NA, length = nrow(dems))
bins_3[dems$Age <= 30] <- "30 and under"
bins_3[(dems$Age > 30 & dems$Age < 50)] <- "btwn 30 and 50"
bins_3[dems$Age >= 50] <- "50+"

#------SAM on age------#
x_age <- x[-2, ]
y_age <- dems$Age[-2]
model <- SAM(t(x_age), y = y_age, resp.type = "Quantitative", genenames = colnames(x), nperms = 1000, fdr.output = 0.05)
# regression yielded nothing

# Bins: SAM multiclass (y  takes 1, 2, 3, ..., k where k is number of classes)
y <- as.numeric(as.factor(bins_3))
y_age <- y[-2]
model <- SAM(t(x_age), y = y_age, resp.type = "Multiclass", genenames = colnames(x), nperms = 1000, fdr.output = 0.05)
# nothing

#----Ethnicity Binning---#
# African
# Asian or Chinese
# Hispanic
# Native (note one Hispanic will get reassigned here)
# else White

ethnicity_group <- rep(NA, length = nrow(demographics))
ethnicity_group[grepl("African", demographics$Ethnicity)] <- "African American"
ethnicity_group[grepl("Asian", demographics$Ethnicity)] <- "Asian"
ethnicity_group[grepl("Chinese", demographics$Ethnicity)] <- "Asian"
ethnicity_group[grepl("Hispanic", demographics$Ethnicity)] <- "Hispanic"
ethnicity_group[grepl("Native", demographics$Ethnicity)] <- "Native"
ethnicity_group[2] <- "Not Reported"
ethnicity_group[is.na(ethnicity_group)] <- "Caucasian"



###--------PCA on Age, Gender, and Ethnicity------###
pr.out <- prcomp(x, scale = FALSE)
transformed_x <- pr.out$x
PC_df <- data.frame(transformed_x, Age_Bins = bins_3, Age = dems$Age, Gender = data$Gender, Ethnicity = demographics$Ethnicity,
                    Ethnicity_group = as.factor(ethnicity_group))

# loadings
loadings1 <- pr.out$rotation[, "PC1"]
loadings2 <- pr.out$rotation[, "PC2"]
loadings3 <- pr.out$rotation[, "PC3"]


# plotting
pca <- ggplot(PC_df, aes(PC1, PC2))
pca + geom_point(size = 3, aes(colour = Age_Bins)) + scale_color_manual(values=c("orange", "red", "blue")) + ggtitle("PCA by Age (Bins)")
pca + geom_point(size = 2, aes(colour = Age)) + ggtitle("PCA by Age")
pca + geom_point(size = 2, aes(colour = Gender)) + scale_color_manual(values=c("firebrick", "blue", "gray")) + ggtitle("PCA by Gender")
pca + geom_point(size = 2, aes(colour = Ethnicity_group)) + ggtitle("PCA by Ethnicity")
# + theme(panel.background = element_rect(fill = "black"), panel.grid.major = element_blank(),
#                                                        panel.grid.minor = element_blank())
# histogram of ethnicities
ggplot(PC_df, aes(Ethnicity_group)) + geom_bar()

##### ---------- Sex Model ---------#######

# training and test sets
set.seed(1)
train <- sample(1:85, 50, replace = FALSE)
x <- dplyr::select(data, -Age, -Donor, -Gender)
x <- as.matrix(x)
y <- data$Gender
x <- x[y!="x", ]
y <- y[y!="x"]
y <- droplevels(y)

# training lasso model
grid <- 10^seq(10,-2, length =100)
lasso.mod <- glmnet(x[train,], y[train], family = "binomial", alpha=1, lambda =grid)

# use cross validation to choose lambda
cv.out <- cv.glmnet(x[train,], y[train], family = "binomial", alpha = 1)
plot(cv.out)
lam <- cv.out$lambda.1se

# prediction
lasso.pred <- predict(lasso.mod, newx = x[-train,], type = "class", s = lam)
pred <- as.factor(lasso.pred[, 1])
table(pred, y[-train])
# 60% correct classification


# extract coefficients
lasso.coefs <- predict(lasso.mod, type = "coefficients", s = lam)
lasso.coefs <- matrix(lasso.coefs, dimnames = list(rownames(lasso.coefs), "coef"))
lasso.coefs[lasso.coefs != 0,]



