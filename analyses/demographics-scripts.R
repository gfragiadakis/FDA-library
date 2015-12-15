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

# plotting
plot(lasso.pred, y[-train])
lasso.cor <- cor.test(lasso.pred[, 1], y[-train])
lasso.df <- data.frame(predicted = lasso.pred[, 1], actual = y[-train])
p <- ggplot(lasso.df, aes(predicted, actual))
p + geom_point() + geom_abline(intercept = 0, slope = 1) + ggtitle("Age prediction based on immune features")

# extract coefficients
lasso.coefs <- predict(lasso.mod, type = "coefficients", s = lam)
lasso.coefs <- matrix(lasso.coefs, dimnames = list(rownames(lasso.coefs), "coef"))
lasso.coefs[lasso.coefs != 0,]


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

# plotting
pca <- ggplot(PC_df, aes(PC1, PC2))
pca + geom_point(size = 3, aes(colour = Age_Bins)) + scale_color_manual(values=c("orange", "red", "blue")) + ggtitle("PCA by Age (Bins)")
pca + geom_point(aes(colour = Age))
pca + geom_point(size = 3, aes(colour = Gender)) + scale_color_manual(values=c("firebrick", "blue", "gray")) + ggtitle("PCA by Gender")
pca + geom_point(aes(colour = Ethnicity_group)) + theme(panel.background = element_rect(fill = "black"), panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank())
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



