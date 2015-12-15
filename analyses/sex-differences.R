library(FDAlibrary)
library(dplyr)
library(ggplot2)

plotDirectory <- "~/Documents/FDAlibrary/analyses/plots/thresholded_asinh_0.2/sex_differences/"

###----------Comparing means and standard deviations of men v. women features--------------###

f_df <- dplyr::filter(compiled_df, Gender == "F")
m_df <- dplyr::filter(compiled_df, Gender == "M")

f_grouped <- group_by(f_df, Condition_Feature)
f_sum <- dplyr::summarise(f_grouped, f_mean = mean(value), f_sd = sd(value))

m_grouped <- group_by(m_df, Condition_Feature)
m_sum <- dplyr::summarise(m_grouped, m_mean = mean(value), m_sd = sd(value))

df_gender <- data.frame(m_sum, f_sum)
ids <- strsplit(as.character(df_gender$Condition_Feature), "_")
df_gender <- data.frame(df_gender, Condition = sapply(ids, "[[", 1),
                        Cell_type = sapply(ids, "[[", 2),
                        Protein = sapply(ids, "[[", 3))

p <- ggplot(df_gender, aes(m_mean, f_mean))
p + geom_point() + geom_text(x = 0.25, y = 2, label = paste("R = ", round(cor(df_gender$m_mean, df_gender$f_mean),digits = 3), sep = "")) + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
ggsave(paste(plotDirectory, "men_v_women_mean.pdf", sep = ""))

p <- ggplot(df_gender, aes(m_sd, f_sd))
p + geom_point() + geom_text(x = 0.25, y = 0.75, label = paste("R = ", round(cor(df_gender$m_sd, df_gender$f_sd),digits = 3), sep = "")) + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
ggsave(paste(plotDirectory, "men_v_women_sd.pdf", sep = ""))

p + geom_point(aes(colour = Condition)) + geom_text(x = 0.25, y = 0.75, label = paste("R = ", round(cor(df_gender$m_sd, df_gender$f_sd),digits = 3), sep = "")) + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
ggsave(paste(plotDirectory, "men_v_women_sd_by_condition.pdf", sep = ""))

###----------SAM comparison of Men v. Women features--------------###

library(samr)
library(reshape2)
df <- dplyr::select(compiled_df, Donor, Gender, Condition_Feature, value)
wide_df <- dcast(df, Donor + Gender ~ Condition_Feature)
x <- as.matrix(dplyr::select(wide_df, -Donor, -Gender))
y <- rep(NA, length(wide_df$Gender))
y[wide_df$Gender == "F"] <- 2
y[wide_df$Gender != "F"] <- 1

model <- SAM(t(x), y, resp.type = "Two class unpaired", genenames = colnames(x), nperms = 1000, fdr.output = 0.05)

# q val 0
up_genes <- model$siggenes.table$genes.up[model$siggenes.table$genes.up[,"q-value(%)"]==0, "Gene ID"]
lo_genes <- model$siggenes.table$genes.lo[model$siggenes.table$genes.lo[,"q-value(%)"]==0, "Gene ID"]
write.csv(up_genes, file = paste(plotDirectory, "SAM_features_women.csv", sep = ""))
write.csv(lo_genes, file = paste(plotDirectory, "SAM_features_men.csv", sep = ""))

# Plot significant features as boxplots

df <- dplyr::filter(df, Gender != "x")
df_up <- df[df[, "Condition_Feature"] %in% up_genes,]
ggplot(df_up, aes(Condition_Feature, value)) + geom_boxplot(aes(colour = Gender)) + scale_color_manual(values=c("firebrick1", "royalblue"))
ggsave(filename = paste(plotDirectory,"gender_sig_up.pdf", sep = ""), width = 30, height = 10)

df_lo <- df[df[, "Condition_Feature"] %in% lo_genes,]
ggplot(df_lo, aes(Condition_Feature, value)) + geom_boxplot(aes(colour = Gender)) + scale_color_manual(values=c("firebrick1", "royalblue"))
ggsave(filename = paste(plotDirectory,"gender_sig_lo.pdf", sep = ""), width = 30, height = 10)

#####--------------Correlation differences----------------------- #####

m_df <- compiled_df %>%
  dplyr::filter(Gender == "M") %>%
  dplyr::select(Donor, Condition_Feature, value) %>%
  reshape2::dcast(., Donor ~ Condition_Feature) %>%
  dplyr::select(-Donor)

f_df <- compiled_df %>%
  dplyr::filter(Gender == "F") %>%
  dplyr::select(Donor, Condition_Feature, value) %>%
  reshape2::dcast(., Donor ~ Condition_Feature) %>%
  dplyr::select(-Donor)

plot_correlations(f_df, cor_threshold = 0.5, output_directory = plotDirectory, background = "black", main_title = "Female Structure")

# Male as female ordering
f_cormat <- cor(f_df)
dd <- as.dist((1-f_cormat)/2)
hc <- hclust(dd)

m_cormat <- cor(m_df)
m_cormat <- m_cormat[hc$order, hc$order]

cormat_melted <- reshape2::melt(m_cormat)

cor_plot <- ggplot(data = cormat_melted, aes(Var1, Var2, fill = value)) + geom_tile(colour = "black") +
  scale_fill_gradient2(low = "red", high = "blue", mid = "black", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(plotDirectory, "Male_correlation_map.pdf", sep = ""), plot = cor_plot,width = 40, height = 40)

# difference between correlations
f_cormat <- cor(f_df)
dd <- as.dist((1-f_cormat)/2)
hc <- hclust(dd)
f_cormat <- f_cormat[hc$order, hc$order]
m_cormat <- cor(m_df)
m_cormat <- m_cormat[hc$order, hc$order]

cormat_diff <- f_cormat - m_cormat
cormat_melted <- reshape2::melt(cormat_diff)

cor_plot <- ggplot(data = cormat_melted, aes(Var1, Var2, fill = value)) + geom_tile(colour = "white") +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(plotDirectory, "Difference_correlation_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)
