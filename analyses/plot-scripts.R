library(dplyr)
library(ggplot2)

plotDirectory <- "~/Documents/FDAlibrary/analyses/plots/"

for (i in levels(asinh_df$Condition)){
  if (i != "Basal"){
    # make this a function (get long format for a given condition)
    df <- asinh_df %>%
      dplyr::filter(Condition == i) %>%
      reshape2::melt(variable.name = "Feature", value.name = "Asinh_ratio")

    p <- ggplot(df, aes(Feature, Asinh_ratio))
    p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(i)
    ggsave(paste(plotDirectory, i,"_feature_boxplots.pdf", sep = ""), width = 40, height = 10)
  }

}


compiled_df <- threshold_features(data = asinh_df, threshold = 0.2)

plotDirectory <- "~/Documents/FDAlibrary/analyses/plots/thresholded_asinh_0.2/"

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

p <- ggplot(compiled_df, aes(Condition_Feature, value))
p + geom_boxplot(aes(colour = Condition)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Thresholded features")
ggsave(paste(plotDirectory, "thresholded_feature_boxplots_by_Condition.pdf", sep = ""), width = 40, height = 10)

p <- ggplot(compiled_df, aes(Condition_Feature, value))
p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Thresholded features")
ggsave(paste(plotDirectory, "thresholded_feature_boxplots.pdf", sep = ""), width = 40, height = 10)


##### Correlation #####
df <- compiled_df %>%
      dplyr::select(Donor, Condition_Feature, value) %>%
      reshape2::dcast(., Donor ~ Condition_Feature) %>%
      dplyr::select(-Donor)

plot_correlations(df, cor_threshold = 0.5, output_directory = paste(plotDirectory, "correlation/", sep = ""), background = "black")

plot_correlations(df, cor_threshold = 0.5, output_directory = paste(plotDirectory, "correlation/", sep = ""), background = "white")


### Looking at stims ###
df <- compiled_df %>%
  dplyr::select(Donor, Condition_Feature, value) %>%
  group_by(Condition_Feature) %>%
  dplyr::summarise(mean = mean(value))

feats <- compiled_df %>%
    dplyr::select(Condition, Feature, Condition_Feature) %>%
    .[!duplicated(.),]

new_df <- df %>%
  merge(., feats) %>%
  dplyr::select(Condition, Feature, mean) %>%
  dcast(., Condition ~ Feature)


new_df[is.na(new_df)] <- 0
dat <- as.matrix(dplyr::select(new_df, -Condition))
rownames(dat) <- new_df$Condition

clustered <- cluster_matrix(dat)
plot(clustered$rc)
plot(clustered$cc)
ordered_data <- clustered$data

qplot(x=Var2, y=Var1, data=melted, fill=value, geom="tile") +
  scale_fill_gradient2(low = "yellow", high = "purple", mid = "white", midpoint = 0, limits=c(-1, 3)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(plotDirectory, "heat_map3.pdf", sep = ""), width = 40, height = 20)



