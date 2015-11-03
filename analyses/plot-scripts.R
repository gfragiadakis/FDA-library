library(dplyr)
library(ggplot2)

for (i in levels(asinh_df$Condition)){
  if (i != "Basal"){
    # make this a function (get long format for a given condition)
    df <- asinh_df %>%
      dplyr::filter(Condition == i) %>%
      reshape2::melt(variable.name = "Feature", value.name = "Asinh_ratio")

    p <- ggplot(df, aes(Feature, Asinh_ratio))
    p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(i)
    ggsave(paste(i,"_feature_boxplots.pdf", sep = ""), width = 40, height = 10)
  }

}


compiled_df <- threshold_features(data = asinh_df, threshold = 0.2)

f_df <- dplyr::filter(compiled_df, Gender == "F")
m_df <- dplyr::filter(compiled_df, Gender == "M")

f_grouped <- group_by(f_df, Feature_Condition)
f_sum <- dplyr::summarise(f_grouped, f_mean = mean(value), f_sd = sd(value))

m_grouped <- group_by(m_df, Feature_Condition)
m_sum <- dplyr::summarise(m_grouped, m_mean = mean(value), m_sd = sd(value))

df_gender <- data.frame(m_sum, f_sum)
ids <- strsplit(as.character(df_gender$Feature_Condition), "_")
df_gender <- data.frame(df_gender, Condition = sapply(ids, "[[", 1),
                        Cell_type = sapply(ids, "[[", 2),
                        Protein = sapply(ids, "[[", 3))


p <- ggplot(df_gender, aes(m_mean, f_mean))
p + geom_point() + geom_text(x = 0.25, y = 2, label = "R = 0.9879372") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
ggsave("men_v_women_mean.pdf")

p <- ggplot(df_gender, aes(m_sd, f_sd))
p + geom_point() + geom_text(x = 0.25, y = 0.75, label = "R = 0.9054926") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
ggsave("men_v_women_sd.pdf")

p + geom_point(aes(colour = Condition)) + geom_text(x = 0.25, y = 0.75, label = "R = 0.9054926") + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
ggsave("men_v_women_sd_by_condition.pdf")

p <- ggplot(compiled_df, aes(Feature_Condition, value))
p + geom_boxplot(aes(colour = Condition)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Thresholded features")
ggsave("thresholded_asinh_0.2/thresholded_feature_boxplots_by_Condition.pdf", width = 40, height = 10)



