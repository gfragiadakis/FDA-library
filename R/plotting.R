#' Generate and plot correlation and adjacency matrices
#'
#'
#' @param df numeric data frame with rows as observations e.g. individuals and columns are features to correlate
#' @param cor_threshold the threshold for a displayed correlation in the adjacency matrix
#' @param output_directory the path to save the generated plots
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @export

# df needs to be a numeric data frame whos rows are individuals and columns are features, also plots
# an adjacency matrix
# ggplot2, reshape2, dplyr
plot_correlations <- function(df, cor_threshold, output_directory){

  # generate correlation matrix
  correlations <- cor(df)

  # function to cluster using correlation between variables as distance
  reorder_cormat <- function(cormat){
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }

  cormat_reordered <- reorder_cormat(correlations)

  make_adjacency <- function(ordered_df, cor_threshold){
    adj <- ordered_df
    adj[adj >= cor_threshold] <- 1
    adj[adj < -cor_threshold] <- -1
    adj[(adj > -cor_threshold) & (adj < cor_threshold)] <- 0
    return(adj)
  }

  adj <- make_adjacency(cormat_reordered, cor_threshold = 0.5)

  cormat_melted <- reshape2::melt(cormat_reordered)

  cor_plot <- qplot(x=Var1, y=Var2, data=cormat_melted, fill=value, geom="tile") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, "correlation_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)

  adj_melted <- reshape2::melt(adj)

  adj_plot <- qplot(x=Var1, y=Var2, data = adj_melted, fill=value, geom="tile") +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, cor_threshold,"_adjacency_map.pdf", sep = ""), plot = adj_plot, width = 40, height = 40)
}






# turn this into stuff for package
# 1. separate functions from analysis
# 2. turn them into loaded functions


library(dplyr)
library(magrittr)
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

# For each condition_df, group by Donors to get medians of asinh_ratio
# filter by features that correspond to abs(medians) >= 0.2
# append all conditions together
# then of this can do medians, variance, etc; plot men and women, etc

threshold_features <- function(data, threshold){

  compiled_df <- data.frame()

  for (i in levels(data$Condition)){

    if (i != "Basal"){

      df <- data %>%
            dplyr::filter(Condition == i) %>%
            reshape2::melt(variable.name = "Feature")

      df_grouped <- dplyr::group_by(df, Feature)
      df_stats <- dplyr::summarise(df_grouped, mean = mean(value))

      thresholded_features <- df_stats$Feature[abs(df_stats$mean) >= threshold]

      new_df <- dplyr::filter(df, Feature %in% thresholded_features)
      compiled_df <- rbind(compiled_df, new_df)
    }

  }

  compiled_df <- dplyr::mutate(compiled_df, Condition_Feature = paste(Condition, Feature, sep = "_"))
  compiled_df$Condition_Feature <- as.factor(compiled_df$Condition_Feature)
  return(compiled_df)
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



