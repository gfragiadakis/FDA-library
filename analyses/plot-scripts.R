library(FDAlibrary)
library(dplyr)
library(ggplot2)

###---------Pre-threshold feature plots---------########

plotDirectory <- "~/Documents/FDAlibrary/analyses/plots/"
asinh_df <- read.csv("~/Documents/FDAlibrary/saved_data_structures/human_signaling_fold_asinh.csv", row.names = 1)

for (i in levels(asinh_df$Condition)){
  if (i != "Basal"){
    df <- asinh_df %>%
      dplyr::filter(Condition == i) %>%
      reshape2::melt(variable.name = "Feature", value.name = "Asinh_ratio")

    p <- ggplot(df, aes(Feature, Asinh_ratio))
    p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(i)
    ggsave(paste(plotDirectory, i,"_feature_boxplots.pdf", sep = ""), width = 40, height = 10)
  }

}

# all features in a single plot
df <- reshape2::melt(asinh_df, variable.name = "Feature", value.name = "Asinh_ratio")
df <- dplyr::mutate(df, Condition_Feature = paste(Condition, Feature, sep = "_"))
p <- ggplot(df, aes(Condition_Feature, Asinh_ratio))
p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p + geom_boxplot() + theme(axis.text.x=element_blank())
ggsave(paste(plotDirectory,"ALL_feature_boxplots.pdf", sep = ""), width = 40, height = 10)


###--------Thresholded at 0.2--------#######

compiled_df <- threshold_features(data = asinh_df, threshold = 0.2)
plotDirectory <- "~/Documents/FDAlibrary/analyses/plots/thresholded_asinh_0.2/"

p <- ggplot(compiled_df, aes(Condition_Feature, value))
p + geom_boxplot(aes(colour = Condition)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Thresholded features")
ggsave(paste(plotDirectory, "thresholded_feature_boxplots_by_Condition.pdf", sep = ""), width = 40, height = 10)

p <- ggplot(compiled_df, aes(Condition_Feature, value))
p + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Thresholded features")
ggsave(paste(plotDirectory, "thresholded_feature_boxplots.pdf", sep = ""), width = 40, height = 10)

df <- compiled_df %>%
  dplyr::select(Donor, Condition_Feature, value) %>%
  reshape2::dcast(., Donor ~ Condition_Feature) %>%
  dplyr::select(-Donor)

plot_correlations(df, cor_threshold = 0.5, output_directory = paste(plotDirectory, "correlation/", sep = ""), background = "black", main_title = "v1")

plot_correlations(df, cor_threshold = 0.5, output_directory = paste(plotDirectory, "correlation/", sep = ""), background = "white", main_title = "v2")


####----------Features x Conditions Heatmap ---------#########
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
melted <- melt(ordered_data)

qplot(x=Var2, y=Var1, data=melted, fill=value, geom="tile") +
  scale_fill_gradient2(low = "yellow", high = "purple", mid = "black", midpoint = 0, limits=c(-1, 3)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(plotDirectory, "thresholded_features_heatmap.pdf", sep = ""), width = 40, height = 20)

###------------ Breaking down the correlations by Condition or by Phospho protein-----#########

ids <- strsplit(as.character(compiled_df$Condition_Feature), "_")
df <- data.frame(compiled_df,
                        Cell_type = sapply(ids, "[[", 2),
                        Protein = sapply(ids, "[[", 3))

# IFN map
IFN_df <- df %>%
  dplyr::filter(Condition == "IFNa2" | Condition == "IFNb" | Condition == "IFNg") %>%
  dplyr::select(Donor, Condition_Feature, value) %>%
  reshape2::dcast(., Donor ~ Condition_Feature) %>%
  dplyr::select(-Donor)

plot_correlations(IFN_df, cor_threshold = 0.5, output_directory = paste(plotDirectory, "/correlation/breakdown/",sep = ""), background = "white", main_title = "Interferons")

# Phospho maps

plotlist = list()
for (p in levels(df$Protein)){
  df_p <- df %>%
    dplyr::filter(Protein == p) %>%
    dplyr::select(Donor, Condition_Feature, value) %>%
    reshape2::dcast(., Donor ~ Condition_Feature) %>%
    dplyr::select(-Donor)
 if (ncol(df_p) > 10){
   cormat <- cor(df_p)
   dd <- as.dist((1-cormat)/2)
   hc <- hclust(dd)
   cormat <- cormat[hc$order, hc$order]
   hist(cormat, main = p)
   cormat_melted <- reshape2::melt(cormat)
   cor_plot <- qplot(x=Var1, y=Var2, data=cormat_melted, fill=value, geom="tile") +
     scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
   ggsave(paste(plotDirectory, "/correlation/breakdown/", p, ".pdf", sep = ""))
 plotlist[[p]] <- cor_plot
 }
}

###------------- Clustered heatmap of individuals v. Condition_features--------##########

df <- compiled_df %>%
  dplyr::select(Donor, Condition_Feature, value) %>%
  dcast(., Donor ~ Condition_Feature)

dat <- as.matrix(dplyr::select(df, -Donor))
rownames(dat) <- df$Donor

clustered <- cluster_matrix(dat)
plot(clustered$rc)
plot(clustered$cc)
ordered_data <- clustered$dat
melted <- melt(ordered_data)


qplot(x=Var2, y=Var1, data=melted, fill=value, geom="tile") +
  scale_fill_gradient2(low = "yellow", high = "purple", mid = "white", midpoint = 0, limits=c(-2, 4)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(plotDirectory, "heat_map_indiv.pdf", sep = ""), width = 40, height = 20)


