#### Input Data ####

# Feed in data (df) from compiled_df where you just have Donor, Condition_Feature, value
# and you cast it, then make it where rownames are Donors (Features are colnames)
# preprocess
df <- compiled_df
df <- dplyr::select(df, Donor, Condition_Feature, value)
df <- reshape2::dcast(df, Donor ~ Condition_Feature)
df <- data.frame(df[,-1], row.names = df$Donor)

# now have df:
# generate correlation matrix
correlations <- cor(df)

# function to cluster using correlation between variables as distance (add to FDA library)
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat_reordered <- reorder_cormat(correlations)
adj <- make_adjacency(cormat_reordered, cor_threshold = 0.5)
adj_melted <- reshape2::melt(adj)
write.csv()

# now have adj_melted that has Var1, Var2, value (Immune Feature 1, Immune Feature 2, Adj value)
# which will be the data frame that will go into Shiny for plotting and for adding to with the brushing
# we will add a column to it called "Module" that corresponds to what group Var2 (y-axis) was put into

# to do in Shiny:
# read in csv of data frame
# plot the data
# display the data frame
# have "new module" button-- observe event
# iterate on brushed strokes (user displays the box, the clicks OK)
# export data frame with module numbers

# will use this for plotting
adj_plot <- ggplot(data = adj_melted, aes(Var1, Var2, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste(output_directory, cor_threshold, main_title, "_adjacency_map.pdf", sep = ""), plot = adj_plot, width = 40, height = 40)



# we want to click a "save module" button that would
# 1. save the selected data frame to a file
# 2. plot a box on top of the plot using the xmin, xmax, ymin, ymax

# making a box in ggplot:
p + annotate("rect", xmin = 3, xmax = 4.2, ymin = 12, ymax = 21,
             alpha = .2)
# you can get them by input$blotPrush$xmin, this should work for display

# for the selected frame, can use download handler:
# http://shiny.rstudio.com/articles/download.html, you would get the brushed points with brushedPoints()

# we need an action button, then observe event and then a thing happens (Ex5 soln)


# now shiny is working, we just need to read in the output and add labels back to compiled_df
library(magrittr)
shinyDirectory <- "~/Documents/FDAlibrary/shiny-modules/"
module_files <- list.files(paste(shinyDirectory, "output", sep = ""), pattern = ".csv")
features_all <- c()
module_assignments <- c()
for (i in module_files){
  module_number <- i %>% strsplit(split = "_") %>% .[[1]] %>% .[2] %>%
                  strsplit(split = "[.]") %>% .[[1]] %>% .[1] %>% as.numeric(.)
  print(module_number)
  module_df <- read.csv(paste(shinyDirectory, "output/", i, sep = ""), row.names = 1)
  features <- as.character(unique(module_df$Var1))
  features_all <- c(features_all, features)
  module_assignments <- c(module_assignments, rep(module_number, length(features)))
}
df <- compiled_df
full_module_df <- data.frame(Condition_Feature = features_all, Module = module_assignments)
df$Condition_Feature <- as.factor(convert_names(compiled_df$Condition_Feature))
df_moduled <- merge(df, full_module_df, all = TRUE)
df_moduled$Module[is.na(df_moduled$Module)] <- 0

