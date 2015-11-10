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