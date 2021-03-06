#' Generate and plot correlation and adjacency matrices
#'
#'
#' @param df numeric data frame with rows as observations e.g. individuals and columns are features to correlate
#' @param cor_threshold the threshold for a displayed correlation in the adjacency matrix
#' @param output_directory the path to save the generated plots
#' @param background background color of correlation map: either "black" or "white"
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @export

plot_correlations <- function(df, cor_threshold, output_directory, background = "white", main_title = "basic"){

  # generate correlation matrix
  correlations <- cor(df)

  # function to cluster using correlation between variables as distance
  reorder_cormat <- function(cormat){
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }

  cormat_reordered <- reorder_cormat(correlations)

  adj <- make_adjacency(cormat_reordered, cor_threshold = 0.5)

  cormat_melted <- reshape2::melt(cormat_reordered)

  cor_plot <- ggplot(data = cormat_melted, aes(Var1, Var2, fill = value)) + geom_tile(colour = background) +
    scale_fill_gradient2(low = "red", high = "blue", mid = background, midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, main_title, "_correlation_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)

  adj_melted <- reshape2::melt(adj)

  adj_plot <- ggplot(data = adj_melted, aes(Var1, Var2, fill = value)) + geom_tile() +
    scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, cor_threshold, main_title, "_adjacency_map.pdf", sep = ""), plot = adj_plot, width = 40, height = 40)
}


#' Generate and plot correlation and adjacency matrices between two different data frames
#'
#' @param df1 one of two numeric data frames with rows as observations e.g. individuals and columns are features to correlate
#' @param df2 one of two numeric data frames with rows as observations e.g. individuals and columns are features to correlate
#' @param cor_threshold the threshold for a displayed correlation in the adjacency matrix
#' @param output_directory the path to save the generated plots
#' @param background background color of correlation map: either "black" or "white"
#' @param main_title title for the start of the pdf put out
#' @import ggplot2
#' @import reshape2
#' @import dplyr
#' @export

plot_correlations_asymmetric <- function(df1, df2, cor_threshold = 0.5, output_directory, background = "white", main_title = "basic"){

  # cluster rows and columns
  cormat <- cor(df1, df2)
  cormat_clustered <- cluster_matrix(cormat)
  clust <- cormat_clustered$data
  cormat_melted <- reshape2::melt(clust)

  # plot correlations
  cor_plot <- qplot(x=Var1, y=Var2, data=cormat_melted, fill=value, geom="tile") +
    scale_fill_gradient2(low = "red", high = "blue", mid = background, midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, main_title, "correlation_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)

  adj <- make_adjacency(clust, cor_threshold)
  adj_melted <- reshape2::melt(adj)

  cor_plot <- qplot(x=Var1, y=Var2, data=adj_melted, fill=value, geom="tile") +
    scale_fill_gradient2(low = "red", high = "blue", mid = background, midpoint = 0, limits=c(-1, 1)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(output_directory, main_title, "adj_map.pdf", sep = ""), plot = cor_plot, width = 40, height = 40)
}
