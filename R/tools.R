#' Cluster rows and columns of a matrix
#' @param data a data matrix
#' @export
#'
cluster_matrix <- function(data){
  rd<-dist(data)
  rc<-hclust(rd)
  cd<-dist(t(data))
  cc<-hclust(cd)
  ordered_data <- data[rc$order, cc$order]
  cluster_object <- list(data = ordered_data, rc = rc, cc =cc)
  return(cluster_object)
}

#' Make an adjacency matrix
#'
#' @param correlation_matrix a matrix of correlated values
#' @param cor_threshold the significance cutoff for the adjacency matrix
#' @export

make_adjacency <- function(correlation_matrix, cor_threshold = 0.5){
  adj <- correlation_matrix
  adj[adj >= cor_threshold] <- 1
  adj[adj < -cor_threshold] <- -1
  adj[(adj > -cor_threshold) & (adj < cor_threshold)] <- 0
  return(adj)
}

