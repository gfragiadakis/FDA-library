#' Cluster rows and columns of a matrix
#' @param data a data matrix
#' @export
#'
cluster_matrix <- function(data){
  rd<-dist(data)
  rc<-hclust(rd)
  cd<-dist(t(data))
  cc<-hclust(cd)
  ordered_data <- dat[rc$order, cc$order]
  cluster_object <- list(data = ordered_data, rc = rc, cc =cc)
  return(cluster_object)
}
