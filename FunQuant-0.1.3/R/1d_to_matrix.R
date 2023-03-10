#' @title Transform every 1D prototype into a column matrix
#'
#' @param gamma a list of prototypes
#'
#' @return A list of matrices prototypes
#' @export
#'
#' @examples
#' oned_to_matrix(list(1:5, runif(5), rep(0,5)))
oned_to_matrix = function(gamma){
  gamma = lapply(1:length(gamma), function(j){
    if(is.null(dim(gamma[[j]])) | length(dim(gamma[[j]])) == 1){t(as.matrix(gamma[[j]]))}
    else{gamma[[j]]}})
  return(gamma)
}
