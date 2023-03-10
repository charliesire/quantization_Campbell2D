#' @title Compute the distance between a point and its nearest centroid, returning this distance and the associated cell number
#'
#' @param x A point in the output space
#' @param gamma A set of prototypes
#' @param distance_func A function computing a distance between two elements in the output spaces
#'
#' @return The distance between a point and its nearest centroid
#' @export
#'
#' @examples
#' distance_to_gamma(array(1:9, dim = c(3,3)), list(array(10, dim = c(3,3)),
#' array(5, dim = c(3,3)), array(6, dim = c(3,3))),
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))})
distance_to_gamma = function(x, gamma, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}){
  distance = Vectorize(function(k){distance_func(x, gamma[[k]])})(1:length(gamma))
  return(list(cellule = which.min(distance), dist = min(distance)))
}
