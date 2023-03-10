#' @title Sorting the prototypes by increasing sum of their elements (absolute value)
#'
#' @param gamma A set of prototypes
#'
#' @return The same set of prototypes but sorted by increasing sum of their elements (absolute value)
#' @export
#'
#' @examples
#'
#' sort_gamma(gamma = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3))))
sort_gamma = function(gamma){
  gamma_sorted = gamma
  sums = Vectorize(function(k){sum(gamma[[k]])})(1:length(gamma))
  for(k in 1:length(gamma)){gamma_sorted[[rank(sums)[k]]] = gamma[[k]]}
  return(gamma_sorted)
}
