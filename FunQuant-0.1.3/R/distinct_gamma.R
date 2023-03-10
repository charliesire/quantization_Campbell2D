#' @title Check that all the elements in gamma are distinct
#'
#' @param gamma A list of prototypes
#'
#' @return A boolean indicating whether the elements in gamma are distinct
#' @export
#'
#' @examples
#' distinct_gamma(list(1,2,34,1))
distinct_gamma = function(gamma){
  for(i in 1:(length(gamma)-1)){
    dist_gamma = distance_to_gamma(gamma[[i]], lapply((i+1):length(gamma), function(j){gamma[[j]]}))$dist
    if(dist_gamma == 0){return(FALSE)}
  }
  return(TRUE)
}
