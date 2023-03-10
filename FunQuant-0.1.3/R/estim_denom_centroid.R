#' Title Compute the estimator which is the denominator of the centroid estimation
#'
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param cell_numbers The output samples that need to be quantized
#' @param cell The cell number of the computed centroid
#' @param bias A number indicating the bias that came out when computing the importance sampling estimators of the membership probabilities of the Voronoi cell. Default is 0.
#'
#' @return A real number which is the denominator of the centroid estimation
#' @export
#' @import abind
#' @examples
#' density_ratio = rep(1,20)
#' cell_numbers = c(1,3,2,1,2,1,1,2,3,3,2,2,2,2,2,3,1,1,3,3)
#' cell = 3
#' estim_denom_centroid(density_ratio = density_ratio, cell_numbers = cell_numbers, cell = cell)
estim_denom_centroid = function(density_ratio, cell_numbers, cell, bias = 0){
  return(sum(density_ratio[as.numeric(cell_numbers) == cell])-bias*length(density_ratio))
}
