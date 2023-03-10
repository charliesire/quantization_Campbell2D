#' @title Compute the centroid and the probability mass of the Voronoï cells
#'
#' @param outputs The output samples that need to be quantized. If method = "percell", a list of output samples must be provided, of length equal to the number of Voronoï cells.
#' @param cell_numbers The voronoi cell number of every output
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell.
#' @param density_ratio A vector indicating the weight fX/g of each output. Default is a vector of 1. If method = "percell", a list of density_ratio must be provided, of length equal to the number of Voronoï cells.
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoi cell. Default is 0 for all Voronoi cells.

#' @return The centroid and the probability mass of each probability cell
#' @export
#' @import abind
#' @examples
#' outputs = array(runif(9*20)*15, dim = c(3,3,20))
#' cell_numbers = c(1,3,2,1,2,1,1,2,3,3,2,2,2,2,2,3,1,1,3,3)
#' density_ratio = rep(1,20)
#' compute_centroids_and_proba(outputs = outputs,cell_numbers = cell_numbers,
#' density_ratio = density_ratio)

compute_centroids_and_proba = function(outputs, cell_numbers, method_IS = "unique", density_ratio = rep(1, dim(outputs)[length(dim(outputs))]), bias = rep(0,length(unique(unlist(cell_numbers)))), batch = FALSE){
  n = length(unlist(cell_numbers))#nb of outputs
  nb_cells = length(unique(unlist(cell_numbers)))
  centroids = list()
  probas = c()
  for(j in 1:nb_cells){
    if(method_IS == "unique"){
      outputs_j = outputs
      cell_numbers_j = cell_numbers
      density_j = density_ratio
    }
    else if(method_IS == "percell"){
      outputs_j = outputs[[j]]
      cell_numbers_j = cell_numbers[[j]]
      density_j = density_ratio[[j]]
    }
    numerator =  estim_num_centroid(outputs = outputs_j, cell_numbers = cell_numbers_j, density_ratio = density_j, cell = j, batch = batch) ## Sum the Y(X)f/nu of the cell
    if(batch){denominator = sum(Vectorize(function(p){estim_denom_centroid(density_ratio = density_j[[p]], cell_numbers = cell_numbers_j[[p]], cell = j, bias = bias[j])})(1:length(density_j)))}
    else{denominator = estim_denom_centroid(density_ratio = density_j, cell_numbers = cell_numbers_j, cell = j, bias = bias[j])}
    centroids[[j]] = numerator/denominator
    probas = c(probas, denominator/n)
  }
  return(list(centroids = centroids, probas = probas))
}
