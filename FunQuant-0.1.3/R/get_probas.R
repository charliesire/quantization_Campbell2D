#' @title Computing the probability masses of each voronoi cells
#'
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell. Default is "unique".
#' @param cell_numbers The voronoi cell number of every output. If NULL, then outputs, gamma and distance_func must be provided.
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma A set of prototypes. Useful only if cell_numbers == NULL. If NULL, "cells" must be provided.
#' @param distance_func A function computing a distance between two elements in the output spaces. Useful only if cell_numbers == NULL.
#' @param cells The cell numbers that are investigated
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoi cell.

#' @return A vector providing the probability masses of each Voronoï cell.
#'
#' @export
#' @import abind
#' @examples
#' gamma = list(array(10, dim = c(3,3)), array(5, dim = c(3,3)), array(6, dim = c(3,3)))
#' outputs = array(runif(9*20)*20, dim = c(3,3,20))
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' density_ratio = rep(1,20)
#' get_probas(density_ratio = density_ratio, outputs = outputs,
#' gamma = gamma, distance_func = distance_func)
get_probas = function(density_ratio, method_IS = "unique", cell_numbers = NULL, outputs = NULL, gamma = NULL, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}, cells = NULL, bias = NULL){
  if(is.null(cell_numbers)){
    if(method_IS == "unique"){
      cell_numbers = get_cell_numbers(outputs = outputs, gamma = gamma, distance_func = distance_func)
    }
    else if(method_IS == "percell"){
      cell_numbers = lapply(cells, function(k){get_cell_numbers(outputs = outputs[[k]], gamma = gamma, distance_func = distance_func)})
    }
  }
  if(is.null(cells)){cells = 1:length(gamma)}
  if(is.null(bias)){bias = rep(0,length(cells))}
  probas = c()
  for(cell in cells){
    if(method_IS == "unique"){
      probas = c(probas, sum(density_ratio[cell_numbers == cell])/length(density_ratio)-bias[cell])
    }
    else if(method_IS == "percell"){
      probas = c(probas, sum(density_ratio[cell_numbers[[cell]] == cell])/length(density_ratio[[cell]])-bias[cell])
    }

  }

  return(probas)
}
