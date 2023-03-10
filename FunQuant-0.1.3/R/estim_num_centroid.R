#' @title Compute the estimator which is the numerator of the centroid estimation
#'
#' @param outputs The output samples that need to be quantized
#' @param cell_numbers The voronoi cell number of every output
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param cell The cell number of the computed centroid
#'
#' @return An array having the same dimension as an output, which is the numerator of the centroid estimator
#' @export
#' @import abind
#' @examples
#' outputs = array(runif(9*20)*15, dim = c(3,3,20))
#' cell_numbers = c(1,3,2,1,2,1,1,2,3,3,2,2,2,2,2,3,1,1,3,3)
#' density_ratio = rep(1,20)
#' cell = 3
#' estim_num_centroid(outputs = outputs,cell_numbers = cell_numbers,
#' density_ratio = density_ratio, cell = cell)
estim_num_centroid = function(outputs, cell_numbers, density_ratio, cell, batch = FALSE){
  if(batch){
    res = 0
    for(batch_i in 1:length(outputs)){
      outputs_cell = asub(x = outputs[[batch_i]], dims = length(dim(outputs[[batch_i]])), idx = which(cell_numbers[[batch_i]] == cell))
      outputs_cell = matrix(outputs_cell, nrow = prod(dim(outputs[[batch_i]])[1:(length(dim(outputs[[batch_i]]))-1)]))
      outputs_cell = t(outputs_cell)
      outputs_cell = outputs_cell*density_ratio[[batch_i]][cell_numbers[[batch_i]] == cell]
      res = res + apply(outputs_cell,2,sum)
    }
    res = array(res, dim = dim(outputs[[1]])[1:(length(dim(outputs[[1]]))-1)])
    
  }
  else{
  outputs_cell = asub(x = outputs, dims = length(dim(outputs)), idx = which(cell_numbers == cell))
  outputs_cell = matrix(outputs_cell, nrow = prod(dim(outputs)[1:(length(dim(outputs))-1)]))
  outputs_cell = t(outputs_cell)
  outputs_cell = outputs_cell*density_ratio[cell_numbers == cell]
  res = apply(outputs_cell,2,sum)
  res = array(res, dim = dim(outputs)[1:(length(dim(outputs))-1)])
  }
  return(res)
}




