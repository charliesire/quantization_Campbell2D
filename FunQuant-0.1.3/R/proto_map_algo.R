
#' @title The prototype maps algorithm with multistart. Providing prototype outputs and their probability masses.

#' @param gamma Optional. If multistart = 1, gamma is a list of initial prototypes. Else, gamma is a list of list of prototypes, which length is equal to multistart.
#' @param nb_cells Required only if gamma is NULL. Indicates the number of prototypes of the quantization.
#' @param outputs The output samples that need to be quantized. If method = "percell", a list of output samples must be provided, of length equal to the number of Voronoï cells.
#' @param multistart Number of starts of the algorithm
#' @param method_IS The method of Importance Sampling : "unique" means there is a unique biased density involved, "percell" means there is one biased density (and then one biased sample) for each cell. Default is "unique".
#' @param density_ratio A vector indicating the weight fX/g of each output. Default is a vector of 1. If method = "percell", a list of density_ratio must be provided, of length equal to the number of Voronoï cells.
#' @param budget The maximum number of iterations of the algorithm. Default is 10^3.
#' @param distance_func A function computing a distance between two elements in the output spaces.
#' @param print_progress A boolean indicating whether to print the progress through the start numbers. Default is FALSE.
#' @param threshold A real positive number. When the distance between the new centroids and the previous ones is lower than this value, then we stop the algorithm.
#' @param trace A boolean. If TRUE, tracing information on the progress of the algorithm is produced. Default is FALSE.
#' @param bias A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoi cell. Default is 0 for all Voronoi cells.
#' @param index_sampling_error Required only if method_IS = "percell". Indicates which of the outputs samples must be used for the computation of the quantization error.
#' @param all_starts A boolean indicating whether the function should return the optimal prototypes obtained for each start.
#'
#' @return A list containing :
#' - gamma : the list of optimal prototypes
#' - probas : a vector indicating the probability mass of the prototypes
#' - cell_numbers : a vector indicating the cell number associated to each output
#' - iterations : an integer indicating the number of iterations performed
#' - record : a list containing all the centroids computed through the iterations of the best start
#' - all_errors : a vector indicating the quantization error of each start
#' - all_starts : a list indicating all the best prototypes obtained for each start. Provided only if all_start = TRUE.
#'
#' @export
#' @import abind
#' @examples
#' set.seed(20)
#' outputs = array(runif(9*20)*15, dim = c(3,3,20))
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' proto_map_algo(nb_cells = 3, multistart = 2,outputs = outputs, distance_func = distance_func)

proto_map_algo = function(gamma = NULL, nb_cells = NULL, outputs, multistart = 1, method_IS = "unique", density_ratio = rep(1, dim(outputs)[length(dim(outputs))]), budget = 10^3, threshold = 0, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},print_progress = FALSE, trace = FALSE, all_starts = FALSE,bias = NULL, index_sampling_error = NULL, seed = NULL, batch = FALSE){
  if(is.null(dim(outputs))){outputs = t(as.matrix(outputs))}
  if(!is.null(seed)){set.seed(seed)}
  best_error = +Inf
  best_list = list()
  all_errors = c()
  all_proto = list()
  for(it in 1:multistart){
    if(print_progress){print(paste("start number", it))}
    if(is.null(gamma)){
      are_distinct = FALSE
      while(are_distinct == FALSE){
        gamma_it = lapply(sample(x = 1:dim(outputs)[length(dim(outputs))], size = nb_cells), function(i){asub(outputs,dims = length(dim(outputs)), idx = i)})
        are_distinct = distinct_gamma(gamma_it)
      }
    }
    else{
      if(multistart == 1){gamma_it = gamma}
      else{gamma_it = gamma[[it]]}
    }
    gamma_it = oned_to_matrix(gamma_it)
    gamma_it = sort_gamma(gamma_it)
    if(is.null(bias)){bias = rep(0, length(gamma_it))}
    record = NULL
    if(trace){record = list(gamma_it)}
    for (i in 1:budget){
      if(print_progress){print(paste("Iteration number", i))}
      if(method_IS == "unique"){
        if(batch){cell_numbers = lapply(outputs,function(x){get_cell_numbers(x,gamma_it,distance_func = distance_func)})}
        else{cell_numbers = get_cell_numbers(outputs,gamma_it,distance_func = distance_func)}
        if(i==1){cccc <<- cell_numbers}
        if(length(table(unlist(cell_numbers)))<length(gamma_it)){stop("One of the initial Voronoi cells is empty")}
      }
      else if(method_IS == "percell"){
        cell_numbers = lapply(1:length(outputs), function(h){get_cell_numbers(outputs[[h]],gamma_it,distance_func = distance_func)})
        if(0 %in% sapply(1:length(cell_numbers), function(y){sum(cell_numbers[[y]] == y)}) ){stop("One of the initial Voronoi cells is empty")}
      }
      centro_probas = compute_centroids_and_proba(outputs = outputs, cell_numbers = cell_numbers,method_IS = method_IS, density_ratio = density_ratio, bias = bias, batch = batch)
      gamma_new = oned_to_matrix(centro_probas$centroids)
      probas = centro_probas$probas
      gamma_evol = sum(Vectorize(function(m){distance_func(gamma_it[[m]], gamma_new[[m]])})(1:length(gamma_it)))
      gamma_it = gamma_new
      if(trace){record[[i]] = gamma_it}
      if (gamma_evol <= threshold){break}
    }
    gamma_it = sort_gamma(gamma_it)
    error = quanti_error(outputs = outputs, gamma = gamma_it, density_ratio = density_ratio, batch = batch, distance_func = distance_func)
    all_errors = c(all_errors, error)
    if(all_starts){all_proto[[it]] = gamma_it}
    if(error < best_error){
      best_error = error
      best_list = list(gamma = gamma_it, probas = probas, cell_numbers = cell_numbers, iterations = i, record = record)
    }
  }
  return(c(best_list, list(all_errors = all_errors, all_starts = all_proto)))
}
