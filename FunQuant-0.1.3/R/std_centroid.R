#' @title Computation of the IS centroid standard deviation for different sets of prototypes.
#'
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma_list A list of gamma on which we want to evaluate the IS centroid standard deviation. Each gamma is a set of prototypes.
#' @param density_ratio density_ratio indicates the weight fX/g of each output.
#' @param distance_func  A function computing a distance between two elements in the output spaces.
#' @param cells The Voronoï cell numbers that we are investigating.
#' @param nv The size of the sample for which we want to estimate the IS centroid standard deviation.
#' @param cell_numbers An optional list providing for each set of prototypes the voronoi cell number of every output.
#' @return A list providing for each set of prototypes a list the IS centroid standard deviation for each voronoi cell
#' @export
#' @import abind
#' @importFrom stats var cov
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-
#' 10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 50))
#' outputs = func2D(design)
#' gamma_list = list(lapply(c(1,3,10,14,18), function(i){outputs[,,i]}))
#' density_ratio = rep(1, 50)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' list_std_centroid = std_centroid(outputs = outputs, gamma_list =
#' gamma_list, density_ratio = density_ratio, distance_func = distance_func
#' , cells = 1:length(gamma_list[[1]]), nv = 50)
std_centroid = function(outputs, gamma_list, density_ratio, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},cells, cell_numbers = NULL, nv){

  weighted_map = t(matrix(outputs, nrow = prod(dim(outputs)[-length(dim(outputs))]),ncol = dim(outputs)[length(dim(outputs))]))*density_ratio
  std_ratio_list = list()
  for(it in 1:length(gamma_list)){
    gamma = gamma_list[[it]]#for all Gamma
    std_ratio_list[[it]] = as.list(rep(0,length(cells)))
    if(is.null(cell_numbers)){cell_numbers_it = get_cell_numbers(outputs = outputs, gamma = gamma, distance_func = distance_func)}
    else{cell_numbers_it = cell_numbers[[it]]}
    for(j in cells){#for all voronoi cells
      map_loop = weighted_map #weighted map is the set of maps multiplied by the weights f_{x}/mu
      density_voronoi = density_ratio*(cell_numbers_it == j) #density_voronoi is the vector of the weights multiplied by the characteristic function of the voronoi cell
      for(i in 1:length(cell_numbers_it)){
        if(cell_numbers_it[i] != j){map_loop[i,] = rep(0,length(map_loop[i,]))} #this
      } #map_loop is now the set of maps multiplied by the weights f_{x}/mu multiplied by the characteristic function of the voronoi cell
      covariance = apply(map_loop,2,function(x){cov(x,density_voronoi)})/nv #this is the covariance between Ai and B for all i
      moy1 = apply(map_loop,2, mean) #This is the empirical expectation of Ai for all i
      moy2 = mean(density_voronoi) #This is the empirical expectation of B
      var1 = apply(map_loop,2,var)/nv #This is the empirical variance of Ai for all i
      var2 = var(density_voronoi)/nv #This is the empirical variance of B
      std_ratio_list[[it]][[j]] = sqrt(1/moy2^2*(var1+var2/moy2^2*moy1^2-2*covariance*moy1/moy2)) #This is the empirical standard error of the ratio

    }
  }
  return(std_ratio_list)
}
