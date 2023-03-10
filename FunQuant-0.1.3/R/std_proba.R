#' @title Computation of the IS coefficients of variation of the membership probability for different set of prototypes.
#'
#' @param outputs The output samples that need to be quantized. Useful only if cell_numbers == NULL.
#' @param gamma_list A list of gamma on which we want to evaluate the IS coefficient of variation of the membership probability. Each gamma is a set of prototypes
#' @param density_ratio density_ratio indicates the weight fX/g of each output
#' @param distance_func  A function computing a distance between two elements in the output spaces.
#' @param cells The Voronoï cell numbers that we are investigating
#' @param nv The size of the sample for which we want to estimate the IS coefficient of variation of the membership probability
#' @param cell_numbers An optional list providing for each set of prototypes the voronoi cell number of every output.
#'
#' @return A list of IS coefficients of variation of the membership probability obtained for each set of prototypes.
#' @export
#' @import abind
#' @importFrom  stats var
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2
#' -10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 50))
#' outputs = func2D(design)
#' gamma_list = list(lapply(c(1,3,10,14,18), function(i){outputs[,,i]}))
#' density_ratio = rep(1, 50)
#' distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
#' list_std_proba = std_proba(outputs = outputs, gamma_list = gamma_list,
#'  density_ratio = density_ratio, distance_func = distance_func,
#'  cells = 1:length(gamma_list[[1]]), nv = 50)
std_proba = function(outputs, gamma_list, density_ratio, distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))},cells, cell_numbers = NULL, nv){

  std_list = list()
  for(it in 1:length(gamma_list)){ #for all Gamma in (Gamma^r)

    if(is.null(cell_numbers)){cell_numbers_it = get_cell_numbers(outputs = outputs, gamma = gamma_list[[it]], distance_func = distance_func)}
    else{cell_numbers_it = cell_numbers[[it]]}
    df_for_std = Vectorize(function(i){density_ratio*(cell_numbers_it == i)})(cells)
    std_list[[it]] = apply(df_for_std, 2, function(x){sqrt(var(x))})/apply(df_for_std,2,mean)/sqrt(nv) #for all voronoi cells of Gamma, we compute the relative standard error
  }
  return(std_list)
}

