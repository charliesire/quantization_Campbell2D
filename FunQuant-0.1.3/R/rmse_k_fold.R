#' @title Computation of the Root Mean Square Error at each pixel of the outputs by kfold Cross Validation, for different values of ncoeff and npc
#'
#' @param outputs The output samples on which the metamodel will be trained and tested by kfold cross validation
#' @param nb_folds Number of folds
#' @param ncoeff_vec A vector providing the different values of ncoeff to be tested. ncoeff fixes the number of coefficients used for PCA.
#' @param npc_vec A vector providing the different numbers of principal components to be tested.
#' @param return_pred A boolean indicating whether the predicted outputs should be returned or not
#' @param formula  an object of class "formula"
#' (or a list of "formula" which the length is equal to the number of modeled principal components)
#' specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
#'  This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
#'  The default is ~1, which defines a constant trend on each principal component.
#' @param design a data frame representing the design of experiments.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param covtype optional character string or vector of character strings
#' specifying the covariance structure to be used on each modeled principal component
#' (see \code{\link{km}} for possible inputs of \code{covtype}).
#' If a vector, the length should be equal to the number of modeled principal components.
#' @param wf name of the wavelet filter to use in the decomposition
#' @param boundary a character string which specifies how boundaries are treated. Only "periodic" is currently implemented (see \code{\link{dwt.2d}}).
#' @param J depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.
#' @param coef.trend,coef.cov,coef.var optional vectors or matrices containing
#' the values for the trend, covariance and variance parameters.
#' If matrices, the number of rows should be equal to the number of modeled principal components.
#' For details, see \code{\link{km}}).
#' @param nugget an optional variance value or vector standing for the homogeneous nugget effect.
#' If vector, the length should be equal to the number of modeled principal components.
#' @param noise.var an optional vector or matrix containing the noise variance
#' at each observation on each modeled principal component.
#' @param lower,upper optional vectors or matrices containing the bounds of the correlation parameters
#' of each principal component for optimization. For details, see \code{\link{km}}).
#' @param parinit an optional vector or matrix containing the initial values for the variables to be optimized over.
#' For details, see \code{\link{km}}).
#' @param multistart an optional integer indicating the number of initial points from which running the BFGS optimizer.
#'  (see \code{\link{km}}).
#' @param kernel an optional function or list of functions containing a new covariance structure
#' for each principal component. At this stage, the parameters must be provided as well, and are not estimated.
#' @param control an optional list of control parameters for optimization. For details, see \code{\link{km}}).
#' @param type A character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK"). Default is "UK.
#' @param seed An optional random seed
#' @param ... other parameters of \code{\link{km}} function from \code{DiceKriging}.
#'
#'
#' @return A list containing several outputs :
#' - grid_cv is the grid made with the pairs (npc, ncoeff) that are tested
#' - outputs_rmse is a list of objects that have the same dimension as an output, obtained for each pair (npc, ncoeff). Each element (called pixel here) of the objects is the RMSE computed between the predicted values of the pixel and the true value of the pixel.
#' - outputs_pred is an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.
#' @export
#' @import waveslim
#' @import foreach
#' @import DiceKriging
#' @import abind
#' @importFrom randomForest randomForest
#' @importFrom dismo kfold
#' @examples
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){X[i,]*exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2
#' -10*X[i,])**2)/(60*X[i,]**2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,]*4)})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' design = data.frame(X = seq(-1,1,l= 20))
#' outputs = func2D(design)

#' list_rmse_k_fold = rmse_k_fold(outputs = outputs, nb_folds = 5,
#' ncoeff_vec = c(50,100,200,400), npc_vec = 2:4, design = design,
#' control = list(trace = FALSE))

rmse_k_fold = function(outputs, nb_folds, ncoeff_vec,npc_vec, return_pred = FALSE,formula = ~1,design, covtype="matern5_2", wf = "d4", boundary = "periodic",J=1,
                         coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                         nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                         parinit = NULL, multistart=1,
                         kernel=NULL,control = NULL,type = "UK",seed = NULL,...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  length_dim = length(dim(outputs))

  grid_cv = expand.grid(ncoeff = ncoeff_vec, npc = npc_vec)
  folds = kfold(dim(outputs)[length(dim(outputs))], nb_folds)
  outputs_rmse = list()
  outputs_pred_list = list()
  outputs_pred = list()
  relative_error_df = data.frame()
  for(k in 1:nb_folds){
    outputs_pred_list[[k]] = rmse_training_test(outputs_train = asub(outputs,dims = length(dim(outputs)), idx = which(folds != k)),outputs_test = asub(outputs,dims = length(dim(outputs)), idx = which(folds == k)), ncoeff_vec = ncoeff_vec,npc_vec = npc_vec, return_pred = TRUE,formula = formula,design_train = as.data.frame(design[folds !=k,]), design_test = as.data.frame(design[folds == k,]), covtype=covtype, wf = wf, boundary = boundary,J=J,
                                                coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                                nugget = nugget, noise.var=noise.var, lower = lower, upper = upper,
                                                parinit = parinit, multistart=multistart,
                                                kernel=kernel,control = control,type = type)$outputs_pred
  }
  for(i in 1:nrow(grid_cv)){
    outputs_pred[[i]] = array(0, dim = dim(outputs))
    dimnames(outputs_pred[[i]]) = lapply(dim(outputs_pred[[i]]), function(i){1:i})
    for(k in 1:nb_folds){
      dimnames(outputs_pred_list[[k]][[i]]) = c(lapply(dim(outputs_pred_list[[k]][[i]])[-length(dim(outputs_pred_list[[k]][[i]]))], function(i){1:i}), list(which(folds == k)))
      afill(outputs_pred[[i]]) = outputs_pred_list[[k]][[i]]
    }
    err = (outputs_pred[[i]] - outputs)^2
    outputs_rmse[[i]] = sqrt(apply(err, 1:(length(dim(err))-1), mean))
  }
  if(return_pred == FALSE){outputs_pred = NULL}
  return(list(grid_cv = grid_cv, outputs_rmse = outputs_rmse, outputs_pred = outputs_pred_list))
}


