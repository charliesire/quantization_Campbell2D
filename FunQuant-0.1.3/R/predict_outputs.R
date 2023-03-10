#' @title Predict outputs at new inputs, based on a training database.
#'
#' @param metamodel_fitted Optional list containing the different metamodel steps fitted on the training data. It contains a classifier, a Fpca2d object, and a list of km objects. Can be obtained with the function fit_metamodel.
#' @param design_train a data frame representing the design of experiments of the training part.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param design_test a data frame representing the design of experiments on which we want to predict outputs.
#' The ith row contains the values of the d input variables corresponding
#' to the ith evaluation.
#' @param outputs_train The training output samples on which the metamodel will be trained
#' @param only_positive A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.
#' @param seed An optional random seed
#' @param ncoeff The number of coefficients used for PCA
#' @param npc The number of principal components
#' @param formula  an object of class "formula"
#' (or a list of "formula" which the length is equal to the number of modeled principal components)
#' specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
#'  This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
#'  The default is ~1, which defines a constant trend on each principal component.#'
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
#' @param classification A boolean indicating whether a classification step should be integrated in the metamodel or not. Default is FALSE.
#' @param control_classification The list of hyperparameters of the classification function. Required only if classfification is TRUE.
#' @param threshold_classification The threshold that creates the two classes of maps for the classification
#' @param threshold_fpca The threshold used for the training of the FPCA. Only the maps for which the sum of the pixel is above this threshold are used for the training. If NULL, this threshold takes the value of threshold_classification.
#' @param ... other parameters of \code{\link{km}} function from \code{DiceKriging}.
#' @return An array containing the predicted outputs.
#' @export
#' @import waveslim
#' @import foreach
#' @import DiceKriging
#' @examples
#'  set.seed(5)
#'  func2D <- function(X){
#'  Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#'  n<-nrow(X)
#'  Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*
#'  exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**
#'  2))*(Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
#'  Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' library(randtoolbox)
#' design = as.data.frame(sobol(300,2))*2-1
#' outputs = func2D(design)
#' design_train = design[1:250,]
#' design_test = design[251:300,]
#' outputs_train = outputs[,,1:250]
#' outputs_test = outputs[,,251:300]
#' outputs_pred = predict_outputs(design_train = design_train,
#' design_test = design_test, outputs_train = outputs_train,
#' ncoeff = 400, npc = 6, control = list(trace = FALSE), classification = TRUE,
#' control_classification = list(nodesize = 4), threshold_classification = 2)

predict_outputs = function(metamodel_fitted = NULL, design_train = NULL, design_test, outputs_train = NULL, only_positive = FALSE, seed = NULL, ncoeff = NULL,npc = NULL, formula = ~1, covtype="matern5_2",wf = "d4", boundary = "periodic",J=1,
                           coef.trend = NULL, coef.cov = NULL, coef.var = NULL,
                           nugget = NULL, noise.var=NULL, lower = NULL, upper = NULL,
                           parinit = NULL, multistart=1,
                           kernel=NULL,control = NULL,type = "UK", classification = FALSE, control_classification = NULL,threshold_classification = NULL,threshold_fpca = NULL,...){
  if(is.null(metamodel_fitted)){
    metamodel_fitted = fit_metamodel(design_train = design_train, outputs_train = outputs_train, seed = seed, ncoeff = ncoeff, npc = npc,
                                     formula = formula, covtype = covtype, wf = wf, boundary = boundary, J=J,  coef.trend = coef.trend, coef.cov = coef.cov, coef.var = coef.var,
                                     nugget = nugget, noise.var=noise.var, lower = lower, upper = upper, parinit = parinit, multistart=multistart,
                                     kernel=kernel,control = control,type = type, classification = classification, control_classification = control_classification,threshold_classification = threshold_classification,threshold_fpca = threshold_fpca,...)
  }
  rf = metamodel_fitted$classifier
  fp = metamodel_fitted$fp
  model = metamodel_fitted$model
  pred_fpca = TRUE
  if(classification){
    rf_pred = as.numeric(predict(rf, design_test)) - 1
    design_test_fpca = design_test[rf_pred == 1,]
    if(sum(rf_pred) == 0){pred_fpca = FALSE}
    }
  else{
    design_test_fpca = design_test
  }

  if(pred_fpca){
    pred =  matrix(sapply(1:length(model), function(g){predict(object = model[[g]], newdata = design_test_fpca, type = type, compute = FALSE, checkNames = FALSE)$mean}), ncol = length(model))
    outputs_pred_draft = inverse_Fpca2d(pred,fp)
  }
  outputs_pred = array(0,dim = c(dim(fp$EigFct)[-length(dim(fp$EigFct))], nrow(design_test)))
  if(classification == FALSE){outputs_pred = outputs_pred_draft}
  else if(pred_fpca){
    outputs_pred = array(0,dim = c(dim(fp$EigFct)[-length(dim(fp$EigFct))], nrow(design_test)))
    dimnames(outputs_pred) = lapply(dim(outputs_pred), function(i){1:i})
    dimnames(outputs_pred_draft) = c(lapply(dim(outputs_pred_draft)[-length(dim(outputs_pred_draft))], function(i){1:i}), list(which(rf_pred == 1)))
    afill(outputs_pred) = outputs_pred_draft
  }
  if(only_positive){outputs_pred = (outputs_pred > 0)*outputs_pred}
  return(outputs_pred)
}
