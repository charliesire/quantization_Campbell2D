% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_proba_training_test.R
\name{rf_probas_training_test}
\alias{rf_probas_training_test}
\title{Computing the relative errors when predicting the membership probabilities of a validation output sample with a classification step, for different values of hyperparameters}
\usage{
rf_probas_training_test(
  design_train,
  design_test,
  outputs_train,
  outputs_test,
  threshold_classification,
  threshold_fpca = NULL,
  list_search,
  density_ratio,
  gamma,
  distance_func = function(A1, A2) {
     return(sqrt(sum((A1 - A2)^2)))
 },
  return_pred = FALSE,
  only_positive = FALSE,
  seed = NULL,
  ncoeff,
  npc,
  formula = ~1,
  covtype = "matern5_2",
  wf = "d4",
  boundary = "periodic",
  J = 1,
  coef.trend = NULL,
  coef.cov = NULL,
  coef.var = NULL,
  nugget = NULL,
  noise.var = NULL,
  lower = NULL,
  upper = NULL,
  parinit = NULL,
  multistart = 1,
  kernel = NULL,
  control = NULL,
  type = "UK",
  bias = NULL,
  ...
)
}
\arguments{
\item{design_train}{a data frame representing the design of experiments of the training part.
The ith row contains the values of the d input variables corresponding
to the ith evaluation.}

\item{design_test}{a data frame representing the design of experiments of the validation part.
The ith row contains the values of the d input variables corresponding
to the ith evaluation.}

\item{outputs_train}{The training output samples on which the metamodel will be trained}

\item{outputs_test}{The validation output samples on which the metamodel performance will be evaluated}

\item{threshold_classification}{The threshold that creates the two classes of maps for the classification}

\item{threshold_fpca}{The threshold used for the training of the FPCA. Only the maps for which the sum of the pixel is above this threshold are used for the training. If NULL, this threshold takes the value of threshold_classification.}

\item{list_search}{A list containing for each hyperparameters to be tested a list of the tested values.}

\item{density_ratio}{density_ratio indicates the weight fX/g of each output}

\item{gamma}{A set of l prototypes defining the Voronoï cells}

\item{distance_func}{A function computing a distance between two elements in the output spaces.}

\item{return_pred}{A boolean indicating whether the predicted outputs should be returned or not}

\item{only_positive}{A boolean indicating whether the predicted outputs should only contained positive values or not. Default is FALSE.}

\item{seed}{An optional random seed}

\item{ncoeff}{The number of coefficients used for PCA}

\item{npc}{The number of principal components}

\item{formula}{an object of class "formula"
(or a list of "formula" which the length is equal to the number of modeled principal components)
specifying the linear trend of the kriging model (see \code{\link{lm}}) on each principal component.
 This formula should concern only the input variables (\code{design}), and not the output (\code{response}).
 The default is ~1, which defines a constant trend on each principal component.#'}

\item{covtype}{optional character string or vector of character strings
specifying the covariance structure to be used on each modeled principal component
(see \code{\link{km}} for possible inputs of \code{covtype}).
If a vector, the length should be equal to the number of modeled principal components.}

\item{wf}{name of the wavelet filter to use in the decomposition}

\item{boundary}{a character string which specifies how boundaries are treated. Only "periodic" is currently implemented (see \code{\link{dwt.2d}}).}

\item{J}{depth of the wavelet decomposition, must be a number less than or equal to log(min(M,N),2). Default is 1.}

\item{coef.trend, coef.cov, coef.var}{optional vectors or matrices containing
the values for the trend, covariance and variance parameters.
If matrices, the number of rows should be equal to the number of modeled principal components.
For details, see \code{\link{km}}).}

\item{nugget}{an optional variance value or vector standing for the homogeneous nugget effect.
If vector, the length should be equal to the number of modeled principal components.}

\item{noise.var}{an optional vector or matrix containing the noise variance
at each observation on each modeled principal component.}

\item{lower, upper}{optional vectors or matrices containing the bounds of the correlation parameters
of each principal component for optimization. For details, see \code{\link{km}}).}

\item{parinit}{an optional vector or matrix containing the initial values for the variables to be optimized over.
For details, see \code{\link{km}}).}

\item{multistart}{an optional integer indicating the number of initial points from which running the BFGS optimizer.
(see \code{\link{km}}).}

\item{kernel}{an optional function or list of functions containing a new covariance structure
for each principal component. At this stage, the parameters must be provided as well, and are not estimated.}

\item{control}{an optional list of control parameters for optimization. For details, see \code{\link{km}}).}

\item{type}{A character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK"). Default is "UK.}

\item{bias}{A vector indicating the bias that came out when computing the importance sampling estimators of the membership probabilities. Each element of the vector is associated to a Voronoï cell.}

\item{...}{other parameters of \code{\link{randomForest}} function from \code{randomForest}.}
}
\value{
A list containing several outputs :
- list_search the list containing for each hyperparameters to be tested a list of the tested values.
- probas_pred_df a dataframe indicating for each combination of hyperparameters values the obtained predicted membership probabilities
- relative_error_df a dataframe indicating for each combination of hyperparameters values the relative error when predicting the membership probabilities
- outputs_pred an array providing the predicted outputs if return_pred is TRUE. If return_pred is FALSE, then outputs_pred is NULL.#' @export
}
\description{
Computing the relative errors when predicting the membership probabilities of a validation output sample with a classification step, for different values of hyperparameters
}
\examples{
 set.seed(5)
 func2D <- function(X){
 Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
 n<-nrow(X)
 Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*exp(-((0.8*Zgrid$z1+0.2*
 Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**2))*(
 Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
 Ymaps<- array(unlist(Y),dim=c(20,20,n))
return(Ymaps)
}
library(randtoolbox)
design = as.data.frame(sobol(300,2))*2-1
outputs = func2D(design)
design_train = design[1:250,]
design_test = design[251:300,]
outputs_train = outputs[,,1:250]
outputs_test = outputs[,,251:300]
df_search = expand.grid(seq(0.1,1,0.3), c(1,5,9,13,17))
list_search = list("nodesize" = as.list(df_search[,2]), "classwt" = lapply(1:nrow(df_search),
function(i){c(df_search[i,1], 1-df_search[i,1])}))
density_ratio = rep(1,50)
gamma = lapply(c(2,3,51,7), function(i){outputs[,,i]})
distance_func = function(A1,A2){return(sqrt(sum((A1-A2)^2)))}
list_rf_prob_train_test = rf_probas_training_test(design_train = design_train,
design_test = design_test, outputs_train = outputs_train,
outputs_test = outputs_test, threshold_classification = 2, threshold_fpca = 0, list_search = list_search,
density_ratio = density_ratio, gamma = gamma, distance_func= distance_func,
ncoeff = 400, npc = 6, control = list(trace = FALSE))
}
