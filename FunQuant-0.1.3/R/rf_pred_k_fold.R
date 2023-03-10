#' @title Predict the output class with random forest classifier by kfold cross validation for different hyperparameters values
#'
#' @param x A dataframe of inputs
#' @param y A vector of categorical output
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param nb_folds Number of folds
#' @param seed An optional random seed
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.
#'
#' @return A list containing a vector of predicted classes for each combination of tested hyperparameters.
#' @export
#' @importFrom randomForest randomForest
#' @importFrom dismo kfold
#' @examples
#' set.seed(5)
#' func2D <- function(X){
#' Zgrid <- expand.grid(z1 = seq(-5,5,l=20),z2 = seq(-5,5,l=20))
#' n<-nrow(X)
#' Y <- lapply(1:n, function(i){(X[i,2] > 0)*X[i,2]*X[i,1]*
#' exp(-((0.8*Zgrid$z1+0.2*Zgrid$z2-10*X[i,1])**2)/(60*X[i,1]**2))*
#' (Zgrid$z1-Zgrid$z2)*cos(X[i,1]*4)^2*sin(X[i,2]*4)^2})
#' Ymaps<- array(unlist(Y),dim=c(20,20,n))
#' return(Ymaps)
#' }
#' library(randtoolbox)
#' x = as.data.frame(sobol(30,2))*2-1
#' outputs = func2D(x)
#' y = as.factor(Vectorize(function(i){sum(outputs[,,i])})(1:dim(outputs)[3]) > 5)
#' list_search = list("nodesize" = as.list(c(1,3,5,7,9,11)))
#' rf_pred_k_fold(x = x,y = y, list_search = list_search, nb_folds = 4)

rf_pred_k_fold = function(x, y, list_search, nb_folds, seed = NULL,...){
  if(is.null(seed)==FALSE){set.seed(seed)}
  folds = kfold(length(y), nb_folds)
  pred = list()
  for(i in 1:length(list_search[[1]])){
    liste_cv = list()
    pred[[i]] = rep(0, length(y))
    for(var in names(list_search)){
      liste_cv = c(liste_cv, var = list_search[[var]][[i]])
    }
      for(k in 1:nb_folds){
        liste_cv_fold = c(liste_cv, list("x" = x[folds!=k, ], "y" = y[folds!=k], "xtest" = x[folds==k,], "ytest" = y[folds == k],...))
        rf = do.call(randomForest, liste_cv_fold)
        pred[[i]][which(folds == k)] = rf$test$predicted
      }
  }
  return(pred)
}
