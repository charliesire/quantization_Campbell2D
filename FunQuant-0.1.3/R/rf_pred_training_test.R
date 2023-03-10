#' @title Predict the output class of a validation sample with random forest classifier for different hyperparameters values
#'
#' @param xtrain A dataframe of training inputs
#' @param xtest A dataframe of test inputs
#' @param ytrain A vector of training outputs for the random forests
#' @param list_search A list containing for each hyperparameters to be tested a list of the tested values.
#' @param ... other parameters of \code{\link{randomForest}} function from \code{randomForest}.
#'
#' @return A list containing a vector of predicted classes for each combination of tested hyperparameters.
#' @export
#' @importFrom randomForest randomForest
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
#' x = as.data.frame(sobol(50,2))*2-1
#' xtrain = x[1:40,]
#' xtest = x[41:50,]
#' ytrain = as.factor(Vectorize(function(i){sum(func2D(xtrain)[,,i])})(1:nrow(xtrain)) > 2)
#' ytest = as.factor(Vectorize(function(i){sum(func2D(xtest)[,,i])})(1:nrow(xtest)) > 2)
#' df_search = expand.grid(seq(0.1,1,0.3), c(1,5,9,13,17))
#' list_search = list("nodesize" = as.list(df_search[,2]), "classwt" =
#' lapply(1:nrow(df_search), function(i){c(df_search[i,1], 1-df_search[i,1])}))
#' list_rf_train_test = rf_pred_training_test(xtrain = xtrain, xtest = xtest,
#' ytrain = ytrain, list_search = list_search)

rf_pred_training_test = function(xtrain, xtest, ytrain, list_search,...){
  pred = list()
  for(i in 1:length(list_search[[1]])){
    liste_cv = list()
    for(var in names(list_search)){
      liste_cv = c(liste_cv, var = list_search[[var]][[i]])
    }
    liste_cv = c(liste_cv,list("x" = xtrain, "y" = ytrain, "xtest" = xtest,...))
    rf = do.call(randomForest, liste_cv)
    pred[[i]] = rf$test$predicted
  }
  return(pred)
}


