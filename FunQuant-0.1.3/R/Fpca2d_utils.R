#' @title Get approximation with scores of Fpca2d
#'
#' @param scores a matrix where columns correpond to scores of principal components.
#' @param fpca an object of class \code{Fpca2d}.
#'
#' @return approximation of scores.
#'
#'
#' @keywords internal
inverse_Fpca2d <-function(scores,fpca){

  # number of principal components
  nPC <- ifelse(is.matrix(scores), # condition
                ncol(scores), # TRUE
                length(scores)) # FALSE

  # number of scores
  n <- ifelse(is.matrix(scores), # condition
              nrow(scores), # TRUE
              1) # FALSE

  # coefficients of wavelet or B-splines basis
  coeff <- attr(fpca,"coeff")
  d<-dim(coeff)[1:2]
  K <- d[1]*d[2] # size of the functional basis

  # number of coefficients on PCA
  ncoeff <-attr(fpca,"ncoeff")

  # poe
  poe <-attr(fpca,"mean_poe")
  # coefficients on PCA
  idx_order <- order(poe, decreasing=TRUE)
  idx_pca <-idx_order[1:ncoeff]


  mu <-as.vector(apply(coeff, c(1,2), mean)) # mean coefficients

  #"---------------------------------"
  # estimation of coefficients on PCA
  #"---------------------------------"
  # rotation matrix
  rot <- attr(fpca,"pca")$rotation
  coeff_pca <- rot%*%t(scores)

  # estimation of all coefficients
  if(isTRUE(fpca$center)){
    # add mean values
    res_coeff <- sapply(1:n, function(i){
      ci <- rep(0,K)
      ci[idx_pca] <- coeff_pca[,i]
      return(ci+mu)
    })# end res_coeff
  }else{
    res_coeff <- sapply(1:n, function(i){
      ci <- rep(0,K)
      # coefficients on PCA
      ci[idx_pca] <- coeff_pca[,i]
      # NOT on PCA
      ci[-idx_pca] <- mu[-idx_pca]
      return(ci)
    })# end res_coeff
  }# end ifelse center

  # coefficients
  if(attr(fpca,"method")=="Bsplines"){
    res_coeff_mat <- res_coeff # keep matrix for B-splines
    res_coeff <- array(res_coeff_mat,dim=c(d[1],d[2],n))
    attr(res_coeff,"matrix")<-res_coeff_mat
    rm(res_coeff_mat) # memory deallocation
  }else{
    res_coeff <- array(res_coeff,dim=c(d[1],d[2],n))
  }# end ifelse



  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  if wavelet basis
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(attr(fpca,"method")=="wavelet"){

    attr(res_coeff,"type")<-attr(fpca,"type")
    attr(res_coeff,"J")<-attr(fpca,"J")  # depth of wavelet decomposition
    attr(res_coeff,"wf")<-attr(fpca,"wf") #wavelet type
    attr(res_coeff,"boundary")<-attr(fpca,"boundary")
    attr(res_coeff,"dim")<-c(d,n)

    class(res_coeff)<-"Wavelet2D"

    # approximation
    return(Inverse2D(res_coeff))

  }# end if wavelet

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  if B-splines basis
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(attr(fpca,"method")=="Bsplines"){
    OPhi<-attr(fpca,"SplinesBasis")
    return(Inverse2D(OPhi,res_coeff))
  }# end if Bsplines

} # end inverse_Fpca2d
