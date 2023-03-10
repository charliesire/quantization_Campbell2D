## Distances

load("NewFitting_Charlie_v090821.RData")

library(sp)
#library(raster)
library(ggplot2)
library(gtools)
library(dplyr)
library(POT)
library(pracma)
library(DiceDesign)
library(randtoolbox)
library(abind)
library(miceadds)
library(waveslim)
library(foreach)
library(DiceKriging)
library(orthogonalsplinebasis)

## Distance function


# distance_func computes the distance between two maps
distance_func = function(A1,A2){
  return(sqrt(sum((A1-A2)^2)))

}

#distance_gamma computes the distance between a map and a set of maps a returns the distance and the closest map of the set

distance_gamma = function(x, gamma){
  card = length(gamma)
  distance = c()
  for (i in 1:card){
    distance = c(distance, distance_func(x, gamma[[i]]))
  }
  return(list(cellule = which.min(distance), dist = min(distance)))

}
# list_loi_exp is a list containing the empirical observation of the variables related to the offshore conditions

list_loi_exp = list()
list_loi_exp[["T"]] = T #tide
list_loi_exp[["S"]] = S #surge
list_loi_exp[["phi"]] = phi #phi
list_loi_exp[["t1"]] = t1 #t1
list_loi_exp[["t2"]] = t2 #t2




proba_for_gpd = mean(list_loi_exp[["S"]] <= 0.55)



## The computation of the density function is directly related to the empirical density.
#For the surge, we add an information : the surge above 0.55 follows a gpd density. The density of the surge is also truncated between 0 and 2.5

density_S_base = function(x){density(list_loi_exp[[2]],  from = x, to = x, n = 1)$y[1]*(x <= 0.55) + as.numeric(dgpd(x,loc = 0.55, scale = fit$fitted.values[1], shape = fit$fitted.values[2]))*(x > 0.55)*(1-proba_for_gpd)}

offset = 0.65 # We will also add an offset representing an increase of the Surge peak. Then the surge is in the interval [offset, 2.5]

cste_trunc_S = stats::integrate(Vectorize(density_S_base), lower = -10, upper = 0)$value + (1 - proba_for_gpd)*stats::integrate(Vectorize(density_S_base), lower = 2.5 - offset, upper = 10)$value


density_i = function(x, i){
  if (i %in% c(1,3,4,5)){
    res = density(list_loi_exp[[i]], from = x, to = x, n = 1)$y[1]}
  else if (i == 2){
    if((x>=offset) & (x <= 2.5)){
    res = density_S_base(x-offset)/(1-cste_trunc_S)
    }
    else{res = 0}
  }
  return(res)
}


## bornes sup and bornes inf are the support of the each of the marginal density functions. Our biased density will be based on uniform density on these supports

bornes_inf = c()
bornes_sup = c()

for(i in 1:5){
  density_df = data.frame(seq(-15,15,l=5000),d = Vectorize(function(x){density_i(x,i)})(seq(-15,15,l=5000)))
  indice_inf = which.min(density_df$d == 0) - 1
  indice_sup = which.min(density_df[(indice_inf+1):nrow(density_df),"d"] != 0)
  bornes_inf = c(bornes_inf, density_df[indice_inf, 1])
  bornes_sup = c(bornes_sup, density_df[indice_sup + indice_inf,1])
}


df_bornes = data.frame(bornes_inf = bornes_inf[1], bornes_sup = bornes_sup[1])
for(i in 2:5){df_bornes = rbind(df_bornes, c(bornes_inf[i], bornes_sup[i]))}


density_xu = function(x){
  res = 1
  for (i in c(1,3,4,5)){
    res = res * density_i(x[i],i)
  }
  res = res *  density_i(x[2],2)
  return(res)
}


## Regarding the variables related to the breach, we consider that the probability of a breach depends on the signals of the tide and the surge.
## If the maximum of the sum of these signals is above 70% of the average size of the embankment, then the probability of breach is 1/2. Else it is 10^-4

#sum_water_signal sums the signals of tide and surge and takes the maximum

sum_water_signal = function(x){
  signal_T = x[1]*cos(2*pi/2.5*seq(-12,12,l=1001))
  function_S = function(t){
    t_pic = x[3]
    if(t > t_pic){return(max(0,x[2]*(1-(t-t_pic)/x[5])))}
    else{return(return(max(0,x[2]*(1+(t_pic-t)/x[4]))))}
  }
  signal_S = Vectorize(function_S)(seq(-12,12,l=1001))
  return(max(signal_S+signal_T))
}

fX = function(x){
  qte_eau = sum_water_signal(x)
  if(qte_eau <= 0.7*4.76){
    p_rup = 10^-4
    if(x[7]==0){
      return(density_xu(x[1:5])*(1-p_rup))
    }
    else{
      return(density_xu(x[1:5])*p_rup)
    }
  }
  if (qte_eau > 0.7*4.76){

    p_rup = 1/2

    if(x[7]==0){
      return(density_xu(x[1:5])*(1-p_rup))
    }
    else{
      return(density_xu(x[1:5])*p_rup*1/10)
    }
  }
}

q_x = 1
for (i in 1:5){q_x = q_x/(bornes_sup[i]-bornes_inf[i])}

g = function(x){
  if(x[7] == 0){return(q_x*5/13)}
  else{return(q_x*8/13*1/10)}
}


## Transform X converts X to the interval from [0,1]^7 to [-1,5]^7

transform_X = function(X){
  for (i in 1:(ncol(X)-1)){
    X[,i] = -1 + 6*X[,i]
  }
  return(X)

}

transform_X_from_phy = function(X){
  for(i in 1:5){
    X[,i] = (X[,i] - df_bornes[i,1])/(df_bornes[i,2]-df_bornes[i,1])
  }
  X[,6] = X[,6]/10
  return(transform_X(X))
}

## Plot the campbell maps

plot_map = function(map, max = NULL, min = NULL){
  gamma_toplot = expand.grid(z1,z2)
  gamma_toplot$f = as.numeric(map)
  if(is.null(max))
  {p = ggplot(gamma_toplot) + geom_raster(aes(x = Var1, y = Var2, fill = f)) + scale_fill_continuous(type = "viridis",direction = -1, name = "h") + theme_bw()}
  else{p = ggplot(gamma_toplot) + geom_raster(aes(x = Var1, y = Var2, fill = f)) + scale_fill_continuous(type = "viridis", direction = -1, limits = c(min, max), name = "h") + theme_bw()  + theme(legend.text = element_text(size=13),legend.title = element_text(size=13))}
  return(p)
}


