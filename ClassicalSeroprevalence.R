# Array with only the VNT from Martinique
xobs = cbind(data_csv$CHIKV_obs,data_csv$ONNV_obs)
N = nrow(xobs)

classical.positivity.chik <- function(X){
  return(X[,1]>=1 & X[,2]<=X[,1]-2)
}
classical.positivity.onnv <- function(S){
  return(X[,2]>=1 & X[,1]<=X[,2]-1)
}



#Samples were considered CHIKV PRNT positive 
#if the titer was >20 and the ONNV titer was < four-fold 
#lower than the CHIKV titer. Because there is a unique one-way antigenic 
#cross-reactivity between CHIKV and ONNV, a sample was designated ONNV positive 
#if its titer was >20 and two-fold or greater than the CHIKV titer.

X=xobs[which(data_csv$location<3),]
CHIKV.positive = which(classical.positivity.chik(X))
ONNV.positive =which(classical.positivity.onnv(X))

print('CHIKV positive Mali')
print(paste0('N = ', length(CHIKV.positive)))
print(paste0('proportion = ', length(CHIKV.positive)/dim(X)[1]))
 
print('ONNV positive Mali')
print(paste0('N = ', length(ONNV.positive)))
print(paste0('proportion = ', length(ONNV.positive)/dim(X)[1]))


# In Martinique
X=xobs[which(data_csv$location==3),]
CHIKV.positive = which(assess.positivity.chik(X))
ONNV.positive = which(assess.positivity.onnv(X))
print('CHIKV positive Martinique')
print(paste0('N = ', length(CHIKV.positive)))
print('ONNV positive Martinique')
print(paste0('N = ', length(ONNV.positive)))
