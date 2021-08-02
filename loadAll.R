
########     ###    ########    ###   
##     ##   ## ##      ##      ## ##  
##     ##  ##   ##     ##     ##   ## 
##     ## ##     ##    ##    ##     ##
##     ## #########    ##    #########
##     ## ##     ##    ##    ##     ##
########  ##     ##    ##    ##     ##



data_csv=read.csv('data.csv', sep=';')
#CHIKV_obs, ONNV_obs # the VNT for CHIKV and ONNV (in a log scale with (20,40, 80, 160, 320) --> (0,1,2,3,4) )
#sex  # sex of the participants (males = 1 or females = 2)
#location # locations of the sampling site (here, 1 = South Mali, 2 = North Mali, 3 = Martinique). 

N=length(data_csv$CHIKV_obs)
m=4


CHIKVColor  = 'orange'
ONNVColor = 'seagreen4'

quan1 <-function(X){
  return(quantile(X,probs = 0.025))
}

quan2 <-function(X){
  return(quantile(X,probs = 0.975))
}

mean.and.ci<-function(X, digits=2){
  print(paste0(round(mean(X), digits=digits) ," (", round(quan1(X), digits=digits), " - ", round(quan2(X), digits=digits), ")"))
}

mean.and.var.ZTP <- function(Lambda){
  m = Lambda/(1-exp(-Lambda))
  v  =(Lambda+Lambda^2)/(1-exp(-Lambda))-(Lambda^2)/(1-exp(-Lambda))^2
  print('Mean')
  mean.and.ci(m)
  print('s.d.')
  mean.and.ci(sqrt(v))
}

# density function for the response of the infecting virus. Zero-truncated Poisson distribution
dresponse <- function(  x,   lambda){
  if(x==0)
    a = 0
  if(x>0){
    a=dpois(x, lambda)/(1-exp(-lambda)) 
  } 
  return(a)
}

# random generation for the distribution of the response of the infecting virus
rresponse <- function( prob){
  a=0
  while(a==0){
    a=rpois(n=1, prob)
  } 
  return(a)
}

# density function for the cross-reactive response. Zero-truncated Poisson distribution
dcrossreactivity <- function(  x,   size,   lambda){
#  a = dresponse(x,lambda*(size)+0.01)
  a = dresponse(x,lambda)
  return(a)
}

# random generation for the distribution of the cross-reactive response
rcrossreactivity <- function( lambda, size ){
  #a = rresponse(lambda*(size)+0.01)
  a = rresponse(lambda)
  return(a)
}

simulate_titers <- function(params,data){
  sex=data$sex
  location=data$location
  N = data$N
  xrealC = c()
  xrealN = c()
  rC = c()
  rN = c()
  xunobsC=c()
  xunobsN = c()
  QC=c()
  QN = c()
  MaxValue= data$m
  
  qC=params$qC
  qN=params$qN
  sexC=params$sexC
  sexN=params$sexN
  locC=params$locC
  locN=params$locN
  qC_Martinique=params$qC_Martinique
  qN_Martinique=params$qN_Martinique
  
  s1C=params$s1C
  s1N=params$s1N
  
  pCN=params$pCN
  pNC=params$pNC
  sCN=params$sCN
  sNC=params$sNC
  
  pC1 = 1-exp(-qC) 
  pC2 = 1-exp(-qC*exp(sexC)) 
  pC3 = 1-exp(-qC*exp(locC)) 
  pC4 = 1-exp(-qC*exp(sexC)*exp(locC)) 
  
  pN1 = 1-exp(-qN) 
  pN2 = 1-exp(-qN*exp(sexN)) 
  pN3 = 1-exp(-qN*exp(locN)) 
  pN4 = 1-exp(-qN*exp(sexN)*exp(locN))   
  
  pC5  = 1-exp(-qC_Martinique)   
  pN5  = 1-exp(-qN_Martinique)   
  
  for(i in 1:N){
    
    LOC = location[i]
    SEX = sex[i] 
    
    if(LOC==1 & SEX==1){
      pC=pC1
      pN=pN1
    }
    if(LOC==1 & SEX==2){
      pC=pC2
      pN=pN2
    }
    if(LOC==2 & SEX==1){
      pC=pC3
      pN=pN3
    }
    if(LOC==2 & SEX==2){
      pC=pC4
      pN=pN4
    }
    if(LOC==3){
      pC=pC5
      pN=pN5
    }
    
    rC  = runif(n=1)<pC
    rN  = runif(n=1)<pN
    
    QC[i] = rC
    QN[i] = rN
    
    if(rC==FALSE & rN  == FALSE){
      xrealC[i] =0
      xrealN[i] = 0
    }  
    if(rC ==TRUE & rN == FALSE){
      xrealC[i] =  rresponse(prob =  s1C) 
      xrealN[i] =  0
      if(runif(n=1)<pCN){ # cross reactivity
        xrealN[i] =  rcrossreactivity( lambda = sCN,  size=  xrealC[i])
      }
    } 
    
    if(rC ==FALSE & rN == TRUE){
      
      xrealN[i] =  rresponse(prob =  s1N) 
      xrealC[i] =  0
      if(runif(n=1)<pNC){ # cross reactivity
        xrealC[i] =  rcrossreactivity( lambda = sNC, size=  xrealN[i])
      }
    }  
    
    if(rC ==TRUE & rN == TRUE){
      
      xC = rresponse(prob =  s1C) 
      xN=  rresponse(prob =  s1N) 
      xunobsC[i] = xC
      xunobsN[i] = xN
      
      xrealC[i]  = xC  
      xrealN[i]  = xN 
      
      if(runif(n=1)<pNC){ # cross reactivity
        xrealC[i]  = xC +   rcrossreactivity(lambda = sNC, size=xN)
      }
      if(runif(n=1)<pCN){ # cross reactivity
        xrealN[i]  = xN +   rcrossreactivity(lambda = sCN, size=xC)
      }
    } 
  }
  unique.titer=seq(0,MaxValue)
  xobsC= unique.titer[apply(as.array(xrealC),1,function(x) which.min(abs(x-unique.titer))) ]
  xobsN= unique.titer[apply(as.array(xrealN),1,function(x) which.min(abs(x-unique.titer))) ]
  #hist(xobsC)
  return(cbind((xobsC),(xobsN),QC,QN))
}

