
loglikehood<-function(pC, 
                      pN,
                      CHIKV,  
                      ONNV, 
                      m, 
                      s1C,
                      s1N,
                      pCN, 
                      pNC, 
                      sCN, 
                      sNC){
  
  V1=0 
  
  P00 =  (1-pC)*(1-pN)
  P10 =  pC*(1-pN)
  P01 =  (1-pC)*pN
  P11 =  pC*pN
  P0=0
  P1=0
  if(CHIKV==0){
    P0=1
  }
  if(ONNV==0){
    P1=1
  }
  A =  P00*P0*P1
  
  
  P2= dresponse(CHIKV, s1C)
  P5= dcrossreactivity(ONNV, CHIKV,  sCN)
  
  A =  A+ P10*P2*(P5*pCN+P1*(1-pCN))
  
  P3 = dcrossreactivity(CHIKV, ONNV, sNC)
  P4= dresponse(ONNV, s1N)
  
  A =  A+ P01*P4*(P3*pNC+P0*(1-pNC))
  
  B=(1-pCN)*(1-pNC)*P2*P4
  
  V1=0
  for(k in 0:ONNV){
    V1=V1+dresponse(k,s1N)*dcrossreactivity(ONNV-k, CHIKV, sCN) 
  }
  B=B + pCN*(1-pNC)*P2*V1
  
  
  V1=0
  for(k in 0:CHIKV){
    V1=V1+dresponse(k,s1C)*dcrossreactivity(CHIKV-k, ONNV,sNC) 
  }
  B=B +(1-pCN)*pNC*P4*V1
  
  V1=0
  for(kC in 0:CHIKV){
    for(kN in 0:ONNV){
      V1=V1+dresponse(kC, s1C)*dresponse(kN, s1N)*dcrossreactivity(CHIKV-kC, kN,  sNC)*dcrossreactivity(ONNV-kN, kC,  sCN)
    }
  }
  
  B=B + pCN*pNC*V1
  A=A+B*P11   
  
  return(A)
}

LogProb <- function(data, params){
  CHIKV_obs= data$CHIKV_obs
  ONNV_obs= data$ONNV_obs
  N= data$N
  sex = data$sex
  location=data$location
  m= data$m
  
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
  
  LP=rep(0,N)
  
  for(i in 1:N){
    A = 0
    CHIKV = CHIKV_obs[i]
    ONNV = ONNV_obs[i]
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
    
    A = loglikehood(pC,pN,
                    CHIKV, ONNV, m,
                    s1C, s1N,
                    pCN,pNC,sCN,sNC)
    
    if(CHIKV==m & ONNV==m){
      A=0
      for(i1 in m:10){
        for(i2 in m:10){
          A = A+ loglikehood(pC,pN, 
                             i1, i2,m,
                             s1C, s1N,
                             pCN,pNC,
                             sCN,sNC)              
        }
      }
    }
    if(CHIKV==m & ONNV<m){
      A=0
      for(i1 in m:10){
        A = A+ loglikehood(pC,pN, i1, ONNV,m,
                           s1C, s1N,
                           pCN,pNC,
                           sCN,sNC)
      }
    }
    if(CHIKV<m & ONNV==m){
      A=0
      for(i2 in m:10){
        A = A+ loglikehood(pC,pN, CHIKV, i2,m,
                           s1C, s1N,
                           pCN,pNC,
                           sCN,sNC)              
      }
    }
    
    if(A==0){
      A  = -10000
    }
    
    LP[i]= log(A)
  } 
  return(LP)
}


mean.params = list(qC=mean(Chains$qC),
                   qN=mean(Chains$qN),
                   qC_Martinique=mean(Chains$qC_Martinique),
                   qN_Martinique=mean(Chains$qN_Martinique),
                   sexC=mean(Chains$sexC),
                   sexN=mean(Chains$sexN),
                   locC=mean(Chains$locC),
                   locN=mean(Chains$locN),
                   pNC=mean(Chains$pNC),
                   pCN=mean(Chains$pCN),
                   sNC=mean(Chains$sNC),
                   sCN=mean(Chains$sCN),
                   s1C=mean(Chains$s1C),
                   s1N=mean(Chains$s1N))


LP1=rep(0,length(Chains$qC))
for(I in 1:length(Chains$qC)){
  print(I)
  params1 = list(qC=Chains$qC[I],
                 qN=Chains$qN[I],
                 qC_Martinique=Chains$qC_Martinique[I],
                 qN_Martinique=Chains$qN_Martinique[I],
                 sexC=Chains$sexC[I],
                 sexN=Chains$sexN[I],
                 locC=Chains$locC[I],
                 locN=Chains$locN[I],
                 pNC=Chains$pNC[I],
                 pCN=Chains$pCN[I],
                 sNC=Chains$sNC[I],
                 sCN=Chains$sCN[I],
                 s1C=Chains$s1C[I],
                 s1N=Chains$s1N[I])
  LP1[I] = sum(LogProb(data,params1))
  
}

LP.Mean = LogProb(data,mean.params)

L1= -2*mean(LP1)
L2=-2*sum(LP.Mean)

pD=L1-L2 # effective number of parameters
DIC = L2+2*pD

