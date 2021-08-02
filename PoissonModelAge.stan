functions{
  
  real gamma_discrete_cdf(int y1, real alpha, real beta, int maxv){
    real a;
    a=  (gamma_cdf(y1, alpha, beta)-gamma_cdf(y1-1,alpha,beta))/(gamma_cdf(maxv, alpha, beta));
    
    return(a);
  }
  
  real dresponse(int x, real lambda){
    real a;
    
    if(x==0)
    a = 0;
    if(x>0){
      a=exp(poisson_lpmf(x| lambda))/(1-exp(-lambda)) ;
    } 
    return(a);
  }
  
  real dcrossreactivity(int x, int size, real lambda){
    real a;
    a = dresponse(x,lambda);
    return(a);
  }
  
  real PNIPI(real pC,real pN,
  int CHIKV, int  ONNV, int m, 
  real s1C, real s1N,
  real pCN, real pNC, real sCN, real sNC){
    
    real A; 
    real B;
    real P00;
    real P10;
    real P01;
    real P11;
    real P0;
    real P1;
    real P2;
    real P3;
    real P4;
    real P5;
    real V1;
    
    V1=0; 
    
    
    P00 =  (1-pC)*(1-pN);
    P10 =  pC*(1-pN);
    P01 =  (1-pC)*pN;
    P11 =  pC*pN;
    P0=0;
    P1=0;
    if(CHIKV==0){
      P0=1;
    }
    if(ONNV==0){
      P1=1;
    }
    A =  P00*P0*P1;
    
    
    P2= dresponse(CHIKV, s1C);
    P5= dcrossreactivity(ONNV, CHIKV,  sCN);
    
    A =  A+ P10*P2*(P5*pCN+P1*(1-pCN));
    
    P3 = dcrossreactivity(CHIKV, ONNV, sNC);
    P4= dresponse(ONNV, s1N);
    
    A =  A+ P01*P4*(P3*pNC+P0*(1-pNC));
    
    B=(1-pCN)*(1-pNC)*P2*P4;
    
    V1=0;
    for(k in 0:ONNV){
      V1=V1+dresponse(k,s1N)*dcrossreactivity(ONNV-k, CHIKV, sCN) ;
    }
    B=B + pCN*(1-pNC)*P2*V1;
    
    V1=0;
    for(k in 0:CHIKV){
      V1=V1+dresponse(k,s1C)*dcrossreactivity(CHIKV-k, ONNV,sNC) ;
    }
    B=B +(1-pCN)*pNC*P4*V1;
    
    V1=0;
    for(kC in 0:CHIKV){
      for(kN in 0:ONNV){
        V1=V1+dresponse(kC, s1C)*dresponse(kN, s1N)*dcrossreactivity(CHIKV-kC, kN,  sNC)*dcrossreactivity(ONNV-kN, kC,  sCN);
      }
    }
    
    B=B + pCN*pNC*V1;
    A=A+B*P11;   
    
    return(A);
  }
  
  
  
}

data {
  
  int N;  // number of data
  int m; // maximal observable titer
  int location[N];
  int sex[N];
  int CHIKV_obs[N];
  int ONNV_obs[N];
  int ages[N];
  
}

parameters {
  // real<lower=0> alpha; 
  real<lower=0> s1C;
  real<lower=0> s1N;
  real<lower=0> sCN;
  real<lower=0> sNC;
  real<lower=0> pCN;
  real<lower=0> pNC;
  real<lower=0> qC;
  real<lower=0> qN;
  real<lower=0> qC_Martinique;
  real<lower=0> qN_Martinique;
  real sexC;
  real sexN;
  real locC;
  real locN;
  
}

transformed parameters{
  
}

model {
  
  real L; 
  int CHIKV;
  int ONNV;
  int LOC;
  int SEX; 
  real A;
  real pC;
  real pN; 
  int age;
  s1C ~ uniform(0,10);
  s1N ~ uniform(0,10);
  sNC ~ uniform(0,10);
  sCN ~ uniform(0,10);
  qC ~ uniform(0,10);
  qN ~ uniform(0,10);
  qC_Martinique ~ uniform(0,10);
  qN_Martinique ~ uniform(0,10);
  sexC ~ normal(0,3);
  sexN ~ normal(0,3);
  locC ~ normal(0,3);
  locN ~ normal(0,3);  
  pNC ~ uniform(0,1); 
  pCN ~ uniform(0,1); 
  
  
  
  
  for(i in 1:N){
    A = 0;
    CHIKV = CHIKV_obs[i];
    ONNV = ONNV_obs[i];
    LOC = location[i];
    SEX = sex[i]; 
    age = ages[i];
    
    
    if(LOC==1 && SEX==1){
      pC= 1-exp(-age*qC);
      pN= 1-exp(-age*qN);
    }
    if(LOC==1 && SEX==2){
      pC=1-exp(-age*qC*exp(sexC));
      pN=1-exp(-age*qN*exp(sexN));
    }
    if(LOC==2 && SEX==1){
      pC=1-exp(-age*qC*exp(locC));
      pN=1-exp(-age*qN*exp(locN));
    }
    if(LOC==2 && SEX==2){
      pC=1-exp(-age*qC*exp(sexC)*exp(locC));
      pN= 1-exp(-age*qN*exp(sexN)*exp(locN)) ; 
    }
    if(LOC==3){
      pC= 1-exp(-qC_Martinique); # no age available for Martinique
      pN= 1-exp(-qN_Martinique); 
    }
    
    A = PNIPI(pC,pN,
    CHIKV, ONNV, m,
    s1C, s1N,
    pCN,pNC,sCN,sNC);
    
    if(CHIKV==m && ONNV==m){
      A=0;
      for(i1 in m:10){
        for(i2 in m:10){
          A = A+ PNIPI(pC,pN, 
          i1, i2,m,
          s1C, s1N,
          pCN,pNC,
          sCN,sNC);              
        }
      }
    }
    if(CHIKV==m && ONNV<m){
      A=0;
      for(i1 in m:10){
        A = A+ PNIPI(pC,pN, i1, ONNV,m,
        s1C, s1N,
        pCN,pNC,
        sCN,sNC);
      }
    }
    if(CHIKV<m && ONNV==m){
      A=0;
      for(i2 in m:10){
        A = A+ PNIPI(pC,pN, CHIKV, i2,m,
        s1C, s1N,
        pCN,pNC,
        sCN,sNC);              
      }
    }
    
    if(A==0){
      A  = -10000;
    }
    
    target+= log(A);
  } 
}
