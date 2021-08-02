
stan_model_baseline <- stan_model(file = 'PoissonModel.stan')
stan_model_proportionalCR <- stan_model(file = 'PoissonModel_CR.stan')
stan_model_ageFOI <- stan_model(file = 'PoissonModelAge.stan')


data=  list(CHIKV_obs =  data_csv$CHIKV_obs,
            ONNV_obs = data_csv$ONNV_obs,
            sex = data_csv$sex,
            location = data_csv$location,  
            m=4,
            N=length(data_csv$CHIKV_obs),
            ages=data_csv$age)

 
init = list(list(s1C=2,
                 s1N=3.9,
                 SCN=1.5,
                 sNC=2.7,
                 pCN=0.8,
                 pNC=0.2,
                 qC=0.23,
                 qN=0.43,
                 qC_Martinique=5,
                 qN_Martinique=0.01,
                 sexC = -0.6,
                 sexN=-0.1,
                 locC = -0.2,
                 locN =-0.3))


fit <- sampling(stan_model_baseline, data = data, chains = 1, cores=1,iter=10000 )

#fit <- sampling(stan_model_ageFOI, data = data, chains = 1, cores=1,iter=10000, init=init)
#fit <- sampling(stan_model_proportionalCR, data = data, chains = 1, cores=1,iter=10000, init=init)

 
