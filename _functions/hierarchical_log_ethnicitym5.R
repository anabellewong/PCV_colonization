hierarchical_model_log <-function(log_gmc_mean, log_gmc_prec){
  model_string<-"
model{
     
for(i in 1:N){

#Likelihood
y[i] ~ dbern(pi[i])

  log(qi[i]) = b0[st.index[i]] + bedouin[i]*b3[st.index[i]] + b1[st.index[i]]*(log_GMC[i] - cp1[st.index[i]])*step(log_GMC[i] - cp1[st.index[i]])

   pi[i] = qi[i]*(1-step(qi[i]-1))  + step(qi[i]-1)*0.999999
}

for(k in 1:N.serotypes){
  b0[k] ~ dnorm(mu_beta0,prec_beta0)
  b3[k] ~ dnorm(mu_beta3,prec_beta3)
  b1[k] ~ dnorm(mu_beta1,prec_beta1)
  cp1[k] ~ dnorm(mu_cp1,prec_cp1)
}

# mu_beta0~dnorm(0,1)
# mu_beta1~dnorm(0,1)
# mu_beta3~dnorm(0,1)
mu_beta0~dnorm(0, 1e-4)
mu_beta1~dnorm(0, 1e-4)
mu_beta3~dnorm(0, 1e-4)
mu_cp1 ~ dnorm(log_gmc_mean,log_gmc_prec)

# prec_beta0 ~ dgamma(4, 2)
# prec_beta1 ~ dgamma(4, 2)
# prec_beta3 ~ dgamma(4, 2)
prec_beta0 ~ dgamma(0.01, 0.01)
prec_beta1 ~ dgamma(0.01, 0.01)
prec_beta3 ~ dgamma(0.01, 0.01)
prec_cp1 ~ dgamma(0.01, 0.01)

}
"

##############################################################
#Model Fitting
##############################################################
inits1=list(".RNG.seed"=c(123), ".RNG.name"='base::Wichmann-Hill')
inits2=list(".RNG.seed"=c(456), ".RNG.name"='base::Wichmann-Hill')
inits3=list(".RNG.seed"=c(789), ".RNG.name"='base::Wichmann-Hill')

##############################################
#Model Organization
##############################################
model_spec<-textConnection(model_string)
model_jags<-jags.model(model_spec, 
                       inits=list(inits1,inits2, inits3),
                       data=list('y' = a1$Colonization,
                                 'N.serotypes'=max(a1$st.index),
                                 'N' = nrow(a1),
                                 'log_GMC'=a1$log_GMC,
                                 'st.index'=a1$st.index,
                                 'bedouin'=a1$bedouin,
                                 'log_gmc_mean'=log_gmc_mean,
                                 'log_gmc_prec'=log_gmc_prec
                                 ),
                       n.adapt=10000, 
                       n.chains=3)

params<-c('b0',
          'b1','cp1','b3')

##############################################
#Posterior Sampling
##############################################
posterior_samples<-coda.samples(model_jags, 
                                params, 
                                n.iter=10000)
posterior_samples.all<-do.call(rbind,posterior_samples)
#post1.summary<-summary(posterior_samples)
#post_means<-colMeans(posterior_samples.all)

post_means<-apply(posterior_samples.all, 2, median)
sample.labs<-names(post_means)
ci<-t(hdi(posterior_samples.all, credMass = 0.95))
ci<-matrix(sprintf("%.1f",round(ci,1)), ncol=2)
row.names(ci)<-sample.labs
post_means<-sprintf("%.1f",round(post_means,1))
names(post_means)<-sample.labs

res.list <- list('post_means'=post_means,'ci'=ci, 'posterior_samples.all'=posterior_samples.all)
return(res.list)

}