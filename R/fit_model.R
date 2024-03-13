library(ggplot2)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)

# set working directory to source file location

# read and process data
# note: this data is not publically available as a .csv file,
# but is printed in the accompanying manuscript.
# the necessary columns are 
# 1. the grouping variable (if you have multiple groups)
# 2. the stress
# 3. the lifetime
# note that for this application, the lifetime is capped at 1e7,
# and there were 3 groups of data.
df = read.csv("data/HT_data.csv")
df = df[,c('HT','MaxEngStress6_35mm','N')]
df$N_recip = -10^5/(as.numeric(df$N))
colnames(df) = c("Type","S","N","N_recip")

## 1 is an event, 0 is censored/runout
df$fracture = as.integer(df$N < 10^7)

# plotting raw data
ggplot(df,aes(y=S,x=N_recip,col=Type,shape=as.logical(fracture))) + geom_point() +
  facet_wrap(.~Type) +
  ylab("S") + 
  xlab("-1/N")
  

types = unique(df$Type)

for(ii in 1:length(types)) {
  
  # fit model for each group of data
  type_index = ii
  which_type = types[type_index]
  print(which_type)
  
  data = df[df$Type == which_type,]
  
  data = arrange(data,S)
  
  mod = lm(N_recip ~ S,data=data) # get reasonable prior centers
  smod = summary(mod)
  
  newdata = data.frame(S=seq(min(data$S),max(data$S),length.out=200),
                       N_recip=NA)
  
  # diagnostic plot to make sure fit looks reasonable
  ggplot(data,aes(x=N_recip,y=S)) + geom_point() +
    geom_vline(xintercept = 0,linetype='dashed') +
    geom_line(data=data.frame(N_recip=predict(mod,newdata=newdata),y=newdata$S),aes(x=N_recip,y=y),inherit.aes=F) +
    geom_line(data=data.frame(N_recip=predict(mod,newdata=newdata,interval='predict')[,2],y=newdata$S),aes(x=N_recip,y=y),inherit.aes=F,linetype='dashed') +
    geom_line(data=data.frame(N_recip=predict(mod,newdata=newdata,interval='predict')[,3],y=newdata$S),aes(x=N_recip,y=y),inherit.aes=F,linetype='dashed') +
    xlab("-1/N")
  
  untrans = function(x) -log10(-(x/100000))
  
  # plot on the original scale, note that the intervals are 
  # a bit conservative based on experience
  ggplot(data,aes(x=untrans(N_recip),y=S)) + geom_point() +
    geom_line(data=data.frame(N_recip=untrans(predict(mod,newdata=newdata)),y=newdata$S),aes(x=N_recip,y=y),inherit.aes=F) +
    geom_line(data=data.frame(N_recip=untrans(predict(mod,newdata=newdata,interval='predict')[,2]),y=newdata$S),aes(x=N_recip,y=y),inherit.aes=F,linetype='dashed') +
    xlim(c(NA,7))
  
  # stan inputs
  x_pred = seq(min(data$S),max(data$S),length.out=200)
  
  stan_data = list(
    N=nrow(data),
    N_cens=sum(data$fracture==0),
    N_obs=sum(data$fracture==1),
    N_pred=length(x_pred),
    x_obs=data$S[data$fracture==1],
    x_cens=data$S[data$fracture==0],
    x_pred=x_pred,
    all_x=data$S,
    all_y=data$N_recip,
    y_obs=data$N_recip[data$fracture==1],
    beta0_prior_loc=mod$coefficients[1],
    beta1_prior_loc=mod$coefficients[2],
    sigma_prior_scale = sd(data$N_recip)/2
  )
  
  stan_data$is_7_obs = rep(0,stan_data$N_obs)
  stan_data$is_7_cens = rep(1,stan_data$N_cens)
  data$is_7 = 1 - data$fracture
  
  # get reasonable priors for logistic model
  res_glm = glm(is_7 ~ S,family=binomial,data=data)
  
  stan_data$logistic_b1_loc = res_glm$coefficients[2]
  stan_data$logistic_b1_scale = abs(res_glm$coefficients[2])
  stan_data$logistic_b0_loc = res_glm$coefficients[1]
  stan_data$logistic_b0_scale = res_glm$coefficients[1]
  
  # run the stan file
  res_stan = stan('linear_model_p_runout.stan',
                  data=stan_data,
                  iter=4000)
  
  res_ext = extract(res_stan)
  
  # make sure fit looks reasonable
  ggplot(data,aes(x=untrans(N_recip),y=S)) + 
    geom_point() +
    geom_line(data=data.frame(logN=untrans(colMeans(res_ext$y_rep)),y=x_pred),
              aes(x=logN,y=y),inherit.aes=F)
  
  saveRDS(list(res_ext=res_ext,data=data,
               res_sum = summary(res_stan)$summary),paste0('model_data/',which_type,"_linear_stan_res.rds"))
  
}

stop = Sys.time()

