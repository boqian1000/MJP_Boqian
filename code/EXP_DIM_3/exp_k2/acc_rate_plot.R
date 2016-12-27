library(mcmcse)
library("rjson")
library(ggplot2)
setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_DIM_10/exp_k3_new")
prior_par <- fromJSON(file="_data_EXPprior_d10")
iter_list <- iter_list[c(-5, -77, -55, -57, -91, -81, -40)]

mu = prior_par[1]
lamb = prior_par[2]
omega = prior_par[3]
theta = prior_par[4]

acc_rate <- function(data){
  n <- length(data)
  return(1 - sum(data[1: n - 1] == data[2: n]) / n)
}

NumRound = 100
time_MH <- c()
time_GBS <- c()
time_oMH <- c()
var <- c()
ess_MH_a <- c()
ESS_MH_a <- c()
ess_GBS_a <- c()
ESS_GBS_a <- c()
ess_oMH_a <- c()
ESS_oMH_a <- c()

ess_MH_persec_a <- c()
ess_oMH_persec_a <- c()
ess_GBS_persec_a <- c()

ess_MH_b <- c()
ESS_MH_b <- c()
ess_oMH_b <- c()
ESS_oMH_b <- c()

ess_GBS_b <- c()
ESS_GBS_b <- c()

ess_MH_persec_b <- c()
ess_oMH_persec_b <- c()
ess_GBS_persec_b <- c()

iter_list <- 1: (NumRound)
for(round in iter_list){
  #  if(round == 23 || round == 58|| round == 42 || round ==31) next
  mh_filename <- paste("__",as.character(round),"_timelist_mh",sep = '')
  mh_time <- fromJSON(file = mh_filename)
  
  gbs_filename <- paste("__",as.character(round),"_timelist_gbs",sep = '')
  gbs_time <- fromJSON(file = gbs_filename)
  
  omh_filename <- paste("__",as.character(round),"_timelist_omh",sep = '')
  omh_time <- fromJSON(file = omh_filename)
  
  var_filename <- paste("__",as.character(round),"var_list",sep = '')
  temp_var <- fromJSON(file = var_filename)
  
  time_MH <- cbind(time_MH, mh_time)
  time_GBS <- cbind(time_GBS, gbs_time)
  time_oMH <- cbind(time_oMH, omh_time)
  var <- cbind(var, temp_var)
}
for(j in iter_list){
  show(j)
  ess_MH_a <- c()
  ess_GBS_a <- c()
  ess_oMH_a <- c()
  ess_MH_b <- c()
  ess_oMH_b <- c()
  ess_GBS_b <- c()
  #if(j == 23 ||j == 58|| j == 42 || j ==31) next
  for(cur_var in var[,1]){
    if(cur_var != 100){
      mh_file_a <- paste("__",as.character(j), "_MH_alpha",as.character(cur_var), sep = '')
    }else{
      mh_file_a <- paste("__",as.character(j), "_MH_alpha",as.character(cur_var),".0", sep = '')
    }
    data <- fromJSON(file = mh_file_a)
    ess_MH_a <- c(ess_MH_a, acc_rate(data))
    
    if(cur_var != 30){
      mh_file_b <- paste("__",as.character(j), "_MH_beta",as.character(cur_var), sep = '')
    }else{
      mh_file_b <- paste("__",as.character(j), "_MH_beta",as.character(cur_var),".0", sep = '')
    }
    
    data <- fromJSON(file = mh_file_b)
    ess_MH_b <- c(ess_MH_b, acc_rate(data))
    
    if(cur_var != 31){
      omh_file_a <- paste("__",as.character(j), "_oMH_alpha",as.character(cur_var), sep = '')
    }else{
      omh_file_a <- paste("__",as.character(j), "_oMH_alpha",as.character(cur_var),".0", sep = '')
    }
    data <- fromJSON(file = omh_file_a)
    ess_oMH_a <- c(ess_oMH_a, acc_rate(data))
    
    if(cur_var != 31){
      omh_file_b <- paste("__",as.character(j), "_oMH_beta",as.character(cur_var), sep = '')
    }else{
      omh_file_b <- paste("__",as.character(j), "_oMH_beta",as.character(cur_var),".0", sep = '')
    }
    data <- fromJSON(file = omh_file_b)
    ess_oMH_b <- c(ess_oMH_b, acc_rate(data))
    if(cur_var != 31){
      gbs_file_a <- paste("__",as.character(j), "_GBS_alpha",as.character(cur_var), sep = '')
    }else{
      gbs_file_a <- paste("__",as.character(j), "_GBS_alpha",as.character(cur_var), ".0", sep = '')
    }
    data <- fromJSON(file = gbs_file_a)
    ess_GBS_a <- c(ess_GBS_a, acc_rate(data))
    if(cur_var != 31){
      gbs_file_b <- paste("__",as.character(j), "_GBS_beta",as.character(cur_var), sep = '')
    }else{
      gbs_file_b <- paste("__",as.character(j), "_GBS_beta",as.character(cur_var),".0", sep = '')
    }
    data <- fromJSON(file = gbs_file_b)
    ess_GBS_b <- c(ess_GBS_b, acc_rate(data))
  }
  ESS_MH_a <- cbind(ESS_MH_a, ess_MH_a)
  ESS_MH_b <- cbind(ESS_MH_b, ess_MH_b)
  ESS_GBS_a <- cbind(ESS_GBS_a, ess_GBS_a)
  ESS_GBS_b <- cbind(ESS_GBS_b, ess_GBS_b)
  ESS_oMH_a <- cbind(ESS_oMH_a, ess_oMH_a)
  ESS_oMH_b <- cbind(ESS_oMH_b, ess_oMH_b)
}

acc_rate_a <- ESS_MH_a
acc_rate_b <- ESS_MH_b
save(acc_rate_a, file = "rda_acc_MH_a")
save(acc_rate_b, file = "rda_acc_MH_b")

#save(var, file = "rda_var")
#save(time_MH, file = "rda_time_MH")
#save(time_oMH, file = "rda_time_oMH")
#save(time_GBS, file = "rda_time_GBS")
#save(ess_GBS_a, file = "rda_ess_GBS_a")
#save(ess_GBS_b, file = "rda_ess_GBS_b")
#save(ess_oMH_a, file = "rda_ess_oMH_a")
#save(ess_oMH_b, file = "rda_ess_oMH_b")


#load("rda_var")
#save(time_MH, file = "rda_time_MH")
#save(time_GBS, file = "rda_time_GBS")
#save(ess_MH_a, file = "rda_ess_MH_a")
#save(ess_MH_b, file = "rda_ess_MH_b")
#save(ess_GBS_a, file = "rda_ess_GBS_a")
#save(ess_GBS_b, file = "rda_ess_GBS_b")

#ess_oMH_persec_a <- ESS_oMH_a / time_oMH
#ess_oMH_persec_b <- ESS_oMH_b / time_oMH
#ess_MH_persec_a <- ESS_MH_a / time_MH
#ess_GBS_persec_a <- ESS_GBS_a / time_GBS
#ess_MH_persec_b <- ESS_MH_b / time_MH
#ess_GBS_persec_b <- ESS_GBS_b / time_GBS
#save(ess_MH_persec_a, file = "rda_ess_MH_persec_a")
#save(ess_MH_persec_b, file = "rda_ess_MH_persec_b")
#save(ess_oMH_persec_a, file = "rda_ess_oMH_persec_a")
#save(ess_oMH_persec_b, file = "rda_ess_oMH_persec_b")
#save(ess_GBS_persec_a, file = "rda_ess_GBS_persec_a")
#save(ess_GBS_persec_b, file = "rda_ess_GBS_persec_b")

rep_n <- dim(ess_MH_persec_a)[1]
NumRound
# mean,  min, 1st quantile, 3rd quantile meadian, max

summarize <- function(mat){
  ret <- c()
  temp <- apply(mat, 1, mean)
  ret <- c(ret, temp)
  temp <- apply(mat, 1, quantile)
  ret <- rbind(ret, temp)
  temp <- apply(mat, 1, quantile, probs = c(0.025, 0.975))
  ret <- rbind(ret, temp)
  
  return(ret)
}

summary_omh_persec_a <- summarize(ESS_oMH_a)
summary_omh_persec_b <- summarize(ESS_oMH_a)
summary_mh_persec_a <- summarize(ESS_MH_b)
summary_mh_persec_b <- summarize(ESS_MH_b)
summary_gbs_persec_a <- summarize(ESS_GBS_a)
summary_gbs_persec_b <- summarize(ESS_GBS_b)

var <- var[,1]

data_alpha_MH_mean <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_a[1,]))), 
                                 ess_alpha_persec = c(summary_mh_persec_a[1,]))
data_alpha_MH_975 <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_a[8,]))), 
                                ess_alpha_persec = c(summary_mh_persec_a[8,]))
data_alpha_MH_025 <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_a[7,]))), 
                                ess_alpha_persec = c(summary_mh_persec_a[7,]))

data_alpha_oMH_mean <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_a[1,]))),
                                  ess_alpha_persec = c(summary_omh_persec_a[1,]))
data_alpha_oMH_975 <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_a[8,]))), 
                                 ess_alpha_persec = c(summary_omh_persec_a[8,]))
data_alpha_oMH_025 <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_a[7,]))), 
                                 ess_alpha_persec = c(summary_omh_persec_a[7,]))

data_alpha_GBS_mean <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_a[1,]))), 
                                  ess_alpha_persec = c(summary_gbs_persec_a[1,]))
data_alpha_GBS_975 <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_a[8,]))), 
                                 ess_alpha_persec = c(summary_gbs_persec_a[8,]))
data_alpha_GBS_025 <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_a[7,]))), 
                                 ess_alpha_persec = c(summary_gbs_persec_a[7,]))

p_ALPHA <- ggplot()
p_ALPHA <- p_ALPHA + 
  geom_point(data = data_alpha_MH_mean, aes(x = var, y =ess_alpha_persec) , colour = "darkblue") +
  geom_line(data = data_alpha_MH_mean, aes(x = var, y =ess_alpha_persec) ,colour = "darkblue") +
  geom_point(data = data_alpha_MH_025, aes(x = var, y =ess_alpha_persec) , colour = "darkblue") +
  geom_line(data = data_alpha_MH_025, aes(x = var, y =ess_alpha_persec) ,linetype = "dotdash", colour = "darkblue") +
  geom_point(data = data_alpha_MH_975, aes(x = var, y =ess_alpha_persec) , colour = "darkblue") +
  geom_line(data = data_alpha_MH_975, aes(x = var, y =ess_alpha_persec) ,linetype = "dotdash",colour = "darkblue") +
  geom_point(data = data_alpha_oMH_mean, aes(x = var, y =ess_alpha_persec) , colour = "green") +
  geom_line(data = data_alpha_oMH_mean, aes(x = var, y =ess_alpha_persec) ,colour = "green") +
  geom_point(data = data_alpha_oMH_025, aes(x = var, y =ess_alpha_persec) , colour = "green") +
  geom_line(data = data_alpha_oMH_025, aes(x = var, y =ess_alpha_persec) ,linetype = "dotdash", colour = "green") +
  geom_point(data = data_alpha_oMH_975, aes(x = var, y =ess_alpha_persec) , colour = "green") +
  geom_line(data = data_alpha_oMH_975, aes(x = var, y =ess_alpha_persec) ,linetype = "dotdash",colour = "green") +
  geom_point(data = data_alpha_GBS_mean, aes(x = var, y =ess_alpha_persec) ,colour = "red") +
  geom_line(data = data_alpha_GBS_mean, aes(x = var, y =ess_alpha_persec) ,colour = "red") +
  geom_point(data = data_alpha_GBS_025, aes(x = var, y =ess_alpha_persec) ,colour = "red") +
  geom_line(data = data_alpha_GBS_025, aes(x = var, y =ess_alpha_persec) ,linetype = "dotdash",colour = "red") +
  geom_point(data = data_alpha_GBS_975, aes(x = var, y =ess_alpha_persec) ,colour = "red") +
  geom_line(data = data_alpha_GBS_975, aes(x = var, y =ess_alpha_persec) ,linetype = "dotdash",colour = "red") +
  labs(x = "variance proposal") + labs(y = "acceptance rate for alpha") 
p_ALPHA


data_beta_MH_mean <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_b[1,]))), 
                                ess_beta_persec = c(summary_mh_persec_b[1,]))
data_beta_MH_975 <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_b[8,]))), 
                               ess_beta_persec = c(summary_mh_persec_b[8,]))
data_beta_MH_025 <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_b[7,]))), 
                               ess_beta_persec = c(summary_mh_persec_b[7,]))

data_beta_oMH_mean <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_b[1,]))), 
                                 ess_beta_persec = c(summary_omh_persec_b[1,]))
data_beta_oMH_975 <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_b[8,]))), 
                                ess_beta_persec = c(summary_omh_persec_b[8,]))
data_beta_oMH_025 <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_b[7,]))), 
                                ess_beta_persec = c(summary_omh_persec_b[7,]))


data_beta_GBS_mean <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_b[1,]))), 
                                 ess_beta_persec = c(summary_gbs_persec_b[1,]))
data_beta_GBS_975 <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_b[8,]))), 
                                ess_beta_persec = c(summary_gbs_persec_b[8,]))
data_beta_GBS_025 <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_b[7,]))), 
                                ess_beta_persec = c(summary_gbs_persec_b[7,]))

p_BETA <- ggplot()
p_BETA <- p_BETA + 
  geom_point(data = data_beta_MH_mean, aes(x = var, y =ess_beta_persec) , colour = "darkblue") +
  geom_line(data = data_beta_MH_mean, aes(x = var, y =ess_beta_persec) ,colour = "darkblue") +
  geom_point(data = data_beta_MH_025, aes(x = var, y =ess_beta_persec) , colour = "darkblue") +
  geom_line(data = data_beta_MH_025, aes(x = var, y =ess_beta_persec) ,linetype = "dotdash", colour = "darkblue") +
  geom_point(data = data_beta_MH_975, aes(x = var, y =ess_beta_persec) , colour = "darkblue") +
  geom_line(data = data_beta_MH_975, aes(x = var, y =ess_beta_persec) ,linetype = "dotdash",colour = "darkblue") +
  geom_point(data = data_beta_oMH_mean, aes(x = var, y =ess_beta_persec) , colour = "darkorange") +
  geom_line(data = data_beta_oMH_mean, aes(x = var, y =ess_beta_persec) ,colour = "darkorange") +
  geom_point(data = data_beta_oMH_025, aes(x = var, y =ess_beta_persec) , colour = "darkorange") +
  geom_line(data = data_beta_oMH_025, aes(x = var, y =ess_beta_persec) ,linetype = "dotdash", colour = "darkorange") +
  geom_point(data = data_beta_oMH_975, aes(x = var, y =ess_beta_persec) , colour = "darkorange") +
  geom_line(data = data_beta_oMH_975, aes(x = var, y =ess_beta_persec) ,linetype = "dotdash",colour = "darkorange") +
  geom_point(data = data_beta_GBS_mean, aes(x = var, y =ess_beta_persec) ,colour = "darkred") +
  geom_line(data = data_beta_GBS_mean, aes(x = var, y =ess_beta_persec) ,colour = "darkred") +
  geom_point(data = data_beta_GBS_025, aes(x = var, y =ess_beta_persec) ,colour = "darkred") +
  geom_line(data = data_beta_GBS_025, aes(x = var, y =ess_beta_persec) ,linetype = "dotdash",colour = "darkred") +
  geom_point(data = data_beta_GBS_975, aes(x = var, y =ess_beta_persec) ,colour = "darkred") +
  geom_line(data = data_beta_GBS_975, aes(x = var, y =ess_beta_persec) ,linetype = "dotdash",colour = "darkred") +
  labs(x = "variance proposal") + labs(y = "acceptance rate for beta")
p_BETA



