library(mcmcse)
library("rjson")
library(ggplot2)

setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_3/q_k2_dim3/")
#setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_10/q_k2_dim10/")
#setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_10_PC_new/q_k2_dim10/")
#setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_3_PC/q_k2_dim3/")
#prior_par <- fromJSON(file="_data_EXPprior_d3")
#iter_list <- iter_list[c(-5, -77, -55, -57, -91, -81, -40)]
#setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/JC_MODEL/jc_k2/")
iter_list <- 1:100
#mu = prior_par[1]
#lamb = prior_par[2]
#omega = prior_par[3]
#theta = prior_par[4]
#exp = "EXP_D10"
exp = "Q_D3"
#exp = "Q_D10"
#exp = "QC_D10"
#exp = "QC_D3"
#exp = "JC"
acc_rate <- function(data){
  n <- length(data)
  return(1 - sum(data[1: n - 1] == data[2: n]) / n)
}

NumRound = 100
time_MH <- c()

time_oMH <- c()
var <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 1.7) 

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


for(j in iter_list){
  show(j)
  ess_MH_a <- c()
  ess_GBS_a <- c()
  ess_oMH_a <- c()
  ess_MH_b <- c()
  ess_oMH_b <- c()
  ess_GBS_b <- c()
  
  for(cur_var in var){
    if(cur_var != 100){
      mh_file_a <- paste("__",as.character(j), "_MH_alpha",as.character(cur_var), sep = '')
    }else{
      mh_file_a <- paste("__",as.character(j), "_MH_alpha",as.character(cur_var),".0", sep = '')
    }
    
    if(file.exists(mh_file_a)){
      data <- fromJSON(file = mh_file_a)
      ess_MH_a <- c(ess_MH_a, ess(data))
    }  
    
    if(cur_var != 30){
      mh_file_b <- paste("__",as.character(j), "_MH_beta",as.character(cur_var), sep = '')
    }else{
      mh_file_b <- paste("__",as.character(j), "_MH_beta",as.character(cur_var),".0", sep = '')
    }
    if(file.exists(mh_file_b)){
      data <- fromJSON(file = mh_file_b)
      ess_MH_b <- c(ess_MH_b, ess(data))
    }
    if(cur_var != 31){
      omh_file_a <- paste("__",as.character(j), "_oMH_alpha",as.character(cur_var), sep = '')
    }else{
      omh_file_a <- paste("__",as.character(j), "_oMH_alpha",as.character(cur_var),".0", sep = '')
    }
    if(file.exists(omh_file_a)){
      data <- fromJSON(file = omh_file_a)
      ess_oMH_a <- c(ess_oMH_a, ess(data))
    }    
    if(cur_var != 31){
      omh_file_b <- paste("__",as.character(j), "_oMH_beta",as.character(cur_var), sep = '')
    }else{
      omh_file_b <- paste("__",as.character(j), "_oMH_beta",as.character(cur_var),".0", sep = '')
    }
    if(file.exists(omh_file_b)){
      data <- fromJSON(file = omh_file_b)
      ess_oMH_b <- c(ess_oMH_b, ess(data))
    }
    
    gbs_file_a <- paste("__",as.character(j), "_GBS_alpha0", sep = '')
    if(file.exists(gbs_file_a)){
      
        data <- fromJSON(file = gbs_file_a)
    ess_GBS_a <- c(ess_GBS_a, ess(data))
    }
    gbs_file_b <- paste("__",as.character(j), "_GBS_beta0", sep = '')
    if(file.exists(gbs_file_b)){
        data <- fromJSON(file = gbs_file_b)
    ess_GBS_b <- c(ess_GBS_b, ess(data))
    }
          }
  ESS_MH_a <- cbind(ESS_MH_a, ess_MH_a)
  ESS_MH_b <- cbind(ESS_MH_b, ess_MH_b)
  ESS_GBS_a <- cbind(ESS_GBS_a, ess_GBS_a)
  ESS_GBS_b <- cbind(ESS_GBS_b, ess_GBS_b)
  ESS_oMH_a <- cbind(ESS_oMH_a, ess_oMH_a)
  ESS_oMH_b <- cbind(ESS_oMH_b, ess_oMH_b)
}

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
summary_omh_persec_b <- summarize(ESS_oMH_b)
summary_mh_persec_a <- summarize(ESS_MH_a)
summary_mh_persec_b <- summarize(ESS_MH_b)
summary_gbs_persec_a <- summarize(ESS_GBS_a)
summary_gbs_persec_b <- summarize(ESS_GBS_b)


data_alpha_MH_mean <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_a[1,]))), 
                                 ess_alpha_persec = c(summary_mh_persec_a[1,]))

data_alpha_oMH_mean <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_a[1,]))),
                                  ess_alpha_persec = c(summary_omh_persec_a[1,]))

data_alpha_GBS_mean <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_a[1,]))),
                                  ess_alpha_persec = c(summary_gbs_persec_a[1,]))


maxx = max(var)
maxy = max(max(data_alpha_MH_mean[3]),max(data_alpha_GBS_mean[3]), max(data_alpha_oMH_mean[3]))
ratio.values <- (maxx)/(maxy)

s1 = 5
s2 =2
p_ALPHA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40)) + coord_fixed(ratio = ratio.values)+ ylim(0, maxy)+
  geom_point(data = data_alpha_MH_mean, aes(x = var, y =ess_alpha_persec) , shape = 0, colour = "skyblue3", size = s1, alpha = 0.8, stroke = 3) +
  geom_line(data = data_alpha_MH_mean, aes(x = var, y =ess_alpha_persec) ,colour = "skyblue3",size =s2, alpha = 0.8, linetype = 'solid') +
  
  geom_point(data = data_alpha_oMH_mean, aes(x = var, y =ess_alpha_persec) , shape = 2, colour = "darkorange", size = s1, alpha = 0.8, stroke = 3) +
  geom_line(data = data_alpha_oMH_mean, aes(x = var, y =ess_alpha_persec) ,colour = "darkorange", size = s2, alpha = 0.8, linetype = 'twodash') +
  
  geom_point(data = data_alpha_GBS_mean, aes(x = var, y =ess_alpha_persec), shape = 1, colour = "tomato1", size = s1, alpha = 0.8, stroke = 3) +
  geom_line(data = data_alpha_GBS_mean, aes(x = var, y =ess_alpha_persec) ,colour = "tomato1", size = s2, alpha = 0.8, linetype = 'longdash') +
  labs(x = expression(paste(sigma^2 ," of MH proposal") )) + labs(y =  expression(paste("ESS for ",alpha))) + theme(legend.position="none")+ 
  theme(
    axis.ticks.x=element_blank(),
    axis.ticks.y = element_blank())
p_ALPHA

setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/ess/")
filename <- paste(exp, "alpha_k2", ".pdf", sep = "")
ggsave(filename, width = 6.8, height = 6.8)




data_beta_MH_mean <- data.frame(var = var, method = factor(rep(c("MH"), each=length(summary_mh_persec_b[1,]))), 
                                ess_beta_persec = c(summary_mh_persec_b[1,]))

data_beta_oMH_mean <- data.frame(var = var, method = factor(rep(c("oMH"), each=length(summary_omh_persec_b[1,]))), 
                                 ess_beta_persec = c(summary_omh_persec_b[1,]))

data_beta_GBS_mean <- data.frame(var = var, method = factor(rep(c("GBS"), each=length(summary_gbs_persec_b[1,]))), 
                                 ess_beta_persec = c(summary_gbs_persec_b[1,]))



maxx = max(var)

maxy = max(max(data_beta_MH_mean[3]), max(data_beta_GBS_mean[3]), max(data_beta_oMH_mean[3])) * 1.1
#maxy = 0.9
ratio.values <- (maxx)/(maxy )

p_BETA <-ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40)) + coord_fixed(ratio = ratio.values) + ylim(0, maxy)
p_BETA <- p_BETA + 
  geom_point(data = data_beta_MH_mean, aes(x = var, y =ess_beta_persec) , shape = 0, colour = "skyblue3", alpha = 0.8, size = s1, stroke = 3) +
  geom_line(data = data_beta_MH_mean, aes(x = var, y =ess_beta_persec) ,colour = "skyblue3", alpha = 0.8, size = s2, linetype = 'solid') +
  geom_point(data = data_beta_oMH_mean, aes(x = var, y =ess_beta_persec) , shape = 2, colour = "darkorange", alpha = 0.8, size = s1, stroke = 3) +
  geom_line(data = data_beta_oMH_mean, aes(x = var, y =ess_beta_persec) ,colour = "darkorange", alpha = 0.8, size = s2, linetype = 'twodash') +
  geom_point(data = data_beta_GBS_mean, aes(x = var, y =ess_beta_persec) , shape = 1, colour = "tomato1", alpha = 0.8, size = s1, stroke = 3) +
  geom_line(data = data_beta_GBS_mean, aes(x = var, y =ess_beta_persec) ,colour = "tomato1", alpha = 0.8, size = s2, linetype = 'longdash') +
  labs(x = expression(paste(sigma^2 ," of MH proposal") )) + labs(y =  expression(paste("ESS for ",beta))) + theme(legend.position="none")+ 
  theme(
    axis.ticks.x=element_blank(),
    axis.ticks.y = element_blank())
p_BETA


filename <- paste(exp, "beta_k2.pdf", sep = "")
ggsave(filename, width = 6.8, height = 6.8)



