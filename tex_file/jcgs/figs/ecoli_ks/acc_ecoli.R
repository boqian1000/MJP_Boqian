library(mcmcse)
library("rjson")
library(ggplot2)
library(scatterplot3d)
setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step/")

acc_rate <- function(data){
  n <- length(data)
  return(1 - sum(data[1: n - 1] == data[2: n]) / n)
}

collect <- function(var_b){
  essAlpha.list <- c()
  essBeta.list <- c()
  essL1.list <- c()
  essL2.list <- c()
  for( round in 1: 100){
    if(round %% 10 == 1){
      #            print(round)
    }
    time_filename <-  paste("__" , as.character(round) , "_timelist_mh_",as.character(var_b),sep = '')
    if(!file.exists(time_filename)){
      next
    }
    alpha_filename <- paste("__" , as.character(round) , "_MH_alpha_",as.character(var_b),sep = '')
    beta_filename <- paste("__" , as.character(round) , "_MH_beta_",as.character(var_b),sep = '')
    l1_filename <- paste("__" , as.character(round) , "_MH_lamb1_",as.character(var_b),sep = '')
    l2_filename <- paste("__" , as.character(round) , "_MH_lamb2_",as.character(var_b),sep = '')
    
    time <- fromJSON(file = time_filename)
    time <- time[2]
    alpha <- fromJSON(file = alpha_filename)
    beta <- fromJSON(file = beta_filename)
    l1 <- fromJSON(file = l1_filename)
    l2 <- fromJSON(file = l2_filename)
    essAlpha.list <- c(essAlpha.list, acc_rate(alpha))
    essBeta.list <- c(essBeta.list, acc_rate(beta))
    essL1.list <- c(essAlpha.list, acc_rate(l1))
    essL2.list <- c(essBeta.list, acc_rate(l2))
  }
  return(c(mean(essAlpha.list), mean(essBeta.list), mean(essL1.list), mean(essL2.list), quantile(essAlpha.list, probs=0.025),
           quantile(essBeta.list, probs=0.025), quantile(essL1.list, probs=0.025), quantile(essL2.list, probs=0.025),
           quantile(essAlpha.list, probs=0.975),
           quantile(essBeta.list, probs=0.975), quantile(essL1.list, probs=0.975), quantile(essL2.list, probs=0.975)
  ))
}

collectAll <- function(){
  
  var_bs <- 0:6
  ess.alpha.list <- c() 
  ess.beta.list <- c() 
  ess.l1.list <- c() 
  ess.l2.list <- c() 
  ess.alpha.listu <- c() 
  ess.beta.listu <- c() 
  ess.l1.listu <- c() 
  ess.l2.listu <- c() 
  ess.alpha.listd <- c() 
  ess.beta.listd <- c() 
  ess.l1.listd <- c() 
  ess.l2.listd <- c() 
  for(var_b in var_bs){
    ess <- collect(var_b)
    ess.alpha.list <- c(ess.alpha.list, ess[1])
    ess.beta.list <- c(ess.beta.list, ess[2])
    ess.l1.list <- c(ess.l1.list, ess[3])
    ess.l2.list <- c(ess.l2.list, ess[4])
    ess.alpha.listd <- c(ess.alpha.listd, ess[5])
    ess.beta.listd <- c(ess.beta.listd, ess[6])
    ess.l1.listd <- c(ess.l1.listd, ess[7])
    ess.l2.listd <- c(ess.l2.listd, ess[8])
    ess.alpha.listu <- c(ess.alpha.listu, ess[9])
    ess.beta.listu <- c(ess.beta.listu, ess[10])
    ess.l1.listu <- c(ess.l1.listu, ess[11])
    ess.l2.listu <- c(ess.l2.listu, ess[12])
    print(ess)
    print(paste("var_b", as.character(var_b)))
    
  }
  return(list(a=ess.alpha.list, b=ess.beta.list, l1=ess.l1.list, l2=ess.l2.list,ad=ess.alpha.listd, bd=ess.beta.listd, l1d=ess.l1.listd, l2d=ess.l2.listd,au=ess.alpha.listu, bu=ess.beta.listu, l1u=ess.l1.listu, l2u=ess.l2.listu))
}

ret <- collectAll()
save(ret, file="retacc.RData")
load("retacc.RData")
ess.alpha <- ret[[1]]
ess.beta <- ret[[2]]
ess.l1 <- ret[[3]]
ess.l2 <- ret[[4]]

k <- c(.1, .3, .5, 1., 3., 5., 10)
plot(k, ess.alpha, type='l')
plot(k, ess.beta, type='l')
plot(k, ess.l1, type='l')
plot(k, ess.l2, type='l')




data <- data.frame(k = k, acc_rate = ess.alpha)

maxx = 11
maxy = max(ess.alpha)
ratio.values <- (maxx)/(maxy)

s1 = 5
s2 =2
p_ALPHA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=30), axis.title=element_text(size=30)) + coord_fixed(ratio = ratio.values)+ ylim(0, maxy)+ xlim(0, 11)+
  geom_point(data = data, aes(x = k, y =acc_rate) , colour = "blue", size = s1, alpha = 0.6) +
  geom_line(data = data, aes(x = k, y =acc_rate) ,colour = "blue",size =s2, alpha = 0.6, linetype = "dashed" ) +
  labs(x = expression(paste("Scale of MH proposal") )) + labs(y =  expression(paste("Acceptance rate for ", alpha))) + theme(legend.position="none")
p_ALPHA

setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/acc/")
filename <- paste("ecoli", "alpha_k2", ".pdf", sep = "")
ggsave(filename, width = 6, height = 6)
