library(rjson)
#setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_seed/distribution_exp/")

setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_seed/EXP_Q_3/q_k2_dim3/")

library(ggplot2)
library(mcmcse)
alpha_list = fromJSON(file = "alpha1")
beta_list = fromJSON(file = "beta1")
loc_alpha_list = fromJSON(file = "loc_alpha1")
loc_beta_list = fromJSON(file = "loc_beta1")
rate_alpha_list = fromJSON(file = "rate_alpha1")
rate_beta_list = fromJSON(file = "rate_beta1")
method <- rep("1",5000)
alpha_list = alpha_list[5001: 10000]
beta_list = beta_list[5001: 10000]
loc_alpha_list = loc_alpha_list[5001: 10000]
loc_beta_list = loc_beta_list[5001: 10000]
rate_alpha_list = rate_alpha_list[5001: 10000]
rate_beta_list = rate_beta_list[5001: 10000]

alpha_list2 = fromJSON(file = "alpha2")
beta_list2 = fromJSON(file = "beta2")
loc_alpha_list2 = fromJSON(file = "loc_alpha2")
loc_beta_list2 = fromJSON(file = "loc_beta2")
rate_alpha_list2 = fromJSON(file = "rate_alpha2")
rate_beta_list2 = fromJSON(file = "rate_beta2")
method2 <- rep("2",5000)
r1 = 6.1/1.1
data <- data.frame(alpha = c(alpha_list,alpha_list2), beta = c(beta_list, beta_list2), method = c(method, method2))
alpha <- ggplot(data, aes(x=alpha, fill=method))+ coord_fixed(ratio = r1) + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face ="bold"))
alpha <- alpha + geom_density(alpha=.3)+ theme(legend.position="none")
alpha
r2 = 6.3 / 1.7
beta <- ggplot(data, aes(x=beta, fill=method))+ coord_fixed(ratio = r2) + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face ="bold"))
beta <- beta + geom_density(alpha=.3)+ theme(legend.position="none")
beta
len = length(alpha_list)
diff_alpha <- alpha_list[1: len - 1] - alpha_list[2:len]
tv_alpha <- sum(diff_alpha * diff_alpha) / len
tv_alpha 
 
diff_alpha2 <- alpha_list2[1: len - 1] - alpha_list2[2:len]
tv_alpha2 <- sum(diff_alpha2 * diff_alpha2) / len
tv_alpha2 
 
 
diff_beta <- beta_list[1: len - 1] - beta_list[2:len]
tv_beta <- sum(diff_beta * diff_beta) / len
tv_beta 
 
diff_beta2 <- beta_list2[1: len - 1] - beta_list2[2:len]
tv_beta2 <- sum(diff_beta2 * diff_beta2) / len
tv_beta2 
 
r3 = 6.1 / 1.1
alpha2 <- ggplot(data, aes(x=alpha, fill=method, y= ..density..))+ coord_fixed(ratio = r3) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face ="bold")) + theme_bw()
alpha2 <- alpha2 + geom_histogram(binwidth=.2, alpha=.5, position="identity")+ theme(legend.position="none")
alpha2

r4 = 6.2 / 1.55
beta2 <- ggplot(data, aes(x=beta, fill=method, y= ..density..))+ coord_fixed(ratio = r4) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face ="bold")) + theme_bw()
beta2 <- beta2 + geom_histogram(binwidth=.2, alpha=.5, position="identity")+ theme(legend.position="none")
beta2

 