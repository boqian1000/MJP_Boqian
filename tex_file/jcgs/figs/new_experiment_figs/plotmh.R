# exp model final plot
setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/data/exp_dim3/")
library(ggplot2)
load("rda_var")
#load("rda_L")
load("gbs_a")
gbs_a$var2 = var[,1] + 0.02
gbs_a$var3 = var[,1] - 0.01

load("gbs_b")
gbs_b$var2 = var[,1] + 0.02
gbs_b$var3 = var[,1] - 0.01

load("mh_a")
mh_a$var2 = var[,1] + 0.025
mh_a$var3 = var[,1] - 0.015

load("mh_b")
mh_b$var2 = var[,1] + 0.025
mh_b$var3 = var[,1] - 0.015

load("omh_a")
omh_a$var2 = var[,1] + 0.023
omh_a$var3 = var[,1] - 0.013

load("omh_b")
omh_b$var2 = var[,1] + 0.023
omh_b$var3 = var[,1] - 0.013
load("p_a")
p_a$var2 = var[,1] + 0.021
p_a$var3 = var[,1] - 0.011
load("p_b")
p_b$var2 = var[,1] + 0.021
p_b$var3 = var[,1] - 0.011


# dim 5
#r1 = 1.7 / 7
#r2 = 1.7 / 8

# dim 3
#r1 = 1.7 / 30
#r2 = 1.7 / 40
#h

# r1 = 10 / 2.5
# r2 = 10 / 2

s1 = 5 
s2 = 1.5





s1 = 5 
s2 = 2
trans_col = 0.6 
ratio.display <- 1/1
maxx = 1.8
# dim3
#maxy = 7.6
#dim10

#maxy = 1.3


maxy = 0
maxy = max(maxy, max(mh_a$mh_a_k1.5u))
maxy = max(maxy, max(mh_a$mh_a_k3u))
maxy = max(maxy, max(mh_a$mh_a_k2u))

miny = 1000
miny = min(miny, min(mh_a$mh_a_k1.5d))
miny = min(miny, min(mh_a$mh_a_k3d))
miny = min(miny, min(mh_a$mh_a_k2d))

miny = 0
ratio.values <- (maxx)/(maxy - miny)

#library(ggplot2)
## square 0, round 1, triangle 2

trans_col = 0.6 
#roberts and rosenthal 2004
p_ALPHA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
p_ALPHA <- p_ALPHA + coord_fixed(ratio = ratio.values / ratio.display)
p_ALPHA <- p_ALPHA +   ylim(miny, maxy) + xlim(0, maxx)+ scale_shape(solid = FALSE) + 
  geom_point(data = mh_a, aes(x = var, y =mh_a_k1.5) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 1, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var,ymax = mh_a_k1.5u, ymin = mh_a_k1.5d), width=0.05, alpha = 0.6,colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var, y =mh_a_k1.5) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_a, aes(x = var2, y =mh_a_k2) , colour = "skyblue3", size = s1, alpha = 0.6,shape = 0, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var2, ymax = mh_a_k2u, ymin = mh_a_k2d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var2, y =mh_a_k2) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_a, aes(x = var3, y =mh_a_k3) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 2, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var3,ymax = mh_a_k3u, ymin = mh_a_k3d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var3, y =mh_a_k3) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  labs(x = expression(paste(sigma^2," of MH proposal") )) + labs(y =  expression(paste("ESS/unit time for ",alpha))) + theme(legend.position="none")

p_ALPHA
setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/new_whole_exp_fitures/")

#ggsave("mh_exp_alpha_dim10.pdf", height = 6.8, width = 6.8)
ggsave("mh_exp_alpha_dim3.pdf", height = 6.8, width = 6.8)

maxx = 1.8
# dim3
maxy = 0
maxy = max(maxy, max(mh_b$mh_b_k1.5u))
maxy = max(maxy, max(mh_b$mh_b_k3u))
maxy = max(maxy, max(mh_b$mh_b_k2u))

miny = 1000
miny = min(miny, min(mh_b$mh_b_k1.5d))
miny = min(miny, min(mh_b$mh_b_k3d))
miny = min(miny, min(mh_b$mh_b_k2d))
miny = 0

ratio.values <- (maxx)/(maxy - miny)


# dim10
#maxy = 1.65

#maxy = 1.65
ratio.values <- (maxx)/(maxy)

p_BETA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
p_BETA <- p_BETA + coord_fixed(ratio = ratio.values / ratio.display)
p_BETA <- p_BETA +   ylim(miny, maxy) + xlim(0, maxx) + scale_shape(solid = FALSE) + 
  geom_point(data = mh_b, aes(x = var, y =mh_b_k1.5) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 1, stroke = 2) +
  geom_errorbar(data = mh_b, aes(x= var,ymax = mh_b_k1.5u, ymin = mh_b_k1.5d), width=0.05, alpha = 0.6,colour = "skyblue3")+
  geom_line(data = mh_b, aes(x = var, y =mh_b_k1.5) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_b, aes(x = var2, y =mh_b_k2) , colour = "skyblue3", size = s1, alpha = 0.6,shape = 0, stroke = 2) +
  geom_errorbar(data = mh_b, aes(x= var2, ymax = mh_b_k2u, ymin = mh_b_k2d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_b, aes(x = var2, y =mh_b_k2) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_b, aes(x = var3, y =mh_b_k3) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 2, stroke = 2) +
  geom_errorbar(data = mh_b, aes(x= var3,ymax = mh_b_k3u, ymin = mh_b_k3d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_b, aes(x = var3, y =mh_b_k3) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  labs(x = expression(paste(sigma^2," of MH proposal") )) + labs(y =  expression(paste("ESS/unit time for ",beta))) + theme(legend.position="none")
p_BETA


#ggsave("mh_exp_beta_dim10.pdf", height = 6.8, width = 6.8)
ggsave("mh_exp_beta_dim3.pdf", height = 6.8, width = 6.8)

