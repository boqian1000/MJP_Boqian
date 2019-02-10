setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/data/q_10/")
library(ggplot2)
load("rda_var")
load("gbs_a")
load("gbs_b")
load("mh_a")
load("mh_b")
load("omh_a")
load("omh_b")

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
#roberts and rosenthal 2004
s1 = 5 
s2 = 2
trans_col = 0.6 
ratio.display <- 1/1
maxx = 1.8
# dim3
#maxy = 7.6
#dim10

#maxy = 1.3


maxy = max(gbs_a$gbs_a_k2u)
maxy = max(maxy, max(gbs_a$gbs_a_k1.5u))
maxy = max(maxy, max(gbs_a$gbs_a_k3u))
maxy = max(maxy, max(mh_a$mh_a_k1.5u))
maxy = max(maxy, max(mh_a$mh_a_k3u))
maxy = max(maxy, max(mh_a$mh_a_k2u))
maxy = max(maxy, max(omh_a$omh_a_k1.5u))
maxy = max(maxy, max(omh_a$omh_a_k2u))
maxy = max(maxy, max(omh_a$omh_a_k3u))


ratio.values <- (maxx)/(maxy)

#library(ggplot2)
## square 0, round 1, triangle 2

trans_col = 0.6 
#roberts and rosenthal 2004
p_ALPHA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
p_ALPHA <- p_ALPHA + coord_fixed(ratio = ratio.values / ratio.display)
p_ALPHA <- p_ALPHA +   ylim(0,maxy) + xlim(0, maxx)+scale_shape(solid = FALSE) + 
  geom_point(data = mh_a, aes(x = var, y =mh_a_k1.5) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 1, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var,ymax = mh_a_k1.5u, ymin = mh_a_k1.5d), width=0.05, alpha = 0.6,colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var, y =mh_a_k1.5) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_a, aes(x = var2, y =mh_a_k2) , colour = "skyblue3", size = s1, alpha = 0.6,shape = 0, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var2, ymax = mh_a_k2u, ymin = mh_a_k2d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var2, y =mh_a_k2) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_a, aes(x = var3, y =mh_a_k3) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 2, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var3,ymax = mh_a_k3u, ymin = mh_a_k3d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var3, y =mh_a_k3) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  labs(x = expression(paste(sigma^2," of MH proposal") )) + labs(y =  expression(paste("ESS per unit time for ",alpha))) + theme(legend.position="none")

p_ALPHA
setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/new_whole_exp_fitures/")
#ggsave("q_alpha_dim5.pdf", height = 8, width = 8)
#ggsave("q_alpha_dim10.pdf", height = 8, width = 8)
#ggsave("q_alpha_dim3.pdf", height = 8, width = 8)

#ggsave("cq_alpha_dim3.pdf", height = 8, width = 8)
#ggsave("cq_alpha_dim5.pdf", height = 8, width = 8)
ggsave("mh_q_alpha_dim10.pdf", height = 8, width = 8)

maxx = 1.8
# dim3
maxy = max(gbs_b$gbs_b_k2u)
maxy = max(maxy, max(gbs_b$gbs_b_k1.5u))
maxy = max(maxy, max(gbs_b$gbs_b_k3u))
maxy = max(maxy, max(mh_b$mh_b_k1.5u))
maxy = max(maxy, max(mh_b$mh_b_k3u))
maxy = max(maxy, max(mh_b$mh_b_k2u))
maxy = max(maxy, max(omh_b$omh_b_k1.5u))
maxy = max(maxy, max(omh_b$omh_b_k2u))
maxy = max(maxy, max(omh_b$omh_b_k3u))
# dim10
#maxy = 1.65

#maxy = 1.65
ratio.values <- (maxx)/(maxy)

p_BETA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
p_BETA <- p_BETA + coord_fixed(ratio = ratio.values / ratio.display)
p_BETA <- p_BETA +   ylim(0,maxy) + xlim(0, maxx)+scale_shape(solid = FALSE) + 
  geom_point(data = mh_b, aes(x = var, y =mh_b_k1.5) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 1, stroke = 2) +
  geom_errorbar(data = mh_b, aes(x= var,ymax = mh_b_k1.5u, ymin = mh_b_k1.5d), width=0.05, alpha = 0.6,colour = "skyblue3")+
  geom_line(data = mh_b, aes(x = var, y =mh_b_k1.5) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_b, aes(x = var2, y =mh_b_k2) , colour = "skyblue3", size = s1, alpha = 0.6,shape = 0, stroke = 2) +
  geom_errorbar(data = mh_b, aes(x= var2, ymax = mh_b_k2u, ymin = mh_b_k2d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_b, aes(x = var2, y =mh_b_k2) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_b, aes(x = var3, y =mh_b_k3) , colour = "skyblue3", size = s1,alpha = 0.6, shape = 2, stroke = 2) +
  geom_errorbar(data = mh_b, aes(x= var3,ymax = mh_b_k3u, ymin = mh_b_k3d), width=0.05,alpha = 0.6, colour = "skyblue3")+
  geom_line(data = mh_b, aes(x = var3, y =mh_b_k3) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  labs(x = expression(paste(sigma^2," of MH proposal") )) + labs(y =  expression(paste("ESS per unit time for ",beta))) + theme(legend.position="none")
p_BETA

#ggsave("q_beta_dim5.pdf", height = 8, width = 8)
#ggsave("q_beta_dim10.pdf", height = 8, width = 8)
#ggsave("q_beta_dim3.pdf", height = 8, width = 8)

#ggsave("cq_beta_dim3.pdf", height = 8, width = 8)
#ggsave("cq_beta_dim5.pdf", height = 8, width = 8)
ggsave("mh_q_beta_dim10.pdf", height = 8, width = 8)



