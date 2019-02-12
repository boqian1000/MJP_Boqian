setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_result/data/jc")
library(ggplot2)
load("rda_var")

load("gbs_a")
load("mh_a")
load("omh_a")

load("gbs_a")
gbs_a$var2 = var[,1] + 0.02
gbs_a$var3 = var[,1] - 0.01

load("mh_a")
mh_a$var2 = var[,1] + 0.025
mh_a$var3 = var[,1] - 0.015

load("omh_a")
omh_a$var2 = var[,1] + 0.023
omh_a$var3 = var[,1] - 0.013
# q 10
#r1 = 1.7 / 1.1
#r2 = 1.7 / 1.25
# q 5
#r1 = 1.7 / 3
#r2 = 1.7 / 3.3
# q 3
#r1 = 1.7 / 12
#s1 = 5 
#s2 = 1.5
#trans_col = 0.6 
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

trans_col = 0.8
#roberts and rosenthal 2004
p_ALPHA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ scale_shape(solid = FALSE)
p_ALPHA <- p_ALPHA + coord_fixed(ratio = ratio.values / ratio.display)
p_ALPHA <- p_ALPHA +   ylim(0,maxy) + xlim(0, maxx)+
  geom_point(data = mh_a, aes(x = var, y =mh_a_k1.5) , colour = "skyblue3", size = s1,alpha = 0.8, shape = 1, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var,ymax = mh_a_k1.5u, ymin = mh_a_k1.5d), width=0.05, alpha = 0.8,colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var, y =mh_a_k1.5) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_a, aes(x = var2, y =mh_a_k2) , colour = "skyblue3", size = s1, alpha = 0.8,shape = 0, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var2, ymax = mh_a_k2u, ymin = mh_a_k2d), width=0.05,alpha = 0.8, colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var2, y =mh_a_k2) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  geom_point(data = mh_a, aes(x = var3, y =mh_a_k3) , colour = "skyblue3", size = s1,alpha = 0.8, shape = 2, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var3,ymax = mh_a_k3u, ymin = mh_a_k3d), width=0.05,alpha = 0.8, colour = "skyblue3")+
  geom_line(data = mh_a, aes(x = var3, y =mh_a_k3) ,colour = "skyblue3", size = s2, alpha = trans_col) +
  labs(x = expression(paste(sigma^2," of MH proposal") )) + labs(y =  expression(paste("ESS/unit time for ",alpha))) + theme(legend.position="none")

p_ALPHA
setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/new_whole_exp_fitures/")
ggsave("mh_jc_alpha.pdf", height = 6.8, width = 6.8)
