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


s1 = 5 
s1 = 5
s2 = 2

trans_col = 0.6 
ratio.display <- 1/1
maxy = max(gbs_a$gbs_a_k2u)
maxy = max(maxy, max(mh_a$mh_a_k2u))
maxy = max(maxy, max(omh_a$omh_a_k2u))

# dim3
maxx = 1.8
#dim10
#maxy = 1.3
ratio.values <- (maxx)/(maxy)



trans_col = 0.6 
#roberts and rosenthal 2004
p_ALPHA <- ggplot() + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
p_ALPHA <- p_ALPHA + coord_fixed(ratio = ratio.values / ratio.display)
p_ALPHA <- p_ALPHA + ylim(0,maxy) + xlim(0, maxx) + scale_shape(solid = FALSE) + 
  geom_point(data = gbs_a, aes(x = var2, y =gbs_a_k2), shape = 1 , colour = "tomato1", alpha = trans_col, size = s1, stroke = 2) +
  geom_errorbar(data = gbs_a, aes(x= var2, ymax = gbs_a_k2u, ymin = gbs_a_k2d),alpha = trans_col, size = 1, width=0.05, colour ="tomato1") +
  geom_line(data = gbs_a, aes(x = var2, y =gbs_a_k2) ,colour = "tomato1", size = s2, alpha = trans_col, linetype = "longdash") +
  geom_point(data = mh_a, aes(x = var2, y =mh_a_k2), shape = 0 , colour ="skyblue3" , alpha = trans_col, size = s1, stroke = 2) +
  geom_errorbar(data = mh_a, aes(x= var2, ymax = mh_a_k2u, ymin = mh_a_k2d), width=0.05, size = 1,alpha = trans_col, colour = "skyblue3") +
  geom_line(data = mh_a, aes(x = var2, y =mh_a_k2) ,colour = "skyblue3", size = s2, alpha = trans_col, linetype = "solid") +
  geom_point(data = omh_a, aes(x = var2, y =omh_a_k2), shape = 2 , colour = "darkorange", alpha = trans_col, size = s1, stroke = 2) +
  geom_errorbar(data = omh_a, aes(x= var2,ymax = omh_a_k2u, ymin = omh_a_k2d), width=0.05, size = 1, colour = "darkorange", alpha = trans_col) +
  geom_line(data = omh_a, aes(x = var2, y =omh_a_k2) ,colour = "darkorange", size = s2, alpha = trans_col, linetype = "twodash") +
  labs(x = expression(paste(sigma^2," of MH proposal") )) + labs(y =  expression(paste("ESS per unit time for ",alpha))) + theme(legend.position="none") + 
  theme(
    axis.ticks.x=element_blank(),
    axis.ticks.y = element_blank())# +
#  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) 

p_ALPHA
setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/NEWFIGURES/")
ggsave("jc_alpha_k2.pdf", width = 8, height = 8)
