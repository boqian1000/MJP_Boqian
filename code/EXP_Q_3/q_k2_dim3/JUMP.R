library(rjson)
setwd( "/Users/Isaac_Zhang/Research/MCMC/simulation_seed/EXP_DIM_3/exp_k2/")
library(ggplot2)
#q
#r = 22 / 320
# exp
r = 50 / 65
s1 = 2
transit <- fromJSON(file = "transit_time_RESULT")
N = length(transit[[1]])
maxx <- c()
minn <- c()
for(j in 1:N){
  maxt = transit[[1]][j]
  mint = transit[[1]][j]
  for(i in 1 : 10){
    if(transit[[i]][j] > maxt) maxt = transit[[i]][j]
    if(transit[[i]][j] < mint) mint = transit[[i]][j]
  }
  maxx <- c(maxx, maxt)
  minn <- c(minn, mint)
}
p <- ggplot()
p <- p + coord_fixed(r) + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=20, face ="bold")) +
  geom_point(x = 1:N, y =transit[[1]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[1]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[2]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[2]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[3]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[3]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[4]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[4]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[5]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[5]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[6]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[6]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[7]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[7]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[8]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[8]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[9]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[9]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =transit[[10]] , colour = "darkblue") +
  geom_line(aes(x = 1:N, y =transit[[10]]) ,colour = "darkblue", size = s1)+
  geom_point(x = 1:N, y =maxx , colour = "black") +
  geom_line(aes(x = 1:N, y =maxx) ,colour = "black", size = s1 + 1)+
  geom_point(x = 1:N, y =minn , colour = "black") +
  geom_line(aes(x = 1:N, y =minn) ,colour = "black", size = s1 + 1)+
  labs(x = "MCMC iteration number") + labs(y = "Number of transitions in MJP path") 
p
