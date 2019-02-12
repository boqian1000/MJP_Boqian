library(mcmcse)
library("rjson")
library(ggplot2)

f <- function(path, j, cur_var, d, name='', save = False){
  setwd( path )
  GBS_file_a <- paste("__",as.character(j), "_GBS_alpha",as.character(0), sep = '')
  MH_file_a <- paste("__",as.character(j), "_MH_alpha",as.character(cur_var), sep = '')
  gbsdata <- fromJSON(file = GBS_file_a)
  mhdata <- fromJSON(file = MH_file_a)
  
  
  id <- (1:10000) * 1 
  print( ks.test(gbsdata[id], mhdata[id]) )
  
  id2 <- (1:1000) * 10
  print( ks.test(gbsdata[id2], mhdata[id2]) )
  
  data <- data.frame(alpha = c(gbsdata, mhdata),method = c(rep("GBS",length(gbsdata)), rep("MH",length(gbsdata))), color = c(rep("tomato1",length(gbsdata)), rep("skyblue3",length(gbsdata)))) 
  
  p1 <- ggplot(data, aes(x=alpha, col = color, linetype = method))#+ coord_fixed(ratio = r1)
  p1 <- p1 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
  p1 <- p1 + geom_line(stat= "density", alpha = .8,  size = 2)+ theme(legend.position="none")+
    labs(x = expression(alpha), y = expression( paste("P(", alpha, "|", X, ")" )) )+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p1 <- p1 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/QC_ks/")

  show(p1)
  filename1 <- paste(name, "hist", as.character(j), cur_var, as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename1, width = 6.8, height = 6.8)
  }
  
  
  # trace plot
  data_gbs = data[data$method == 'GBS', ]
  maxx = dim(data_gbs)[1]
  maxy = max(data_gbs$alpha)
  miny = min(data_gbs$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p2 <- ggplot(data_gbs, aes((1: length(alpha)), alpha))
  p2 <- p2 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p2 <- p2 + geom_point(alpha=0.8, col = "tomato1") +  labs(x = "Iteration") +
    labs(y = expression(alpha))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p2 <- p2 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  show(p2)

  filename2GBS <- paste(name, "traceGBS", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename2GBS, width = 6.8, height = 6.8)
  }
  
  data_mh = data[data$method == 'MH', ]
  maxx = dim(data_mh)[1]
  maxy = max(data_mh$alpha)
  miny = min(data_mh$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p3 <- ggplot(data_mh, aes((1: length(alpha)), alpha))
  p3 <- p3 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p3 <- p3 + geom_point(alpha=0.8, col = "skyblue3") +  labs(x = "Iteration") +
    labs(y = expression(alpha))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p3 <- p3 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))
  
  show(p3)

  filenameMH <- paste(name, "traceMH", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filenameMH, width = 6.8, height = 6.8)
  }
  
  
  bacf <- acf(gbsdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) + coord_fixed(ratio = ratio.values)+
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0), color = "tomato1", alpha = 0.8) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  q <- q + labs(y = "Autocorr", x = "Lag") #+ theme(axis.title.y = element_text(vjust = 0.9)) 
  q <- q + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))  +theme(
    axis.ticks.x=element_blank(),
    axis.ticks.y = element_blank())
  
  show(q)
  
  filename3 <- paste(name, "gbsacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename3, width = 6.8, height = 6.8)
  }
  bacf <- acf(mhdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  ratio.values <- (maxx)/(maxy)
  
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) + coord_fixed(ratio = ratio.values)+
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0), color = "skyblue3", alpha = 0.8) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  q <- q + labs(y = "Autocorr", x = "Lag") #+ theme(axis.title.y = element_text(vjust = 0.9)) 
  q <- q + theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)), axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)))  +theme(
    axis.ticks.x=element_blank(),
    axis.ticks.y = element_blank())
  
  show(q)
  
  filename4 <- paste(name, "mhacf", as.character(j), cur_var,as.character(d) ,".pdf", sep = "_")
  if(save){
    ggsave(filename4, width = 6.8, height = 6.8)
  }
}
rej_rate <- function(data){
  n <- length(data)
  return(sum(data[1: n - 1] == data[2: n]) / n)
}

#in paper
path <- "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_10_PC_new/q_k2_dim10/"
f(path, 4, 0.3,10, 'qc', TRUE)


path <- "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_3_PC/q_k2_dim3"

f(path, 6, 0.3,3, 'qc', TRUE)




f(path, file, 1, 0.1,3, 'qc')


f(path, file, 11, 0.2,3, 'qc')

#

#


f(path, file, 11,  0.5,3, 'qc')


f(path, file, 1, 0.8,3, 'qc')


f(path, file, 1, 1.2,3, 'qc')

f(path, file, 6, 1.7,3, 'qc')


path <- "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_10_PC_new/q_k2_dim10/"
file <- "_data_EXPprior_immi_d10"

f(path, file, 9, 0.1,10, 'qc')


f(path, file, 4, 0.2,10, 'qc')
#
f(path, 4, 0.3,10, 'qc', TRUE)
#f(path, file, 9, 0.3,10, 'qc')


f(path, file, 3, 0.5,10, 'qc')


f(path, file, 61, 0.8, 10, 'qc')

f(path, file, 25, 1.2, 10, 'qc')

f(path, file, 28, 1.7, 10, 'qc')



