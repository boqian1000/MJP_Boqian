library(mcmcse)
library("rjson")
library(ggplot2)

f <- function(path, file, j, cur_var, d, name=''){
  setwd( path )
  prior_par <- fromJSON(file = file)
  mu = prior_par[1]
  lamb = prior_par[2]
  omega = prior_par[3]
  theta = prior_par[4]
  GBS_file_a <- paste("__",as.character(j), "_GBS_alpha",as.character(0), sep = '')
  MH_file_a <- paste("__",as.character(j), "_MH_alpha",as.character(cur_var), sep = '')
  gbsdata <- fromJSON(file = GBS_file_a)
  mhdata <- fromJSON(file = MH_file_a)
  length(gbsdata)
  
  
  id <- (1:10000) * 1 
  print( ks.test(gbsdata[id], mhdata[id]) )
  
  id2 <- (1:1000) * 10
  print( ks.test(gbsdata[id2], mhdata[id2]) )
  
  data <- data.frame(alpha = c(gbsdata, mhdata),method = c(rep("GBS",length(gbsdata)), rep("MH",length(gbsdata)))) 
  
  p1 <- ggplot(data, aes(x=alpha, col = method, linetype = method))#+ coord_fixed(ratio = r1)
  p1 <- p1 + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=20))
  p1 <- p1 + geom_line(stat= "density", alpha = .75,  size = 2)+ theme(legend.position="none")+
    labs(x = expression(alpha))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/QC_ks/")
  p1
  filename1 <- paste(name, "hist", as.character(j), cur_var, as.character(d), ".pdf", sep = "_")
  ggsave(filename1, width = 6, height = 6)
  show(p1)
  
  # trace plot
  maxx = dim(data)[1]
  maxy = max(data$alpha)
  miny = min(data$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p2 <- ggplot(data, aes((1: length(alpha)), alpha, col = method))
  p2 <- p2 + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=20))+ coord_fixed(ratio = ratio.values)
  p2 <- p2 + geom_point() +  labs(x = "Iteration") +
    labs(y = expression(alpha))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(legend.position="none")
  p2
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/QC_ks/")
  filename2 <- paste(name, "trace", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  ggsave(filename2)
  
  
  
  bacf <- acf(gbsdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=30))+ coord_fixed(ratio = ratio.values)
  q
  filename3 <- paste(name, "gbsacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  ggsave(filename3)
  
  bacf <- acf(mhdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=20), axis.title=element_text(size=30))+ coord_fixed(ratio = ratio.values)
  q
  filename4 <- paste(name, "mhacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  ggsave(filename4)
  
}

rej_rate <- function(data){
  n <- length(data)
  return(sum(data[1: n - 1] == data[2: n]) / n)
}

path <- "/Users/Isaac_Zhang/Research/MCMC/simulation_result/EXP_Q_3_PC/q_k2_dim3"
file <- "_data_EXPprior_immi_d3"



f(path, file, 1, 0.1,3, 'qc')


f(path, file, 11, 0.2,3, 'qc')

#
f(path, file, 6, 0.3,3, 'qc')

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
f(path, file, 4, 0.3,10, 'qc')
#f(path, file, 9, 0.3,10, 'qc')


f(path, file, 3, 0.5,10, 'qc')


f(path, file, 61, 0.8, 10, 'qc')

f(path, file, 25, 1.2, 10, 'qc')

f(path, file, 28, 1.7, 10, 'qc')



