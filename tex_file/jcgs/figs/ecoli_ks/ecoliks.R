library(mcmcse)
library("rjson")
library(ggplot2)
j = 1
cur_var = 1
d= 0
#burning_size = 2000
# SIZE 20000
f_a <- function(j, cur_var,  name='', save = False){
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step/")
  MH_file_a <- paste("__",as.character(j), "_MH_alpha_",as.character(cur_var), sep = '')
  mhdata <- fromJSON(file = MH_file_a)
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step_gbs")
  GBS_file_a <- paste("__",as.character(j), "_GBS_alpha_", sep = '')
  gbsdata <- fromJSON(file = GBS_file_a)
  #bmh <- fromJSON(file = "__2_burnin_MH_alpha_0")
  #length(bmh)
  id <- (1:10000) * 1 
  print( ks.test(gbsdata[id], mhdata[id]) )
  
  id2 <- (1:1000) * 10
  print( ks.test(gbsdata[id2], mhdata[id2]) )
  

  data <- data.frame(alpha = c(gbsdata, mhdata),method = c(rep("GBS",length(gbsdata)), rep("MH",length(gbsdata)))) 

  p1 <- ggplot(data, aes(x=alpha, col = method, linetype = method))#+ coord_fixed(ratio = r1)
  p1 <- p1 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
  p1 <- p1 + geom_line(stat= "density", alpha = .75,  size = 2)+ theme(legend.position="none")+
    labs(x = expression(alpha), y = expression( paste("P(", alpha, "|", X, ")" )) )+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/ecoli_ks/")
  p1
  show(p1)
  filename1 <- paste(name, "alphahist", as.character(j), cur_var, as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename1, width = 6, height = 6)
  }
  
  
  # trace plot
  data_gbs = data[data$method == 'GBS', ]
  maxx = dim(data_gbs)[1]
  maxy = max(data_gbs$alpha)
  miny = min(data_gbs$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p2 <- ggplot(data_gbs, aes((1: length(alpha)), alpha))
  p2 <- p2 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p2 <- p2 + geom_point(alpha=0.6, col = "firebrick1") +  labs(x = "Iteration") +
    labs(y = expression(alpha))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p2
  show(p2)
  filename2GBS <- paste(name, "alphatraceGBS", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename2GBS, width = 6, height = 6)
  }
  
  data_mh = data[data$method == 'MH', ]
  maxx = dim(data_mh)[1]
  maxy = max(data_mh$alpha)
  miny = min(data_mh$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p3 <- ggplot(data_mh, aes((1: length(alpha)), alpha))
  p3 <- p3 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p3 <- p3 + geom_point(alpha=0.6, col = "skyblue3") +  labs(x = "Iteration") +
    labs(y = expression(alpha))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p3
  
  show(p3)
  
  filenameMH <- paste(name, "alphatraceMH", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filenameMH, width = 6, height = 6)
  }
  
  
  bacf <- acf(gbsdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  q
  show(q)
  filename3 <- paste(name, "alphagbsacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename3, width = 6, height = 6)
  }
  bacf <- acf(mhdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  
  show(q)
  filename4 <- paste(name, "alphamhacf", as.character(j), cur_var,as.character(d) ,".pdf", sep = "_")
  if(save){
    ggsave(filename4, width = 6, height = 6)
  }
}

f_b <- function(j, cur_var,  name='', save = False){
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step/")
  MH_file_a <- paste("__",as.character(j), "_MH_beta_",as.character(cur_var), sep = '')
  mhdata <- fromJSON(file = MH_file_a)
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step_gbs")
  GBS_file_a <- paste("__",as.character(j), "_GBS_beta_", sep = '')
  gbsdata <- fromJSON(file = GBS_file_a)
  #bmh <- fromJSON(file = "__2_burnin_MH_alpha_0")
  #length(bmh)
  id <- (1:10000) * 1 
  print( ks.test(gbsdata[id], mhdata[id]) )
  
  id2 <- (1:1000) * 10
  print( ks.test(gbsdata[id2], mhdata[id2]) )
  
  
  data <- data.frame(alpha = c(gbsdata, mhdata),method = c(rep("GBS",length(gbsdata)), rep("MH",length(gbsdata)))) 
  
  p1 <- ggplot(data, aes(x=alpha, col = method, linetype = method))#+ coord_fixed(ratio = r1)
  p1 <- p1 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
  p1 <- p1 + geom_line(stat= "density", alpha = .75,  size = 2)+ theme(legend.position="none")+
    labs(x = expression(alpha), y = expression( paste("P(", beta, "|", X, ")" )) )+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/ecoli_ks/")
  p1
  show(p1)
  filename1 <- paste(name, "betahist", as.character(j), cur_var, as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename1, width = 6, height = 6)
  }
  
  
  # trace plot
  data_gbs = data[data$method == 'GBS', ]
  maxx = dim(data_gbs)[1]
  maxy = max(data_gbs$alpha)
  miny = min(data_gbs$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p2 <- ggplot(data_gbs, aes((1: length(alpha)), alpha))
  p2 <- p2 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p2 <- p2 + geom_point(alpha=0.6, col = "firebrick1") +  labs(x = "Iteration") +
    labs(y = expression(beta))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p2
  show(p2)
  filename2GBS <- paste(name, "betatraceGBS", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename2GBS, width = 6, height = 6)
  }
  
  data_mh = data[data$method == 'MH', ]
  maxx = dim(data_mh)[1]
  maxy = max(data_mh$alpha)
  miny = min(data_mh$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p3 <- ggplot(data_mh, aes((1: length(alpha)), alpha))
  p3 <- p3 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p3 <- p3 + geom_point(alpha=0.6, col = "skyblue3") +  labs(x = "Iteration") +
    labs(y = expression(beta))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p3
  
  show(p3)
  
  filenameMH <- paste(name, "betatraceMH", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filenameMH, width = 6, height = 6)
  }
  
  
  bacf <- acf(gbsdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  q
  show(q)
  filename3 <- paste(name, "betagbsacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename3, width = 6, height = 6)
  }
  bacf <- acf(mhdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  
  show(q)
  filename4 <- paste(name, "betamhacf", as.character(j), cur_var,as.character(d) ,".pdf", sep = "_")
  if(save){
    ggsave(filename4, width = 6, height = 6)
  }
}


f_l1 <- function(j, cur_var,  name='', save = False){
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step/")
  MH_file_a <- paste("__",as.character(j), "_MH_lamb1_",as.character(cur_var), sep = '')
  mhdata <- fromJSON(file = MH_file_a)
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step_gbs")
  GBS_file_a <- paste("__",as.character(j), "_GBS_lamb1_", sep = '')
  gbsdata <- fromJSON(file = GBS_file_a)
  #bmh <- fromJSON(file = "__2_burnin_MH_alpha_0")
  #length(bmh)
  id <- (1:10000) * 1 
  print( ks.test(gbsdata[id], mhdata[id]) )
  
  id2 <- (1:1000) * 10
  print( ks.test(gbsdata[id2], mhdata[id2]) )
  
  
  data <- data.frame(alpha = c(gbsdata, mhdata),method = c(rep("GBS",length(gbsdata)), rep("MH",length(gbsdata)))) 
  
  p1 <- ggplot(data, aes(x=alpha, col = method, linetype = method))#+ coord_fixed(ratio = r1)
  p1 <- p1 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
  p1 <- p1 + geom_line(stat= "density", alpha = .75,  size = 2)+ theme(legend.position="none")+
    labs(x = expression(alpha), y = expression( paste("P(", lambda, "1|", X, ")" )) )+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/ecoli_ks/")
  p1
  show(p1)
  filename1 <- paste(name, "l1hist", as.character(j), cur_var, as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename1, width = 6, height = 6)
  }
  
  
  # trace plot
  data_gbs = data[data$method == 'GBS', ]
  maxx = dim(data_gbs)[1]
  maxy = max(data_gbs$alpha)
  miny = min(data_gbs$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p2 <- ggplot(data_gbs, aes((1: length(alpha)), alpha))
  p2 <- p2 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p2 <- p2 + geom_point(alpha=0.6, col = "firebrick1") +  labs(x = "Iteration") +
    labs(y = paste(lambda, "1"))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p2
  show(p2)
  filename2GBS <- paste(name, "l1traceGBS", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename2GBS, width = 6, height = 6)
  }
  
  data_mh = data[data$method == 'MH', ]
  maxx = dim(data_mh)[1]
  maxy = max(data_mh$alpha)
  miny = min(data_mh$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p3 <- ggplot(data_mh, aes((1: length(alpha)), alpha))
  p3 <- p3 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p3 <- p3 + geom_point(alpha=0.6, col = "skyblue3") +  labs(x = "Iteration") +
    labs(y = paste(lambda, "1"))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p3
  
  show(p3)
  
  filenameMH <- paste(name, "l1traceMH", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filenameMH, width = 6, height = 6)
  }
  
  
  bacf <- acf(gbsdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  q
  show(q)
  filename3 <- paste(name, "l1gbsacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename3, width = 6, height = 6)
  }
  bacf <- acf(mhdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  
  show(q)
  filename4 <- paste(name, "l1mhacf", as.character(j), cur_var,as.character(d) ,".pdf", sep = "_")
  if(save){
    ggsave(filename4, width = 6, height = 6)
  }
}

f_l2 <- function(j, cur_var,  name='', save = False){
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step/")
  MH_file_a <- paste("__",as.character(j), "_MH_lamb2_",as.character(cur_var), sep = '')
  mhdata <- fromJSON(file = MH_file_a)
  setwd("/Users/Isaac_Zhang/Research/MCMC/simulation_seed/ECOLI_lag_inner_result/ECOLI_2step_gbs")
  GBS_file_a <- paste("__",as.character(j), "_GBS_lamb2_", sep = '')
  gbsdata <- fromJSON(file = GBS_file_a)
  #bmh <- fromJSON(file = "__2_burnin_MH_alpha_0")
  #length(bmh)
  id <- (1:10000) * 1 
  print( ks.test(gbsdata[id], mhdata[id]) )
  
  id2 <- (1:1000) * 10
  print( ks.test(gbsdata[id2], mhdata[id2]) )
  
  
  data <- data.frame(alpha = c(gbsdata, mhdata),method = c(rep("GBS",length(gbsdata)), rep("MH",length(gbsdata)))) 
  
  p1 <- ggplot(data, aes(x=alpha, col = method, linetype = method))#+ coord_fixed(ratio = r1)
  p1 <- p1 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))
  p1 <- p1 + geom_line(stat= "density", alpha = .75,  size = 2)+ theme(legend.position="none")+
    labs(x = expression(alpha), y = expression( paste("P(", lambda, "2|", X, ")" )) )+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  setwd("/Users/Isaac_Zhang/Research/MCMC/revision/New_figures/ecoli_ks/")
  p1
  show(p1)
  filename1 <- paste(name, "l2hist", as.character(j), cur_var, as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename1, width = 6, height = 6)
  }
  
  
  # trace plot
  data_gbs = data[data$method == 'GBS', ]
  maxx = dim(data_gbs)[1]
  maxy = max(data_gbs$alpha)
  miny = min(data_gbs$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p2 <- ggplot(data_gbs, aes((1: length(alpha)), alpha))
  p2 <- p2 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p2 <- p2 + geom_point(alpha=0.6, col = "firebrick1") +  labs(x = "Iteration") +
    labs(y = paste(expression(lambda), "2"))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p2
  show(p2)
  filename2GBS <- paste(name, "l2traceGBS", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename2GBS, width = 6, height = 6)
  }
  
  data_mh = data[data$method == 'MH', ]
  maxx = dim(data_mh)[1]
  maxy = max(data_mh$alpha)
  miny = min(data_mh$alpha)
  ratio.values <- (maxx)/(maxy - miny)
  
  p3 <- ggplot(data_mh, aes((1: length(alpha)), alpha))
  p3 <- p3 + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  p3 <- p3 + geom_point(alpha=0.6, col = "skyblue3") +  labs(x = "Iteration") +
    labs(y = paste(expression(lambda), "2"))+theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) + theme(legend.position="none")
  p3
  
  show(p3)
  
  filenameMH <- paste(name, "l2traceMH", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filenameMH, width = 6, height = 6)
  }
  
  
  bacf <- acf(gbsdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  q
  show(q)
  filename3 <- paste(name, "l2gbsacf", as.character(j), cur_var,as.character(d), ".pdf", sep = "_")
  if(save){
    ggsave(filename3, width = 6, height = 6)
  }
  bacf <- acf(mhdata, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  maxx = dim(bacfdf)[1]
  maxy = max(bacfdf$acf)
  ratio.values <- (maxx)/(maxy)
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw()+ theme(axis.text=element_text(size=40), axis.title=element_text(size=40))+ coord_fixed(ratio = ratio.values)
  
  show(q)
  filename4 <- paste(name, "l2mhacf", as.character(j), cur_var,as.character(d) ,".pdf", sep = "_")
  if(save){
    ggsave(filename4, width = 6, height = 6)
  }
}







f_a(31, 3, 'ecoli', TRUE)

f_b(31, 3, 'ecoli', TRUE)


f_l1(31, 3, 'ecoli', TRUE)
f_l2(31, 3, 'ecoli', TRUE)

