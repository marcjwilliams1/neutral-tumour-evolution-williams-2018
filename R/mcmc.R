library(coda)
library(adaptMCMC)
library(VGAM)


LLbetabin <- function(data, theta){
  p <- theta[1]
  rho <- theta[2]
  R <- dbetabinom(x = data[,2], size = data[,1], prob = p, rho = rho)
  return(sum(log(R)))
}

LLbin <- function(data, theta){
  p <- theta[1]
  rho <- theta[2]
  R = dbetabinom(x = data[,1], size = data[,2], prob = p)
  R <- R[R>0]
  return(sum(log(R)))
}

prior <- function(theta, plims){
  p1 <- dunif(theta[1],min=plims[1,1],max=plims[1,2])
  p2 <- dunif(theta[2],min=plims[2,1],max=plims[2,2])
  return( p1*p2 )
}

priorbin <- function(theta, plims){
  p1 <- dunif(theta[1],min=plims[1],max=plims[2])
  return( p1  )
}

randbb <- function (num, pars){
  depth <- pars[1]
  p <- pars[2]
  rho <- pars[3]
  return(rbetabinom(num, size = depth, prob = p, rho = rho) / depth)
}

randbin <- function (num, pars){
  depth <- pars[1]
  p <- pars[2]
  return(rbinom(num, size = depth, prob = p) / depth)
}

makeplots <- function(post, data, depth, lims){
  
  VAFs <- randbb(100000, c(depth, median(post[,1]), median(post[,2])))
  dffit <- data.frame(VAF = VAFs, 
                   type = paste("Beta Binomial (rho = ", round(median(post[,2]),4), ")", sep=""))
  
  VAFs <- randbin(100000, c(depth, median(post[,1])))
  dffit2 <- data.frame(VAF = VAFs, 
                      type = "Binomial")
  
  k = normdensity(data, lims)
  
  g1 <- ggplot(data, aes(VAF)) + 
    geom_histogram(aes(k=k,y=k*(..density..)), bins = 100, fill = "azure4",alpha = 0.8) + 
    geom_line(data = dffit, stat="density",aes(VAF, col = type), size = 2.0, adjust = 5) +
    geom_line(data = dffit2, stat="density",aes(VAF, col = type), size = 2.0, adjust = 5, linetype = 2) +
    #geom_density(data = df, aes(VAF), size = 2.0, col = "deepskyblue4") +
    scale_colour_manual(values=c(alpha("deepskyblue4",0.6), alpha("firebrick4",0.6))) +
    theme_bw(base_family = 'Helvetica') + 
    xlab("VAF") +
    ylab("Density") +
    xlim(c(0.0,0.7)) +
    #ylim(c(0,60)) +
    theme(legend.title = element_blank())

  dfpost <- data.frame(rho = post[,2])
  g2 <- ggplot(dfpost, aes(rho)) + 
    geom_histogram(fill = "firebrick4", alpha = 0.8, bins = 100) +
    theme_bw(base_family = 'Helvetica') + 
    xlab("rho") +
    ylab("Frequency")

  dfpost <- data.frame(p = post[,1])
  g3 <- ggplot(dfpost, aes(p)) + 
    geom_histogram(fill = "deepskyblue4", alpha = 0.8, bins = 100) +
    theme_bw(base_family = 'Helvetica') + 
    xlab("Frequency of clonal mutations") +
    ylab("Frequency")

  out <- list()
  out[[1]] <- g1
  out[[2]] <- g2
  out[[3]] <- g3
  return(out)
  
}

normdensity <- function (df, lims){
  d <- density(df$VAF, n = 100, from = 0.0, to = 1.0)
  totd <- sum(d$y)
  d <- data.frame(x = d$x, y = d$y)
  d <- dplyr::filter(d, x > lims[1], x < lims[2])
  return(1/(sum(d$y)/totd))
}

makeplotsbin <- function(post, data, depth, lims){
  
  VAFs <- randbin(100000, c(depth, median(post[,1])))
  dffit <- data.frame(VAF = VAFs, 
                      type = "Binomial")
  k = normdensity(data, lims)
  
  g1 <- ggplot(data, aes(VAF)) + 
    geom_histogram(aes(k=k,y=k*(..density..)), bins = 100, fill = "azure4",alpha = 0.8) + 
    geom_line(data = dffit, stat="density",aes(VAF, col = type), size = 2.0, adjust = 5) +
    #geom_density(data = df, aes(VAF), size = 2.0, col = "deepskyblue4") +
    scale_colour_manual(values=alpha("deepskyblue4",0.6)) +
    theme_bw(base_family = 'Helvetica') + 
    xlab("VAF") +
    ylab("Density") +
    xlim(c(0.0,0.7)) +
    #ylim(c(0,60)) +
    theme(legend.title = element_blank())
  
  
  dfpost <- data.frame(p = post[,1])
  g2 <- ggplot(dfpost, aes(p)) + 
    geom_histogram(fill = "deepskyblue4", alpha = 0.8, bins = 100) +
    theme_bw(base_family = 'Helvetica') + 
    xlab("Frequency of clonal mutations") +
    ylab("Frequency")
  
  out <- list()
  out[[1]] <- g1
  out[[2]] <- g2
  return(out)
  
}