# ===== Workspace preparation

rm(list = ls()) # clear all variables
set.seed(1) # fix random seed for debugging and reproduction of results

library(randtoolbox)    # QMC sequence generation
library(pracma)         # tic-toc (for now)

source('File1_Helpers.R')      
source("File1_Utilities.R")

# ===== Model parameters
xmin <- 10 
xmax <- 12   # range of upperbound of U(0,x) r.v.'s
gamma <- 8   # threshold valud to calculate P(U>gamm)

s <- 1   # problem dimension
c1 <- 1/(xmax - gamma)   # 1-rho^2 <= c1|x-x'|^alpha
c2 <- 2*(xmax-xmin)      # this is only true for 1D Halton sequence
alpha <- 1               # this is only true for 1D Halton sequence
c <- c1*(c2^alpha)

r <- 200     # no. of replications for per-period estimation
r1 <- 1000    # no. of replications for per-period DB investment

# ===== Experiment design parameters
NRep <- 1e4                         # no. of (macro) replications
NCycles <- 1e3                      # no. of periods
NDB.Est <- 10                       # no. of DBs for estimating c1 and alpha

# ===== Beginning of experiment

# ***** Placeholders of results
Est.True <- mat.or.vec(NRep, NCycles)
Est.SMC <- mat.or.vec(NRep, NCycles)
Est.SMC.r1 <- mat.or.vec(NRep, NCycles)
Est.GDB <- mat.or.vec(NRep, NCycles)
Est.GDB.Est <- mat.or.vec(NRep, NCycles)

c1.Est.vec <- mat.or.vec(NRep, 1)
alpha.Est.vec <- mat.or.vec(NRep, 1)

for(rep in 1:NRep){
    tic()
    
    # prediction points for the current macro replication
    x.cur <- xmin + (xmax-xmin)*runif(NCycles)
    
    Est.True[rep, ] <- Uniform.True(gamma, x.cur)  # true mu
    Est.SMC[rep, ] <- Uniform.SMC(gamma, x.cur, r) # SMC with r reps.
    Est.SMC.r1[rep, ] <- Uniform.SMC(gamma, x.cur, r1) # SMC with r1 reps.
    
    Est.GDB[rep, ] <-
        Uniform.CV.GDB(xmin, xmax, gamma, x.cur, beta.est = FALSE,
                         s, alpha, c1, c2, r, r1)
    
    c1.alpha <- Uniform.est.c1.alpha(10, xmin, xmax, gamma, r)
    c1.Est <- c1.alpha$c1
    alpha.Est <- c1.alpha$alpha
    Est.GDB.Est[rep, ] <-
        Uniform.CV.GDB(xmin, xmax, gamma, x.cur, beta.est = TRUE,
                         s, alpha.Est, c1.Est, c2, r, r1)
    
    c1.Est.vec[rep] <- c1.Est
    alpha.Est.vec[rep] <- alpha.Est
    print(sprintf("This is replication %d", rep))
    toc()
}

save.image(file = "Uniform_Ind.RData")

# calculate MSEs for various estimators
MSE.SMC <- (Est.True - Est.SMC)^2
MSE.SMC.r1 <- (Est.True - Est.SMC.r1)^2
MSE.GDB <- (Est.True - Est.GDB)^2
MSE.GDB.Est <- (Est.True - Est.GDB.Est)^2

SMC <- colMeansSds(MSE.SMC)
SMC.r1 <- colMeansSds(MSE.SMC.r1)
GDB <- colMeansSds(MSE.GDB)
GDB.Est <- colMeansSds(MSE.GDB.Est)

# mean MSEs
NCycles <- length(SMC$ColMean)
dat.SMC <- data.frame(Estimators = "SMC.r", period = 1:NCycles, Val = mean(SMC$ColMean))
dat.SMC.r1 <- data.frame(Estimators = "SMC.r1", period = 1:NCycles, Val = mean(SMC.r1$ColMean))
dat.GDB <- data.frame(Estimators = "GDB", period = 1:NCycles, Val = GDB$ColMean)
dat.GDB.Est <- data.frame(Estimators = "GDB.Est", period = 1:NCycles, Val = GDB.Est$ColMean)
dat.ref <- data.frame(Estimators = "Slope.Ref", period = 1:NCycles, Val = mean(SMC$ColMean)/sqrt(1:NCycles))

dat <- rbind(dat.SMC, dat.SMC.r1, dat.GDB, dat.GDB.Est, dat.ref)

library(ggplot2)
library(scales)
dev.new()
ggplot(data = dat, aes(x = period, y = Val, colour = Estimators, linetype = Estimators)) +
    scale_colour_manual(values = c("#7b3294","#c2a5cf","#008837","#a6dba0","black")) + 
    geom_line(size = 2) +
    scale_linetype_manual(values = c("solid","solid","dotted","longdash","solid")) +
    scale_x_continuous(name = "Number of Investment-Delivery Cycles", trans = "log10", 
                       breaks = c(1,10,100,1000),
                       labels = c("1", "10", "100", "1000")) +
    scale_y_continuous(name ="Mean Squared Error", trans = "log10", 
                       breaks = c(1e-3,1e-4,1e-5,1e-6),
                       labels = trans_format("log10", math_format(10^.x)),
                       limits = c(1e-5, 1.5e-3)) +theme_bw()+
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1), 
          legend.position = c(0,0),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.background = element_rect(colour ="black")) + 
    guides(colour = guide_legend(keywidth = 5, keyheight = 1.5))

dat.c1.alpha <- data.frame(c1 = c1.Est.vec, alpha = alpha.Est.vec)
p1 <- ggplot(dat.c1.alpha, aes(x = c1)) +
    geom_histogram(bins = 80, fill="#7fbf7b", colour="#f7f7f7") +
    geom_vline(aes(xintercept = 0.25, color = "mean"),
               color="#af8dc3", linetype="longdash", size = 2) +
    scale_x_continuous(name = "Value", limits = c(0.07,0.5),
                       breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5)) +
    scale_y_continuous(name ="Count") +
    ggtitle(expression(paste("Histogram of estimated ", "c"[1]))) +
    theme(plot.title = element_text(size = 16, face="bold", hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14,face="bold"),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"))

p2 <- ggplot(dat.c1.alpha, aes(x = alpha)) +
    geom_histogram(bins = 80, fill="#7fbf7b", colour="#f7f7f7") +
    geom_vline(aes(xintercept = 1, color = "mean"),
               color="#af8dc3", linetype="longdash", size = 2) +
    scale_x_continuous(name = "Value", limits = c(0.65, 1.1),
                       breaks = c(0.7, 0.8, 0.9, 1.0, 1.1)) +
    scale_y_continuous(name ="Count") +
    ggtitle(expression(paste("Histogram of estimated ", alpha))) +
    theme(plot.title = element_text(size = 16, face="bold", hjust = 0.5),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 14,face="bold"),
          panel.background = element_rect(fill = "grey95"),
          axis.line = element_line(size = 1, colour = "black"))
dev.new()
multiplot(p1, p2, cols = 2)

mean(SMC$ColMean)
mean(SMC.r1$ColMean)
GDB$ColMean[c(1,10,100,1000)]
GDB.Est$ColMean[c(1,10,100,1000)]
