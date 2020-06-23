# ===== Workspace preparation
rm(list = ls()) # clear all variables
set.seed(0) # fix random seed for debugging and reproduction of results

library(randtoolbox)    # QMC sequence generation
library(pracma)         # tic-toc (for now)
source("File1_Helpers.R")
source("File1_Utilities.R")

# ===== Model parameters
# bounds for Poisson mean (lambda) & exponential mean (theta)
LambdaRange <- c(2, 4)
ThetaRange <- c(4, 6)

Trigger <- 16     # triggering point of CAT bond
Recovery <- 0.5    # recovery rate in CAT bond

s <- 2   # problem dimension
c2 <- 2*3  # formula for c2, 2D Halton w/ primes 2 and 3, also scaling factor 2

r <- 200     # no. of replications for per-period estimation
r1 <- 1000    # no. of replications for per-period DB investment

# ===== Experiment design parameters
NRep <- 1e4                      # no. of (macro) replications
NPeriods <- 1e2                      # no. of periods
NDB.est <- 25                        # no. of DBs for estimating c1 and alpha

# ===== Beginning of experiment

# ***** Placeholders of results
Est.True <- mat.or.vec(NRep, NPeriods)
Est.SMC <- mat.or.vec(NRep, NPeriods)
Est.SMC.r1 <- mat.or.vec(NRep, NPeriods)
Est.GDB.Est <- mat.or.vec(NRep, NPeriods)

for(rep in 1:NRep){
    tic()
    
    # prediction points for the current macro replication
    x.cur <- matrix(runif(2*NPeriods), NPeriods, 2)
    lambda <- LambdaRange[1] + (LambdaRange[2] - LambdaRange[1]) * x.cur[, 1]
    theta <- ThetaRange[1] + (ThetaRange[2] - ThetaRange[1]) * x.cur[, 2]
    
    # true reinsurance liability and SMC estimates
    Est.True[rep, ] <- CatBond.True(lambda, theta, Trigger, Recovery)  
    Est.SMC[rep, ] <- CatBond.SMC(lambda, theta, Trigger, Recovery, r) 
    Est.SMC.r1[rep, ] <- CatBond.SMC(lambda, theta, Trigger, Recovery, r1) 
    
    c1.alpha <- CatBond.est.c1.alpha(NDB.est, LambdaRange, ThetaRange, 
                                 Trigger, Recovery, r)
    c1 <- c1.alpha$c1
    alpha <- c1.alpha$alpha
    
    Est.GDB.Est[rep, ] <- 
        CatBond.CV.GDB(cbind(lambda, theta), Trigger, Recovery, 
                           LambdaRange, ThetaRange,
                           s, alpha, c1, c2, r, r1)
    
    print(sprintf("This is replication %d", rep))
    toc()
}

save.image(file = "Reinsurance_Ind_100Periods.RData")

# calculate MSEs for various estimators
MSE.SMC <- (Est.True - Est.SMC)^2
MSE.SMC.r1 <- (Est.True - Est.SMC.r1)^2
MSE.GDB.Est <- (Est.True - Est.GDB.Est)^2

# mean MSEs and their standard deviations
SMC <- colMeansSds(MSE.SMC)
SMC.r1 <- colMeansSds(MSE.SMC.r1)
GDB.Est <- colMeansSds(MSE.GDB.Est)

# mean MSEs
NPeriods <- length(SMC$ColMean)
dat.SMC <- data.frame(Estimators = "SMC.r", period = 1:NPeriods, Val = mean(SMC$ColMean))
dat.SMC.r1 <- data.frame(Estimators = "SMC.r1", period = 1:NPeriods, Val = mean(SMC.r1$ColMean))
dat.GDB.Est <- data.frame(Estimators = "GDB.Est", period = 1:NPeriods, Val = GDB.Est$ColMean)
slope <- mean(alpha)/(mean(alpha) + s)
dat.ref <- data.frame(Estimators = "Slope.Ref", period = 1:NPeriods, Val = mean(SMC$ColMean)/(1:NPeriods)^slope)

dat <- rbind(dat.SMC, dat.SMC.r1, dat.GDB.Est, dat.ref)

library(ggplot2)
library(scales)
ggplot(data = dat, aes(x = period, y = Val, colour = Estimators, linetype = Estimators)) +
    scale_colour_manual(values = c("#7b3294","#c2a5cf","#008837","black")) + 
    geom_line(size = 2) +
    scale_linetype_manual(values = c("solid","solid","longdash","solid")) +
    scale_x_continuous(name = "Periods (No. Repeated Experiments)", trans = "log10", 
                       breaks = c(1,10,100),
                       labels = c("1", "10", "100")) +
    scale_y_continuous(name ="Mean Squared Error", trans = "log10", 
                       breaks = c(2e-4, 7e-5, 4e-5),
                       limits = c(4e-5, 3e-4)) +theme_bw()+
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16,face="bold"),
          axis.line = element_line(size = 1, colour = "black"),
          legend.justification = c(-0.1,-0.1), legend.position = c(0,0),
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black"),
          legend.direction = "horizontal") + 
    guides(colour = guide_legend(keywidth = 5, keyheight = 1.5))

mean(SMC$ColMean)
mean(SMC.r1$ColMean)
GDB.Est$ColMean[c(1,10,100)]











