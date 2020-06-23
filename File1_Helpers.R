###########################################
# File1_Helpers.R                         #
# Helper functions for numerical examples #
###########################################

#####################################
# Part 1: helpers for both examples #
#####################################

k.prime <- function(n, r1, r, s, alpha, c1, c2){
    # Computes the optimal number of databases k.prime in the n-th cycle, given 
    # the per-cycle simulation investment budget r1, delivery budget r, and 
    # problem dimension s.
    # Parameters s, alpha, c1, and c2 satisfy the following relationships
    # 1-rho^2 <= c1 |x-x'|^alpha, <= c1 d^alpha <= ck^{-alpha/s}
    
    c = c1*(c2^alpha)*(r-2)/(r-3)
    
    # k.prime according to the Optimal Solutio in the paper
    return((n*r1)^(s/(alpha + s)) * ((c*alpha)/(r*s))^(s / (alpha+s)))
}

visit.schedule <- function(N, r1, r, s, alpha, c1, c2){
    # compute the visiting schedule and quasi-CV candidates for a total of N 
    # investment cycles, given the per-cycle simulation investment budget r1, 
    # delivery budget r, and problem dimension s.Parameters s, alpha, c1, and c2
    # satisfy the following relationships
    # 1-rho^2 <= c1 |x-x'|^alpha, <= c1 d^alpha <= ck^{-alpha/s}
    
    if(N<= 1){
        stop("Argument N should be a positive number.")
    }
    N = as.integer(N) # coerce N to be an integer

    # Placeholders for results
    kn.prime = rep(0, N)    # optimal number of DBs in cycle n
    sn = rep(0, N)          # visit schedule in each of the N cycles
    kn.hat = rep(0, N)      # number of initiated DB in each of the N cycles
    mn = rep(0, N)          # QCV size requirement in each of the N cycles
    kn = rep(0, N)          # number of QCVs in each of the N cycles
    init.ind = rep(F,N)     # initiation indicators
    
    # no. of DBs needed in N cycles according to k.prime
    NDB = ceiling(k.prime(N, r1, r, s, alpha, c1, c2))
    Rk = rep(0, NDB) # database sizes, needed to determine QCVs
    
    # main loop
    for(n in 1:N){
        if(n == 1){
            # first cycle, initiate the first database
            init.ind[n] = T
            kn.prime[n] = k.prime(n, r1, r, s, alpha, c1, c2)
            sn[n] = 1
            Rk[sn[n]] = r1
            kn.hat[n] = 1
            mn[n] = 0
            kn[n] = 1
        } else {
            # all other cycles
            kn.prime[n] = k.prime(n, r1, r, s, alpha, c1, c2)
            if(kn.hat[n-1] < kn.prime[n]){
                # database initiation
                init.ind[n] = T
                sn[n] = kn.hat[n-1] + 1
                Rk[sn[n]] = r1
                kn.hat[n] = kn.hat[n-1] + 1
            } else {
                # database augmentation
                init.ind[n] = F
                kn.hat[n] = kn.hat[n-1]
                ind = which.min(Rk[1:kn.hat[n]])
                sn[n] = ind
                Rk[ind] = Rk[ind] + r1
            }
            
            mn[n] = Rk[1] - r1
            kn[n] = sum(Rk[1:kn.hat[n]] >= mn[n])
        }
    }
    
    return(list(NDB = NDB, Schedule = sn, C.QCV = kn, Init.Ind = init.ind))
}

#######################################
# Part 2: helpers for uniform example #
#######################################

Uniform.est.c1.alpha <- function(NDsn, xmin, xmax, gamma, r){
    # Estimate c1 and alpha via linear regression based on NDsn design points.
    # Uniform r.v. simulation model with parameters xmin, xmax, & threshold 
    # value gamma. Each design point has r replications (CRNs).
    
    # location of design points, QMC points by default
    QMC.est = halton(NDsn, dim = 1, init = TRUE, normal = FALSE, usetime = FALSE)
    DB.Points.est = xmin + (xmax-xmin)*QMC.est
    
    # common random numbers
    CRN = runif(r)
    
    # placeholders and temporary varables
    log.dist = rep(0, NDsn*(NDsn-1)/2)
    log.rho2 = rep(0, NDsn*(NDsn-1)/2)
    ind = 1
    
    for(i in 1:(NDsn-1)){
        x.cur = DB.Points.est[i]
        est.cur = x.cur*CRN > gamma
        
        for(j in (i+1):NDsn){
            x.new = DB.Points.est[j]
            est.new = x.new*CRN > gamma
            
            if(all(est.cur==est.new)){
                cor.cur.new = 1
            } else {
                cor.cur.new = cor(est.cur, est.new)
            }
            log.rho2[ind] = log(1-cor.cur.new^2)
            log.dist[ind] = log(abs(x.cur - x.new))
            ind = ind + 1
        }
    }
    
    # remove Inf entries
    log.dist = log.dist[!is.infinite(log.rho2)]
    log.rho2 = log.rho2[!is.infinite(log.rho2)]
    
    xbar <- mean(log.dist)
    ybar <- mean(log.rho2)
    n <- length(log.dist)
    
    SXY <- sum(log.dist*log.rho2) - n*xbar*ybar
    SXX <- sum(log.dist*log.dist) - n*xbar^2
    
    alpha <- SXY/SXX
    ln.c1 <- ybar - alpha*xbar
    c1 <- exp(ln.c1)
    
    return(list(c1 = c1, alpha = alpha))
}

Uniform.True <- function(gamma, x.pred){
    # True value of the expected value to be estimated by Monte Carlo
    return(1 - gamma/x.pred)
}

Uniform.SMC <- function(gamma, x.pred, NSim){
    # standard MC of uniform r.v. simulation given prediction points x.pred and 
    # threshold value gamma using NSim replications
    
    NPred = length(x.pred)
    results = matrix(runif(NSim*NPred), NSim, NPred)
    results = sweep(results, MARGIN = 2, x.pred, '*') > gamma
    
    mean.est = .colMeans(results, NSim, NPred)
    
    return(mean.est)
}

Uniform.CV.Classic <- function(xmin, xmax, gamma, x.pred, beta.est = FALSE,
                               s, alpha, c1, c2, r, r1){
    # Classic CV (true mean) of uniform r.v. simulation given model parameters 
    # xmin, xmax, and threshold value gamma.
    # Prediction points are specified in x.pred.
    # User decides if CV parameter beta is estimated or calculated (default)
    
    # Parameters N, s, alpha, c1, c2, r, and r1 are used to compute databases 
    # for comparison with GDB)
    
    # compute visiting schedule and nearest neighbor DBs to be used as CV
    N = length(x.pred)
    
    Visit.Schedule = visit.schedule(N, r1, r, s, alpha, c1, c2)
    NDB = Visit.Schedule$NDB    # no. of DBs
    C.QCV = Visit.Schedule$C.QCV        # quasi-CV candidates in each cycle
    # the actual visiting schedule does not matter in classical CV because 
    # all DBs are assumed an accurate mean (no visit required at all)
    
    QMC = halton(NDB, dim = 1, init = TRUE, normal = FALSE, usetime = FALSE)
    DBs = xmin + (xmax - xmin) * QMC    # location of DBs
    
    if(beta.est){
        CRN <- runif(r)             # common random numbers (CRNs)
        CV.CRN = mat.or.vec(r, NDB) # CRNs in each DB for implementation convinience
        for(id in 1:NDB){
            CV.CRN[, id] <- DBs[id] * CRN > gamma
        }
    }
    
    # placeholder of result
    CV.mean = rep(0, N)
    
    # Classical CV in each cycle
    for(cycleID in 1:N){
        
        # current prediction point and the closest DB point
        x.cur = x.pred[cycleID]
        IDX = which.min(abs(x.cur - DBs[1:C.QCV[cycleID]]))
        x.DB = DBs[IDX]
        CV.True = 1 - gamma / x.DB  # true mu at each DB
        CV.cur = x.cur * CRN > gamma
        
        # Calculate beta
        if(beta.est){   
            # estimate beta
            CV.DB = CV.CRN[, IDX]
            var.DB = var(CV.DB)
            if(var.DB==0){
                stop("DB CRN variance is zero.")
            }
            beta = cov(CV.cur, CV.DB) / var.DB
        } else {
            # calculate true beta in classical CV
            if(x.cur < x.DB){
                beta = x.DB*(x.cur-gamma)/(x.cur*(x.DB-gamma))
            } else {
                beta = x.DB/x.cur
            }
        }
        
        # CV estimator
        CV.mean[cycleID] <- mean(CV.cur) - beta*(mean(CV.DB) - CV.True)
    }
    
    return(CV.mean)
}

Uniform.CV.GDB <- function(xmin, xmax, gamma, x.pred, beta.est = FALSE,
                           s, alpha, c1, c2, r, r1){
    # Green CV (DB mean) of uniform r.v. simulation given model parameters 
    # xmin, xmax, and threshold value gamma.
    # Prediction points are specified in x.pred.
    # User decides if CV parameter beta is estimated or calculated (default)
    
    # Experiment design parameters N, s, alpha, c1, c2, r, and r1 are used
    # to compute databases and their visiting schedule (could be approx.)
    
    # compute visiting schedule and nearest neighbor DBs to be used as CV
    N = length(x.pred)
    Visit.Schedule = visit.schedule(N, r1, r, s, alpha, c1, c2)
    NDB = Visit.Schedule$NDB    # no. of DBs
    
    QMC = halton(NDB, dim = 1, init = TRUE, normal = FALSE, usetime = FALSE)
    DBs = xmin + (xmax-xmin)*QMC        # location of DBs
    C.QCV = Visit.Schedule$C.QCV        # quasi-CV candidates in each cycle
    Schedule = Visit.Schedule$Schedule  # DB visited in each cycle
    Init.Ind = Visit.Schedule$Init.Ind  # Initiation indicators
    
    # placeholder of result
    CV.mean = rep(0, N) # GDB estimates    
    
    CRN = runif(r)              # common random numbers
    CV.CRN = mat.or.vec(r, NDB) # r CRNs, 1st story in DB
    CV.GDB = rep(0, NDB)        # DB mean estimate, 2nd story in DB
    Rk = rep(0, NDB)            # DB sizes
    
    # GDB procedure in each cycle
    for(cycleID in 1:N){
        
        # ===== Investment phase
        # DB visited at current cycle
        sn <- Schedule[cycleID]
        DB.cur <- DBs[sn]
        
        if(Init.Ind[cycleID]){
            # invest CRNs in the first story of visited DB for initiations
            CV.CRN[, sn] <- DB.cur * CRN > gamma
        }
        New.Mean <- mean(DB.cur * runif(r1) > gamma)
        CV.GDB[sn] <- CV.GDB[sn] * Rk[sn]/(Rk[sn] + r1) + New.Mean * r1/(Rk[sn] + r1)
        Rk[sn] <- Rk[sn] + r1
        
        # ===== Observe xn
        xn = x.pred[cycleID]
        
        # ===== Delivery phase
        CV.cur = xn * CRN > gamma   # run replications at xn
        
        # closest candidate DB as the quasi-CV
        IDX = which.min(abs(xn - DBs[1:C.QCV[cycleID]]))
        x.DB = DBs[IDX]
        CV.DB = CV.CRN[, IDX]
        
        # Calculate beta
        if(beta.est){   
            # estimate beta
            var.DB = var(CV.DB)
            if(var.DB == 0){
                stop("DB CRN variance is zero.")
            }
            beta = (cov(CV.cur, CV.DB) / var.DB)
        } else {
            # calculate true beta in classical CV
            if(xn < x.DB){
                beta = x.DB*(xn-gamma)/(xn*(x.DB-gamma))
            } else {
                beta = x.DB/xn
            }
        }
        
        DBCV.factor = Rk[IDX]/(Rk[IDX]+r)   # beta adj. factor for using DB as CV
        beta = beta*DBCV.factor
        
        # CV estimators
        CV.mean[cycleID] <- mean(CV.cur) - beta*(mean(CV.DB) - CV.GDB[IDX])
    }
    
    return(CV.mean)
}

###########################################
# Part 3: helpers for reinsurance example #
###########################################

CatBond.est.c1.alpha <- function(NDsn, LambdaRange, ThetaRange, 
                                     K, p, r){
    # Estimate c1 and alpha using linear regression based on NDsn DBs
    # Function if model specific: used for reinsurance example with model 
    # parameters LambdaRange, ThetaRange, Trigger, and Recovery
    
    # location of DBs
    DBs = halton(NDsn, dim = 2, init = TRUE, normal = FALSE, usetime = FALSE)
    DBs[, 1] = LambdaRange[1] + (LambdaRange[2] - LambdaRange[1])*DBs[, 1]  
    DBs[, 2] = ThetaRange[1] + (ThetaRange[2] - ThetaRange[1])*DBs[, 2]
    
    # common random numbers
    CRN.poi = runif(r)  # CDF of common random Poisson r.v.
    CRN.gam = runif(r)  # CDF of common random Gamma r.v.
    
    # placeholders and temporary varables
    log.dist = rep(0, NDsn*(NDsn-1)/2)
    log.rho2 = rep(0, NDsn*(NDsn-1)/2)
    ind = 1
    
    for(i in 1:(NDsn - 1)){
        xn = DBs[i, ]
        
        CRN.poi.cur = qpois(CRN.poi, xn[1])
        est.cur = p + (1-p) * (qgamma(CRN.gam, shape = CRN.poi.cur, 
                                      scale = xn[2]) <= K)
        
        for(j in (i+1):NDsn){
            x.new = DBs[j, ]
            CRN.poi.new = qpois(CRN.poi, x.new[1]) 
            est.new = p + (1-p) * (qgamma(CRN.gam, shape = CRN.poi.new, 
                                          scale = x.new[2]) <= K)
            
            if(all(est.cur == est.new)){
                cor.cur.new = 1
            } else {
                cor.cur.new = cor(est.cur, est.new)
            }
            
            log.rho2[ind] = log(1-cor.cur.new^2)
            log.dist[ind] = log(sqrt(sum((xn - x.new)^2)))
            ind = ind + 1
        }
    }
    
    # remove Inf entries
    log.dist = log.dist[!is.infinite(log.rho2)]
    log.rho2 = log.rho2[!is.infinite(log.rho2)]
    
    model <- lm(log.rho2~log.dist)
    c1 <- exp(model$coefficients[1])
    alpha <- model$coefficients[2]
    
    return(list(c1 = c1, alpha = alpha))
}

CatBond.True <- function(lambda, theta, K, p){
    # compute true CAT bond price in the reinsurance example with lambda and 
    # theta being the means of the involved Poisson and Exponential random 
    # variables. K and p are the trigger value and the recovery rate.
    
    # number of terms is approximating true price
    nTerms = max(100, 2*max(lambda))
    
    # unify inputs as matrices
    lambda = as.matrix(lambda)
    theta = as.matrix(theta)
    
    # Poisson pmf and Gamma cdf
    Pj = apply(lambda, MARGIN = 1, function(x) dpois(0:nTerms, x))
    Fj = apply(theta, MARGIN = 1, function(x) pgamma(K, 0:nTerms, scale = x))
    price = p + (1 - p) * colSums(Pj*Fj)
    return(price)
}

CatBond.SMC <- function(lambda, theta, K, p, NSim){
    # standard MC of reinsurance simulation given trigger value K, recovery 
    # rate p, and prediction points x.pred using NSim replications
    
    # unify inputs as matrices
    lambda = as.matrix(lambda)
    theta = as.matrix(theta)
    
    NPt = length(lambda)   # no. of inputs
    results = rep(0, NPt)   # placeholder of results
    for(pt in 1:NPt){
        
        lambda.cur = lambda[pt]
        theta.cur = theta[pt]
        
        temp = rpois(NSim, lambda.cur)  # pred frequencies
        temp = rgamma(NSim, shape = temp, scale = theta.cur)  # pred losses
        temp = p + (1-p) * (temp <= K)
        results[pt] = mean(temp)
    }
    
    return(results)
}

CatBond.CV.Classic <- function(PredPt, K, p, LambdaRange, ThetaRange,
                               s, alpha, c1, c2, r, r1){
    # Classic CV (true mean) of CatBond simulation given model parameters 
    # LambdaRange, ThetaRange, trigger K, and recovery p.
    # Prediction points are specified in PredPt, a 2-columns dataframe
    # Experiment design parameters s, alpha, c1, c2, r, and r1 are used
    # to compute databases (for comparison with GreenDB)
    
    # compute visiting schedule and nearest neighbor DBs to be used as CV
    NCycles = dim(PredPt)[1]
    Visit.Schedule = visit.schedule(NCycles, r1, r, s, alpha, c1, c2)
    NDB = Visit.Schedule$NDB        # no. of DBs
    C.QCV = Visit.Schedule$C.QCV    # quasi-CV candidates
    
    # compute location of DBs
    DBs = halton(NDB, dim = 2, init = TRUE, normal = FALSE, usetime = FALSE)
    DBs[, 1] = LambdaRange[1] + (LambdaRange[2] - LambdaRange[1]) * DBs[, 1]  
    DBs[, 2] = ThetaRange[1] + (ThetaRange[2] - ThetaRange[1]) * DBs[, 2]   
    
    # placeholder of result
    CV.mean = rep(0, NCycles)
    
    # common random numbers
    CRN.poi = runif(r)  # CDF of common random Poisson r.v.
    CRN.gam = runif(r)  # CDF of common random Gamma r.v.
    
    # For implementation convenience, fill in all the CRNs in each DB first
    CV.CRN = mat.or.vec(r, NDB) # CRNs in each DB for implementation convinience
    for(id in 1:NDB){
        x.DB = DBs[id, ]
        CRN.poi.DB = qpois(CRN.poi, x.DB[1])
        CV.CRN[, id] = p + (1-p) * (qgamma(CRN.gam, shape = CRN.poi.DB, 
                                           scale = x.DB[2]) <= K)
    }
    
    # Classical CV
    for(cycleID in 1:NCycles){
        
        # current prediction point and nearst DB point
        xn = PredPt[cycleID, ]
        DBs.cur = DBs[1:C.QCV[cycleID], , drop = F]
        IDX = which.min(apply(DBs.cur, 1, function(pt) sum((pt - xn)^2)))
        x.DB = DBs[IDX, ]
        
        CV.True = CatBond(x.DB[1], x.DB[2], K, p)  # true mu at chosen DB
        
        CRN.poi.cur = qpois(CRN.poi, xn[1])
        CV.cur = p + (1-p) * (qgamma(CRN.gam, shape = CRN.poi.cur, 
                                     scale = xn[2]) <= K)
        CV.DB = CV.CRN[, IDX]
        
        # estimate beta
        var.DB = var(CV.DB)
        if(var.DB == 0){
            stop("DB CRN variance is zero.")
        }
        beta = cov(CV.cur, CV.DB) / var.DB
        
        # CV estimator
        CV.mean[cycleID] <- mean(CV.cur) - beta*(mean(CV.DB) - CV.True)
    }
    
    return(CV.mean)
}

CatBond.CV.GDB <- function(PredPt, K, p, LambdaRange, ThetaRange,
                           s, alpha, c1, c2, r, r1){
    # GDB CV (DB mean) of CatBond simulation given model parameters 
    # LambdaRange, ThetaRange, trigger K, and recovery p.
    # Prediction points are specified in PredPt, a 2-columns dataframe
    # Experiment design parameters s, alpha, c1, c2, r, and r1 are used
    # to compute databases
    
    # compute visiting schedule and nearest neighbor DBs to be used as CV
    NCycles = dim(PredPt)[1]
    Visit.Schedule = visit.schedule(NCycles, r1, r, s, alpha, c1, c2)
    NDB = Visit.Schedule$NDB        # no. of DBs
    C.QCV = Visit.Schedule$C.QCV        # quasi-CV candidates in each cycle
    Schedule = Visit.Schedule$Schedule  # DB visited in each cycle
    Init.Ind = Visit.Schedule$Init.Ind  # Initiation indicators
    
    # compute location of DBs
    DBs = halton(NDB, dim = 2, init = TRUE, normal = FALSE, usetime = FALSE)
    DBs[, 1] = LambdaRange[1] + (LambdaRange[2] - LambdaRange[1]) * DBs[, 1]  
    DBs[, 2] = ThetaRange[1] + (ThetaRange[2] - ThetaRange[1]) * DBs[, 2]   
    
    # placeholder of result
    CV.mean = rep(0, NCycles)
    
    # common random numbers
    CRN.poi = runif(r)  # CDF of common random Poisson r.v.
    CRN.gam = runif(r)  # CDF of common random Gamma r.v.
    
    CV.CRN = mat.or.vec(r, NDB)    # r CRNs, 1st story in DB
    CV.GDB = rep(0, NDB)      # DB estimate, 2nd story in DB
    Rk = rep(0, NDB)            # DB sizes
    
    # GDB procedure
    for(cycleID in 1:NCycles){
        
        # ===== Investment phase
        # DB visited at current cycle
        sn <- Schedule[cycleID]
        x.DB <- DBs[sn, ]
        
        if(Init.Ind[cycleID]){
            # invest CRNs in the first story of visited DB for initiations
            CRN.poi.DB = qpois(CRN.poi, x.DB[1])
            CV.CRN[, sn] = p + (1-p) * (qgamma(CRN.gam, shape = CRN.poi.DB, 
                                               scale = x.DB[2]) <= K)
        }
        New.Mean <- CatBond.SMC(x.DB[1], x.DB[2], K, p, r1)
        CV.GDB[sn] <- CV.GDB[sn] * Rk[sn]/(Rk[sn] + r1) + New.Mean * r1/(Rk[sn] + r1)
        Rk[sn] <- Rk[sn] + r1
        
        # ===== Observe xn
        xn <- PredPt[cycleID, ]
        DBs.cur <- DBs[1:C.QCV[cycleID], , drop = F]
        IDX = which.min(apply(DBs.cur, 1, function(pt) sum((pt - xn)^2)))
        x.DB = DBs[IDX, ]
        
        CRN.poi.cur = qpois(CRN.poi, xn[1])
        CV.cur = p + (1-p) * (qgamma(CRN.gam, shape = CRN.poi.cur, 
                                     scale = xn[2]) <= K)
        
        CRN.poi.DB = qpois(CRN.poi, x.DB[1])
        CV.DB = CV.CRN[, IDX]
        
        # estimate beta
        var.DB = var(CV.DB)
        if(var.DB == 0){
            stop("DB CRN variance is zero.")
        }
        beta = (cov(CV.cur, CV.DB) / var.DB)*Rk[IDX]/(Rk[IDX]+r)

        # CV estimators
        CV.mean[cycleID] <- mean(CV.cur) - beta*(mean(CV.DB) - CV.GDB[IDX])
        
    }
    
    return(CV.mean)
}













