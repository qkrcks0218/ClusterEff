############################################
# Last update : October 09 2021
# Generate Simulation Data
############################################

############################################
# Load Packages
############################################

library(rlearner)
library(grf)
library(lme4)
library(foreach)
library(doParallel)

PATH <- getwd()

registerDoParallel(10) 

foreach (ii=1:420, .combine = rbind)  %dopar% {
    
    setwd(PATH)
    
    
    BATCH <- ii    
    source("../ClusterFtSource.R")
    pos.X <- 3:8
    pos.AX <- 2:8
    CPS.B <- 5
    
    PARA.GRID <- expand.grid(c(500),1:200,c(0,1.5))
    N <- PARA.GRID[BATCH,1]
    New.BATCH <- PARA.GRID[BATCH,2]
    sV <- PARA.GRID[BATCH,3]
    sVType <- if(sV==0){"No"} else {"Strong"}
    
    set.seed( 1000*N+New.BATCH )
    
    M <- round( runif(N, min=2500/N*0.8, max=2500/N*1.2) )
    M.start <- cumsum(M)-M+1
    M.end <- cumsum(M)
    
    
    ATE <- list()
    
    GP <- rep(1:N, M)
    
    V.GP <- rnorm(N)*sV
    V <- rep(V.GP,M)
    
    Cl.X <- cbind( rep(rnorm(N),M), rep(rbinom(N,1,0.7),M) )
    X <- cbind( rnorm(sum(M)), rnorm(sum(M)), rbinom(sum(M),1,0.3), Cl.X)
    
    
    
    PS.function <- function(X,V){
        return( expit(-0.5 + 0.5*X[,1] - as.numeric(X[,2]>1) + 0.5*X[,3] - 0.25*X[,4] + X[,5] + V )  )
    }
    
    PS <- PS.function(X,V)
    A <- rbinom(sum(M),1,PS)
    
    Y <- list()
    true.g <- list()
    true.g.1 <- list()
    true.g.0 <- list()
    
    sU <- list()
    sY <- list()
    U <- list()
    
    sU[[1]] <- 0.5
    sY[[1]] <- 1
    U.GP <- rnorm(N)*sU[[1]]
    U[[1]] <- rep(U.GP,M)
    eps <- rnorm(sum(M))*sY[[1]]
    
    eps <- rep(0,sum(M))
    Cl.eps <- rep(0,N)
    for(jj in 1:N){
        
        pos <- M.start[jj]:M.end[jj]
        TYPE <- runif(1)
        if(TYPE <= 0.25){
            Cl.eps[jj] <- rnorm(1)*0.5
            eps[pos] <- rnorm(M[jj]) + rnorm(1)*0.5
        } else if (TYPE <= 0.5){
            Cl.eps[jj] <- rt(1,5)*0.5
            eps[pos] <- rnorm(M[jj]) + rt(1,5)*0.5
        } else if (TYPE <= 0.75){
            Cl.eps[jj] <- sample(c(-1,1),1)*rexp(1)*0.5
            eps[pos] <- rnorm(M[jj]) + sample(c(-1,1),1)*rexp(1)*0.5
        } else {
            Cl.eps[jj] <- (runif(1)-0.5)*0.5
            eps[pos] <- rnorm(M[jj]) + (runif(1)-0.5)*0.5
        }
        
    }
    
    Y[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*A + 2*X[,1] - X[,4]^2 + X[,2]*X[,5] + eps
    true.g[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*A + 2*X[,1] - X[,4]^2 + X[,2]*X[,5]
    true.g.1[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*1 + 2*X[,1] - X[,4]^2 + X[,2]*X[,5]
    true.g.0[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*0 + 2*X[,1] - X[,4]^2 + X[,2]*X[,5]
    
    ATE[[1]] <- 4
    
    
    sU[[2]] <- 1.5
    sY[[2]] <- 1
    U.GP <- rnorm(N)*sU[[2]]
    U[[2]] <- rep(U.GP,M)
    eps <- rnorm(sum(M))*sY[[2]]
    
    Y[[2]] <- true.g[[2]] <- true.g.1[[2]] <- true.g.0[[2]] <- rep(0,sum(M))
    Y[[3]] <- true.g[[3]] <- true.g.1[[3]] <- true.g.0[[3]] <- rep(0,sum(M))
    Y[[4]] <- true.g[[4]] <- true.g.1[[4]] <- true.g.0[[4]] <- rep(0,sum(M))
    
    Data <- data.frame(cbind(Y[[1]],Y[[2]],Y[[3]],Y[[4]],
                             A,X,rep(M,M),GP,
                             true.g[[1]],true.g.1[[1]],true.g.0[[1]],
                             true.g[[2]],true.g.1[[2]],true.g.0[[2]],
                             true.g[[3]],true.g.1[[3]],true.g.0[[3]],
                             true.g[[4]],true.g.1[[4]],true.g.0[[4]]) )
    colnames(Data) <- c("Y1","Y2","Y3","Y4",
                        "A","X1","X2","X3","X4","X5","X6M","GP",
                        "g1A","g11","g10","g2A","g21","g20",
                        "g3A","g31","g30","g4A","g41","g40")
    write.csv(Data,sprintf("GenData/Data_N%0.3d_%sPS_B%0.4d.csv",N,sVType,New.BATCH),row.names=FALSE)
}

stopImplicitCluster()
