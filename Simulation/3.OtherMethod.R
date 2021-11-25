############################################
# Last update : October 09 2021
# Estimation from Other Method
############################################

PATH <- "F:/Dropbox/Chan/Research/2021/ClusterEff_Code_Submit/Simulation"
setwd(PATH)

library(rlearner)
library(grf)
library(lme4)
library(geepack)


source("../MySL.R")
source("../ClusterFtSource.R")

for(BATCH in 1:400){
    
    
    Tot.NF <- 1
    pos.X <- 3:8
    
    
    
    Para.Mat <- expand.grid(c(500),1:Tot.NF,1:210,c(0,1.5))
    sV <- Para.Mat[BATCH,4]
    Sim.Seed <- Para.Mat[BATCH,3]
    Num.NF <- Para.Mat[BATCH,2]
    N <- Para.Mat[BATCH,1]
    if(sV==0){
        sVType <- "No"
    } else {
        sVType <- "Strong"
    }
    
    RawData <- read.csv(sprintf("GenData/Data_N%0.3d_%sPS_B%0.4d.csv",N,sVType,Sim.Seed))
    varname <- colnames(RawData)
    N <- length(unique(RawData$GP))
    M <- as.numeric( table(RawData$GP) )
    M.start <- cumsum(M)-M+1
    M.end <- cumsum(M)
    GP <- RawData$GP
    A <- RawData$A
    X <- as.matrix( RawData[,substring(varname,1,1)=="X"] )
    
    
    Y <- list()
    true.g <- list()
    true.g.1 <- list()
    true.g.0 <- list()
    
    for(jj in 1:4){
        Y[[jj]] <- RawData[,varname==paste("Y",jj,sep="")]
        true.g[[jj]] <- RawData[,varname==paste("g",jj,"A",sep="")]
        true.g.1[[jj]] <- RawData[,varname==paste("g",jj,"1",sep="")]
        true.g.0[[jj]] <- RawData[,varname==paste("g",jj,"0",sep="")]
    }
    sU <- list()
    sY <- list()
    sU[[1]] <- sU[[3]] <- 0.5
    sU[[2]] <- sU[[4]] <- 1.5
    sY[[1]] <- sY[[2]] <- sY[[3]] <- sY[[4]] <- 1
    
    ATE <- list()
    ATE[[1]] <- ATE[[2]] <- 4
    ATE[[3]] <- 0.1252752 
    ATE[[4]] <- 0.1016713 
    
    DataY <- list()
    DataY[[1]] <- data.frame( cbind(Y[[1]],A,X,GP) )
    colnames(DataY[[1]]) <- c("Y","A","X1","X2","X3","X4","X5","X6M","GP")
    
    DataY[[2]] <- data.frame( cbind(Y[[2]],A,X,GP) )
    colnames(DataY[[2]]) <- c("Y","A","X1","X2","X3","X4","X5","X6M","GP")
    
    DataY[[3]] <- data.frame( cbind(Y[[3]],A,X,GP) )
    colnames(DataY[[3]]) <- c("Y","A","X1","X2","X3","X4","X5","X6M","GP")
    
    DataY[[4]] <- data.frame( cbind(Y[[4]],A,X,GP) )
    colnames(DataY[[4]]) <- c("Y","A","X1","X2","X3","X4","X5","X6M","GP")
    
    
    RESULT.GLMM <- rep(0,8)
    RESULT.GEE <- rep(0,8)
    RESULT.RB <- rep(0,8)
    RESULT.RL <- rep(0,8)
    RESULT.UB <- rep(0,8)
    RESULT.UL <- rep(0,8)
    RESULT.GRF <- rep(0,8)
    
    ##########################################
    # GLMM/GLM
    ##########################################
    
    GLMM.OR <- list()
    
    GLMM.OR[[1]] <- lmer(Y~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5+(1|GP),data=DataY[[1]])
    GLMM.OR[[2]] <- lmer(Y~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5+(1|GP),data=DataY[[2]])
    
    Design.Mat <- Design.Mat0 <- Design.Mat1 <- list()
    for(bb in 1:2){
        Design.Mat[[bb]] <- model.matrix(GLMM.OR[[bb]])
        Design.Mat0[[bb]] <- model.matrix(GLMM.OR[[bb]])
        Design.Mat0[[bb]][,c(2,5,6)] <- 0
        Design.Mat1[[bb]] <- model.matrix(GLMM.OR[[bb]])
        Design.Mat1[[bb]][,2] <- 1
        Design.Mat1[[bb]][,5] <- DataY[[1]]$X2^2
        Design.Mat1[[bb]][,6] <- DataY[[1]]$X3
    }
    
    
    OR.Est.GLMM <- OR.Est0.GLMM <- OR.Est1.GLMM <- list()
    for(bb in 1:2){
        OR.Est.GLMM[[bb]] <- Design.Mat[[bb]]%*%fixef(GLMM.OR[[bb]])
        OR.Est0.GLMM[[bb]] <- Design.Mat0[[bb]]%*%fixef(GLMM.OR[[bb]])
        OR.Est1.GLMM[[bb]] <- Design.Mat1[[bb]]%*%fixef(GLMM.OR[[bb]])
    }
    
    GLMM.PS <- glmer(A~X1+I(X2>1)+X3+X4+X5+(1|GP),data=DataY[[1]],family="binomial")
    coef.PS <- fixef(GLMM.PS)
    sV.GLMM <- as.data.frame(VarCorr(GLMM.PS))[1,5]
    
    PS.function.GLMM <- function(X,V){
        return( expit( coef.PS[1] + coef.PS[2]*X[,1] + coef.PS[3]*as.numeric(X[,2]>1) + coef.PS[4]*X[,3] + coef.PS[5]*X[,4] + coef.PS[6]*X[,5] + V )  )
    }
    
    GPS.GLMM <- JointPS(A,X,PS.function.GLMM,GP,sV.GLMM) 
    IPS.GLMM <- JointPS(A,X,PS.function.GLMM,1:sum(M),sV.GLMM)
    
    IF.LIST.GLMM <- list()
    
    for(bb in 1:2){
        IF.LIST.GLMM[[bb]] <- findIFwithoutB(DataY[[bb]]$Y,DataY[[bb]]$A,DataY[[bb]]$GP,OR.Est.GLMM[[bb]],OR.Est1.GLMM[[bb]],OR.Est0.GLMM[[bb]],IPS.GLMM)
        RESULT.GLMM[bb] <- mean( IF.LIST.GLMM[[bb]] ) 
        RESULT.GLMM[bb+4] <- sd( IF.LIST.GLMM[[bb]] )/sqrt(N)
    }
    
    
    ##########################################
    # GEE
    ##########################################
    
    GEE.OR <- list()
    GEE.OR[[1]] <- geeglm(Y~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5,id=GP,data=DataY[[1]],corstr="exchangeable")
    GEE.OR[[2]] <- geeglm(Y~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5,id=GP,data=DataY[[2]],corstr="exchangeable")
    
    Design.Mat <- Design.Mat0 <- Design.Mat1 <- list()
    for(bb in 1:2){
        Design.Mat[[bb]] <- model.matrix(~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5,DataY[[bb]])
        Design.Mat0[[bb]] <- model.matrix(~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5,DataY[[bb]])
        Design.Mat0[[bb]][,c(2,5,6)] <- 0
        Design.Mat1[[bb]] <- model.matrix(~A+I(X2^2):A+X3:A+X1+I(X4^2)+X2:X5,DataY[[bb]])
        Design.Mat1[[bb]][,2] <- 1
        Design.Mat1[[bb]][,5] <- DataY[[1]]$X2^2
        Design.Mat1[[bb]][,6] <- DataY[[1]]$X3
    }
    
    OR.Est.GEE <- OR.Est0.GEE <- OR.Est1.GEE <- list()
    for(bb in 1:2){
        OR.Est.GEE[[bb]] <-  Design.Mat[[bb]]%*%GEE.OR[[bb]]$coefficients
        OR.Est0.GEE[[bb]] <- Design.Mat0[[bb]]%*%GEE.OR[[bb]]$coefficients
        OR.Est1.GEE[[bb]] <- Design.Mat1[[bb]]%*%GEE.OR[[bb]]$coefficients
    }
    
    GEE.PS <- geeglm(A~X1+I(X2>1)+X3+X4+X5,id=GP,data=DataY[[1]],family=binomial,corstr="exchangeable")
    coef.PS <- GEE.PS$coefficients
    
    PS.function.GEE <- function(X,V){
        return( expit( coef.PS[1] + coef.PS[2]*X[,1] + coef.PS[3]*as.numeric(X[,2]>1) + coef.PS[4]*X[,3] + coef.PS[5]*X[,4] + coef.PS[6]*X[,5] + V)  )
    }
    
    IPS.GEE <- JointPS(A,X,PS.function.GEE,1:sum(M),0)
    
    IF.LIST.GEE <- list()
    
    for(bb in 1:2){
        IF.LIST.GEE[[bb]] <- findIFwithoutB(DataY[[bb]]$Y,DataY[[bb]]$A,DataY[[bb]]$GP,OR.Est.GEE[[bb]],OR.Est1.GEE[[bb]],OR.Est0.GEE[[bb]],IPS.GEE)
        RESULT.GEE[bb] <- mean( IF.LIST.GEE[[bb]] ) 
        RESULT.GEE[bb+4] <- sd( IF.LIST.GEE[[bb]] )/sqrt(N)
    }
    
    
    ##########################################
    # R/U learner
    ##########################################
    
    RB <- RL <- UB <- UL <- list()
    
    RB[[1]] <- rboost(as.matrix(DataY[[1]][,pos.X]),DataY[[1]][,2],DataY[[1]][,1],k_folds=5)
    RB[[2]] <- rboost(as.matrix(DataY[[2]][,pos.X]),DataY[[2]][,2],DataY[[2]][,1],k_folds=5,p_hat=RB[[1]]$p_hat)
    
    RL[[1]] <- rlasso(as.matrix(DataY[[1]][,pos.X]),DataY[[1]][,2],DataY[[1]][,1],k_folds=5)
    RL[[2]] <- rlasso(as.matrix(DataY[[2]][,pos.X]),DataY[[2]][,2],DataY[[2]][,1],k_folds=5,p_hat=RL[[1]]$p_hat)
    
    UB[[1]] <- uboost(as.matrix(DataY[[1]][,pos.X]),DataY[[1]][,2],DataY[[1]][,1],k_folds=5,p_hat=RB[[1]]$p_hat,m_hat=RB[[1]]$m_hat)
    UB[[2]] <- uboost(as.matrix(DataY[[2]][,pos.X]),DataY[[2]][,2],DataY[[2]][,1],k_folds=5,p_hat=RB[[1]]$p_hat,m_hat=RB[[2]]$m_hat)
    
    UL[[1]] <- ulasso(as.matrix(DataY[[1]][,pos.X]),DataY[[1]][,2],DataY[[1]][,1],k_folds=5,p_hat=RL[[1]]$p_hat,m_hat=RL[[1]]$m_hat)
    UL[[2]] <- ulasso(as.matrix(DataY[[2]][,pos.X]),DataY[[2]][,2],DataY[[2]][,1],k_folds=5,p_hat=RL[[1]]$p_hat,m_hat=RL[[2]]$m_hat)
    
    ##########################################
    # GRF
    ##########################################
    
    GRF <- list()
    
    GRF[[1]] <- causal_forest(DataY[[1]][,pos.X],DataY[[1]]$Y,DataY[[1]]$A, sample.weights=1/DataY[[1]]$X6M, clusters=DataY[[1]]$GP)
    GRF[[2]] <- causal_forest(DataY[[2]][,pos.X],DataY[[2]]$Y,DataY[[2]]$A, sample.weights=1/DataY[[2]]$X6M, clusters=DataY[[2]]$GP)
    
    for(bb in 1:2){
        RESULT.RB[c(0,4)+bb] <- ATEcalculate(DataY[[bb]], RB[[bb]], weight.vector=1/DataY[[bb]]$X6M, clusters=DataY[[bb]]$GP)
        RESULT.RL[c(0,4)+bb] <- ATEcalculate(DataY[[bb]], RL[[bb]], weight.vector=1/DataY[[bb]]$X6M, clusters=DataY[[bb]]$GP)
        RESULT.UB[c(0,4)+bb] <- ATEcalculate(DataY[[bb]], UB[[bb]], weight.vector=1/DataY[[bb]]$X6M, clusters=DataY[[bb]]$GP)
        RESULT.UL[c(0,4)+bb] <- ATEcalculate(DataY[[bb]], UL[[bb]], weight.vector=1/DataY[[bb]]$X6M, clusters=DataY[[bb]]$GP)
        RESULT.GRF[c(0,4)+bb] <- average_treatment_effect(GRF[[bb]])
    }
    
    RRR <- matrix(c(RESULT.GLMM[1:4],RESULT.GEE[1:4],RESULT.RB[1:4],RESULT.RL[1:4],RESULT.UB[1:4],RESULT.UL[1:4],RESULT.GRF[1:4],
                    RESULT.GLMM[5:8],RESULT.GEE[5:8],RESULT.RB[5:8],RESULT.RL[5:8],RESULT.UB[5:8],RESULT.UL[5:8],RESULT.GRF[5:8]),1,56)
    
    COVER <- matrix(as.numeric( RRR[,1:28] - qnorm(0.975)*RRR[,29:56] <= rep(sapply(ATE,max),7) &
                                    rep(sapply(ATE,max),7) <= RRR[,1:28] + qnorm(0.975)*RRR[,29:56] ),1,28)
    
    RRR <- data.frame(cbind(RRR,COVER))
    colnames(RRR) <- c( paste("Est_",c(paste("GLMM_",1:4,sep=""),
                                       paste("GEE_",1:4,sep=""),
                                       paste("RBoost_",1:4,sep=""),
                                       paste("RLasso_",1:4,sep=""),
                                       paste("UBoost_",1:4,sep=""),
                                       paste("ULasso_",1:4,sep=""),
                                       paste("GRF_",1:4,sep="")),sep=""),
                        paste("SE_",c(paste("GLMM_",1:4,sep=""),
                                      paste("GEE_",1:4,sep=""),
                                      paste("RBoost_",1:4,sep=""),
                                      paste("RLasso_",1:4,sep=""),
                                      paste("UBoost_",1:4,sep=""),
                                      paste("ULasso_",1:4,sep=""),
                                      paste("GRF_",1:4,sep="")),sep=""),
                        paste("COVER_",c(paste("GLMM_",1:4,sep=""),
                                         paste("GEE_",1:4,sep=""),
                                         paste("RBoost_",1:4,sep=""),
                                         paste("RLasso_",1:4,sep=""),
                                         paste("UBoost_",1:4,sep=""),
                                         paste("ULasso_",1:4,sep=""),
                                         paste("GRF_",1:4,sep="")),sep=""))
    write.csv(RRR,sprintf("Others/Others_N%0.4d_sV_%s_B%0.4d.csv",N,sVType,Sim.Seed),row.names=FALSE)
    
}




