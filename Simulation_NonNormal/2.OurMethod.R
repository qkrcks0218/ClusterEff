############################################
# Last update : October 09 2021
# Estimation from Our Method
############################################

PATH <- getwd()

library(lme4)

source("../MySL.R")
source("../ClusterFtSource.R")


CPS.B <- 5                                  # Number of CPS estimation
Tot.NF <- 5                                 # Number of cross-fitting procedures
SL.hpara <- list()                          # Super Learner parameters
SL.hpara$SLL <- c(1,2,3,4,5,6,7,9,10)
# Superlearner basic learning algorithms: 
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 9: gbm
# 10: 1-layer MLP
SL.hpara$MTRY <- c(2,4)                # random forest parameters
SL.hpara$MLPL <- c(2,4)                # number of nodes in MLP
SL.hpara$NMN <- 25                     # gbm parameter
SL.hpara$MLPdecay <- 10^c(-3,-4,-5)    # MLP decay parameter



Para.Mat1 <- expand.grid(c(500),c(0,1.5)) # 

for(BATCH1 in 1:2){          # Recommend Parallel/Cluster Computing
    
    sV <- Para.Mat1[BATCH1,2]                 # \sigma_V
    N <- Para.Mat1[BATCH1,1]                  # Cluster size
    if(sV==0){
        sVType <- "No"
    } else {
        sVType <- "Strong"
    }
    
    if(N>10){
        pos.X <- 3:8
        pos.AX <- 2:11
    } else {
        pos.X <- c(3,4,5,6,8)               # Remove one binary variable due to singularity
        pos.AX <- c(2,3,4,5,6,8)            # Remove one binary variable due to singularity
    }
    
    Para.Mat2 <- expand.grid(1:Tot.NF,1:200)
    RRR1 <- matrix(0,1000,10)
    RRR1 <- RRR2 <- data.frame(RRR1)
    colnames(RRR1) <- colnames(RRR2) <- c( "N","sV","Data.Seed","SS.Seed",
                                           paste("Eff_Est_",c("IPS","CPS_G"),sep=""),
                                           paste("SE_Est_",c("IPS","CPS_G"),sep=""),
                                           paste("Cover_Est_",c("IPS","CPS_G"),sep="") )
    
    for(BATCH2 in 1:1000){
        
        Sim.Seed <- Para.Mat[BATCH2,2]           # Data Number
        Num.NF <- Para.Mat[BATCH2,1]             # Cross-fitting 
        
        RawData <- read.csv(sprintf("GenData/Data_N%0.3d_%sPS_B%0.4d.csv",N,sVType,Sim.Seed))
        varname <- colnames(RawData)
        N <- length(unique(RawData$GP))
        M <- as.numeric( table(RawData$GP) )
        M.start <- cumsum(M)-M+1
        M.end <- cumsum(M)
        GP <- RawData$GP
        A <- RawData$A
        X <- RawData[,substring(varname,1,1)=="X"]
        
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
        
        Data <- DataY[[1]]
        
        DataA <- data.frame( cbind(A,
                                   as.numeric((rep(aggregate(A~GP,data=Data,FUN="sum")$A,M)-A)/(rep(M,M)-1)),
                                   X,
                                   as.numeric((rep(aggregate(X[,1]~GP,data=Data,FUN="sum")[,2],M)-X[,1])/(rep(M,M)-1)),
                                   as.numeric((rep(aggregate(X[,2]~GP,data=Data,FUN="sum")[,2],M)-X[,2])/(rep(M,M)-1)),
                                   as.numeric((rep(aggregate(X[,3]~GP,data=Data,FUN="sum")[,2],M)-X[,3])/(rep(M,M)-1)),
                                   GP) )
        
        colnames(DataA) <- c("A","OA","X1","X2","X3","X4","X5","X6M","OX1","OX2","OX3","GP")
        
        
        ##########################################
        # Split Sample
        ##########################################
        
        set.seed( Num.NF+Sim.Seed*1000 )
        
        temp <- sample(N,N)
        
        SSMS.C <- sort( temp[1:round(N/2)] )
        SSAS.C <- (1:N)[-SSMS.C]
        
        SSM <- SSA <- NULL
        for(ii in 1:N){
            pos <- which(GP==ii)
            if(sum(SSMS.C==ii)>0){ SSM <- c(SSM,pos) } else if (sum(SSAS.C==ii)) { SSA <- c(SSA,pos) } 
        }
        
        SS <- list()
        SS$MS <- SSMS.C
        SS$AS <- SSAS.C
        
        SSI <- list()
        SSI$MS <- SSM
        SSI$AS <- SSA
        
        ##########################################
        # PS Estimation
        ##########################################
        
        PS.fit <- list()
        IPS.est <- list()
        GPS.est <- list()
        OPS.est <- list()
        
        PS.fit$MS <- IPS.Estimation(Data$A[SSI$MS],Data[SSI$MS,pos.X],Data$GP[SSI$MS],type="ML",SL.hpara=SL.hpara)
        PS.fit$AS <- IPS.Estimation(Data$A[SSI$AS],Data[SSI$AS,pos.X],Data$GP[SSI$AS],type="ML",SL.hpara=SL.hpara)
        
        IPS.est$ML2 <- rep(0,sum(M))
        IPS.est$ML2[SSI$AS] <- IPS.Prediction(PS.fit$MS,Data$A[SSI$AS],DataA[SSI$AS,pos.X],Data$GP[SSI$AS])$IPS
        IPS.est$ML2[SSI$MS] <- IPS.Prediction(PS.fit$AS,Data$A[SSI$MS],DataA[SSI$MS,pos.X],Data$GP[SSI$MS])$IPS
        
        
        IPS.est$CPS <- rep(0,sum(M))
        CPS.temp.mat <- matrix(0,sum(M),CPS.B)
        
        for(cps.iter in 1:CPS.B){
            
            SSI.temp <- list()
            SSI.temp$MS <- SSI.temp$AS <- NULL
            for(ii in 1:N){
                pos <- which(rep(1:N, M)==ii)
                pos <- sample(pos,min(M))
                if(sum(SSMS.C==ii)>0){ 
                    SSI.temp$MS <- c(SSI.temp$MS,pos) 
                } else { 
                    SSI.temp$AS <- c(SSI.temp$AS,pos) 
                }
            }
            SSI.temp$MS <- sort(SSI.temp$MS)
            SSI.temp$AS <- sort(SSI.temp$AS)
            
            DataB <- DataA
            SL.hpara2 <- SL.hpara
            
            CVL.MS <- CVL.AS <- list()
            
            if(N>10){
                MSind <- c(round(quantile(1:length(SS$MS),(1:5)/5)))
                rMS <- SS$MS[ order(runif(length(SS$MS))) ]
                MSCind <- list()
                MSCind[[1]] <- rMS[1:MSind[1]]
                MSCind[[2]] <- rMS[(MSind[1]+1):MSind[2]]
                MSCind[[3]] <- rMS[(MSind[2]+1):MSind[3]]
                MSCind[[4]] <- rMS[(MSind[3]+1):MSind[4]]
                MSCind[[5]] <- rMS[(MSind[4]+1):MSind[5]]
                
                for(cvind in 1:5){
                    TempCV <- NULL
                    for(jjj in 1:length(MSCind[[cvind]])){
                        TempCV <- c(TempCV,which( DataB$GP[SSI.temp$MS]==MSCind[[cvind]][jjj] ))
                    }
                    CVL.MS[[cvind]] <- TempCV
                }
                
                ASind <- c(round(quantile(1:length(SS$AS),(1:5)/5)))
                rAS <- SS$AS[ order(runif(length(SS$AS))) ]
                ASCind <- list()
                ASCind[[1]] <- rAS[1:ASind[1]]
                ASCind[[2]] <- rAS[(ASind[1]+1):ASind[2]]
                ASCind[[3]] <- rAS[(ASind[2]+1):ASind[3]]
                ASCind[[4]] <- rAS[(ASind[3]+1):ASind[4]]
                ASCind[[5]] <- rAS[(ASind[4]+1):ASind[5]]
                
                for(cvind in 1:5){
                    TempCV <- NULL
                    for(jjj in 1:length(ASCind[[cvind]])){
                        TempCV <- c(TempCV,which( DataB$GP[SSI.temp$AS]==ASCind[[cvind]][jjj] ))
                    }
                    CVL.AS[[cvind]] <- TempCV
                } 
            } else {
                MSCind <- list()
                MSCind[[1]] <- SS$MS[1]
                MSCind[[2]] <- SS$MS[2]
                MSCind[[3]] <- SS$MS[3]
                MSCind[[4]] <- SS$MS[4]
                MSCind[[5]] <- SS$MS[5]
                
                for(cvind in 1:5){
                    TempCV <- NULL
                    for(jjj in 1:length(MSCind[[cvind]])){
                        TempCV <- c(TempCV,which( DataB$GP[SSI.temp$MS]==MSCind[[cvind]][jjj] ))
                    }
                    CVL.MS[[cvind]] <- TempCV
                }
                
                ASCind <- list()
                ASCind[[1]] <- SS$AS[1]
                ASCind[[2]] <- SS$AS[2]
                ASCind[[3]] <- SS$AS[3]
                ASCind[[4]] <- SS$AS[4]
                ASCind[[5]] <- SS$AS[5]
                
                for(cvind in 1:5){
                    TempCV <- NULL
                    for(jjj in 1:length(ASCind[[cvind]])){
                        TempCV <- c(TempCV,which( DataB$GP[SSI.temp$AS]==ASCind[[cvind]][jjj] ))
                    }
                    CVL.AS[[cvind]] <- TempCV
                }
            }
            
            
            
            PS.fit$CPS.MS <- IPS.Estimation(DataB$A[SSI.temp$MS],DataB[SSI.temp$MS,pos.AX],DataB$GP[SSI.temp$MS],type="ML",SL.hpara=SL.hpara2,CVlist=CVL.MS)
            PS.fit$CPS.AS <- IPS.Estimation(DataB$A[SSI.temp$AS],DataB[SSI.temp$AS,pos.AX],DataB$GP[SSI.temp$AS],type="ML",SL.hpara=SL.hpara2,CVlist=CVL.AS)
            
            CPS.temp.mat[SSI$AS,cps.iter] <- IPS.Prediction(PS.fit$CPS.MS,DataA$A[SSI$AS],DataA[SSI$AS,pos.AX],DataA$GP[SSI$AS])$IPS
            CPS.temp.mat[SSI$MS,cps.iter] <- IPS.Prediction(PS.fit$CPS.AS,DataA$A[SSI$MS],DataA[SSI$MS,pos.AX],DataA$GP[SSI$MS])$IPS
            
            print(cps.iter)
        }
        
        IPS.est$CPS <- apply(CPS.temp.mat,1,median)
        
        ##########################################
        # OR Estimation 
        ##########################################
        
        for(or.iter in 1:2){
            
            if(or.iter<=2){
                OT <- gaussian()
            } else {
                OT <- binomial()
            }
            
            Data <- DataY[[or.iter]]
            TRUE.G <- true.g[[or.iter]]
            TRUE.G.1 <- true.g.1[[or.iter]]
            TRUE.G.0 <- true.g.0[[or.iter]]
            
            OR.fit <- list()
            OR.est <- OR.est1 <- OR.est0 <- list()
            
            OR.fit$MS <- OR.Estimation(Data$Y[SSI$MS],Data$A[SSI$MS],Data[SSI$MS,pos.X],Data$GP[SSI$MS],type="ML",outcome.type=OT,SL.hpara=SL.hpara)
            OR.fit$AS <- OR.Estimation(Data$Y[SSI$AS],Data$A[SSI$AS],Data[SSI$AS,pos.X],Data$GP[SSI$AS],type="ML",outcome.type=OT,SL.hpara=SL.hpara)
            
            OR.est$ML2 <- rep(0,sum(M))
            OR.est$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,Data$A[SSI$AS],Data[SSI$AS,pos.X])$output
            OR.est$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,Data$A[SSI$MS],Data[SSI$MS,pos.X])$output
            
            OR.est1$ML2 <- rep(0,sum(M))
            OR.est1$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,rep(1,length(SSI$AS)),Data[SSI$AS,pos.X])$output
            OR.est1$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,rep(1,length(SSI$MS)),Data[SSI$MS,pos.X])$output
            
            OR.est0$ML2 <- rep(0,sum(M))
            OR.est0$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,rep(0,length(SSI$AS)),Data[SSI$AS,pos.X])$output
            OR.est0$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,rep(0,length(SSI$MS)),Data[SSI$MS,pos.X])$output
            
            
            ##########################################
            # Estimate
            ##########################################
            
            RESULT <- list()
            RESULT$Effect <- rep(0,2)
            RESULT$SE <- rep(0,2)
            RESULT$COVER <- rep(0,2)
            
            IF.LIST <- list()
            
            IF.LIST$ML2 <- findIFwithoutB(Data$Y,Data$A,Data$GP,OR.est$ML2,OR.est1$ML2,OR.est0$ML2,IPS.est$ML2)
            RESULT$Effect[1] <- mean( IF.LIST$ML2 )
            RESULT$SE[1] <- sqrt(obtainV(IF.LIST$ML2,split=2)/N)
            
            
            Clustering <- rep(M,M)
            if(N<200){
                Tnum <- 1*as.numeric(N==25) + 2*as.numeric(N==50) + 3*as.numeric(N==100)  
                QM <- c(0,quantile(M,(1:Tnum)/Tnum))
                for(qn in 1:Tnum){
                    pos.qn <- which(Clustering<=QM[qn+1] & Clustering>QM[qn])
                    Clustering[pos.qn] <- qn
                }
            }
            
            
            Beta.MS <- findbetaDisc(Data$Y[SSI$MS],Data$A[SSI$MS],Data$GP[SSI$MS],OR.est$ML2[SSI$MS],IPS.est$CPS[SSI$MS],Clustering=Clustering[SSI$MS])$beta.nobdd
            Beta.AS <- findbetaDisc(Data$Y[SSI$AS],Data$A[SSI$AS],Data$GP[SSI$AS],OR.est$ML2[SSI$AS],IPS.est$CPS[SSI$AS],Clustering=Clustering[SSI$AS])$beta.nobdd
            
            IF.LIST$CPSG[SS$MS] <- findIFwithB.Disc(Data$Y[SSI$MS],Data$A[SSI$MS],Data$GP[SSI$MS],OR.est$ML2[SSI$MS],OR.est1$ML2[SSI$MS],OR.est0$ML2[SSI$MS],IPS.est$CPS[SSI$MS],Beta.MS,Clustering=Clustering[SSI$MS])
            IF.LIST$CPSG[SS$AS] <- findIFwithB.Disc(Data$Y[SSI$AS],Data$A[SSI$AS],Data$GP[SSI$AS],OR.est$ML2[SSI$AS],OR.est1$ML2[SSI$AS],OR.est0$ML2[SSI$AS],IPS.est$CPS[SSI$AS],Beta.AS,Clustering=Clustering[SSI$AS])
            
            
            RESULT$Effect[2] <- mean( IF.LIST$CPSG )
            RESULT$SE[2] <- sqrt(obtainV(IF.LIST$CPSG,split=2)/N)
            
            RESULT$COVER <- as.numeric(apply( cbind(RESULT$Effect - qnorm(0.975)*RESULT$SE - ATE[[or.iter]] , RESULT$Effect + qnorm(0.975)*RESULT$SE-ATE[[or.iter]] ), 1, prod )<0)
            
            
            
            if(or.iter<=2){
                RRR1[BATCH2,] <- c(N, sV, Sim.Seed, Num.NF, RESULT$Effect,RESULT$SE,RESULT$COVER)
            } else {
                RRR2[BATCH2,] <- c(N, sV, Sim.Seed, Num.NF, RESULT$Effect,RESULT$SE,RESULT$COVER)
            }
            
        }
    }
    
    write.csv(RRR1,sprintf("Ours/N%0.4d_%sPS_G1.csv",N,sVType),row.names=FALSE)
    write.csv(RRR2,sprintf("Ours/N%0.4d_%sPS_G2.csv",N,sVType),row.names=FALSE)
    
}

