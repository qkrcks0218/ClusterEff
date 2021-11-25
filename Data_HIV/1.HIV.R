############################################
# Last update : October 09 2021
# Application: ECLS-K
############################################

####################################################
# Cleaning 
# Download data from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/CVOPZL
# Use bio_for_analysis.dta
####################################################

PATH <- "F:/Dropbox/Chan/Research/2021/ClusterEff_Code_Submit/Data_HIV"

setwd(PATH)

library(readstata13)

# RD <- read.dta13("bio_for_analysis.dta")

RD <- RD[RD$eligible_ITsampling==0,]
Vars <- c( "an2_ba_b9_205_more1", # Y
           "sampledVCT", # A
           "SCH_sdkcpe", 
           "SCH_situation",
           "s_Utreat",
           "s_HIVtreat",
           "sex",
           "an2_age2009",
           "an_years_educ",
           "an_curently_maried",
           "ever_children",
           "schoolid") 
loc <- rep(0,length(Vars))
for(jj in 1:length(loc)){
    loc[jj] <- which(colnames(RD)==Vars[jj])
}

Data <- RD[,loc]
Data <- Data[apply(is.na(Data),1,sum)==0,]
Data <- Data[order(Data$schoolid),]

colnames(Data) <- c("Y","A",
                    "C_score","C_loc","C_T1","C_T2",
                    "X_gender","X_age","X_education",
                    "X_married","X_child",
                    "Cid")

Data$C_urban <- as.numeric(Data$C_loc==1)
Data$C_rural <- as.numeric(Data$C_loc==4)
Data$X_gender <- as.numeric(Data$X_gender=="1 Male")
Data$C_size <- rep(table(Data$Cid),table(Data$Cid))

Data <- Data[,c("Y","A",
                "C_score","C_urban","C_rural","C_size","C_T1","C_T2",
                "X_gender","X_age","X_education",
                "X_married","X_child",
                "Cid")]
Data <- Data[order(Data$Cid),]

write.csv(Data,"Data_Cleaned.csv",row.names=FALSE)


####################################################
# Nuisance Function Estimation
####################################################

source("../MySL.R")
source("../ClusterFtSource.R")

TrueData <- read.csv("Data_Cleaned.csv")

N <- length(unique(TrueData$Cid))
M <- table(TrueData$Cid)

TrueData0 <- TrueData
TrueData0[,2] <- 0
TrueData1 <- TrueData
TrueData1[,2] <- 1

index <- list()
for(ii in 1:N){
    index[[ii]] <- which( TrueData$Cid==unique(TrueData$Cid)[ii] )
}

Cov.Pos <- 3:(dim(TrueData)[2]-1)

SizeQ <- c(1,quantile(unique(cbind(TrueData$C_size,TrueData$Cid))[,1],(1:50)/50))
Size <- rep(0,sum(M))
for(jjj in 1:length(Size)){
    Size[jjj] <- sum(TrueData$C_size[jjj]>SizeQ)
}

Sized <- rep(0,N)
for(jjj in 1:length(Sized)){
    Sized[jjj] <- sum(M[jjj]>SizeQ)
}


SL.hpara <- list()
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
SL.hpara$MTRY <- c(5,10)               # random forest parameters
SL.hpara$MLPL <- c(2,4)                # number of nodes in MLP
SL.hpara$NMN <- 25                     # gbm parameter
SL.hpara$MLPdecay <- 10^c(-1,-2)       # MLP decay parameter

for(NF.num in 1:100){
    
    set.seed(NF.num)
    
    SS <- list()
    SMS <- NULL
    US <- sort( unique(Sized) )
    ind.MSAS <- 1
    for(jjj in 1:length(US)){
        POS <- which(Sized==US[jjj])
        if(ind.MSAS==1){
            NMS <- ceiling(length(POS)/2)
            ind.MSAS <- 1-ind.MSAS 
        } else {
            NMS <- floor(length(POS)/2)
            ind.MSAS <- 1-ind.MSAS 
        }
        
        SMS <- c(SMS,sample(POS,NMS))
        
    }
    
    
    SS$MS <- sort( SMS )
    SS$AS <- (1:N)[-SS$MS]
    
    SSI <- list()
    SSI$MS <- NULL
    for(ii in 1:length(SS$MS)){
        SSI$MS <- c(SSI$MS, index[[ SS$MS[ii] ]])
    }
    SSI$AS <- NULL
    for(ii in 1:length(SS$AS)){
        SSI$AS <- c(SSI$AS, index[[ SS$AS[ii] ]])
    }
    
    
    SSIcsv <- matrix(0,sum(M),1)
    SSIcsv[SSI$MS,1] <- 1
    SSIcsv[SSI$AS,1] <- 2
    
    SSIcsv <- data.frame(SSIcsv)
    colnames(SSIcsv) <- c("SS2")
    
    write.csv(SSIcsv,sprintf("Outcome/SSlist_NF%0.4d.csv",NF.num),row.names=FALSE)
    
    OR.fit <- list()
    OR.fit$MS <- OR.Estimation(TrueData$Y[SSI$MS],TrueData$A[SSI$MS],TrueData[SSI$MS,Cov.Pos],TrueData$Cid[SSI$MS],type="ML",outcome.type=binomial(),SL.hpara=SL.hpara)
    OR.fit$AS <- OR.Estimation(TrueData$Y[SSI$AS],TrueData$A[SSI$AS],TrueData[SSI$AS,Cov.Pos],TrueData$Cid[SSI$AS],type="ML",outcome.type=binomial(),SL.hpara=SL.hpara)
    
    OR.est <- list()
    OR.est$ML2 <- rep(0,sum(M))
    OR.est$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,TrueData$A[SSI$AS],TrueData[SSI$AS,Cov.Pos])$output
    OR.est$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,TrueData$A[SSI$MS],TrueData[SSI$MS,Cov.Pos])$output
    
    OR.est1 <- list()
    OR.est1$ML2 <- rep(0,sum(M))
    OR.est1$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,TrueData1$A[SSI$AS],TrueData1[SSI$AS,Cov.Pos])$output
    OR.est1$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,TrueData1$A[SSI$MS],TrueData1[SSI$MS,Cov.Pos])$output
    
    OR.est0 <- list()
    OR.est0$ML2 <- rep(0,sum(M))
    OR.est0$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,TrueData0$A[SSI$AS],TrueData0[SSI$AS,Cov.Pos])$output
    OR.est0$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,TrueData0$A[SSI$MS],TrueData0[SSI$MS,Cov.Pos])$output
    
    RES <- cbind(rep((1:N),M),(1:sum(M)),rep(NF.num,sum(M)),OR.est$ML2,OR.est1$ML2,OR.est0$ML2)
    RES <- data.frame(RES)
    colnames(RES) <- c("Cluster.ID","Ind.ID","BATCH","OR","OR1","OR0")
    
    
    write.csv(RES,sprintf("Outcome/Dupas_OR2_NF%0.4d.csv",NF.num),row.names=FALSE)
    
}

####################################################
# Calculate IF
####################################################


TrueData <- read.csv("Data_Cleaned.csv")

N <- length(unique(TrueData$Cid))
M <- table(TrueData$Cid)

TrueData$X_age <- scale(TrueData$X_age)
TrueData$X_education <- scale(TrueData$X_education)

TrueData0 <- TrueData
TrueData0[,2] <- 0
TrueData1 <- TrueData
TrueData1[,2] <- 1

Cov.Pos <- 3:(dim(TrueData)[2]-1)

OR.fit <- OR.est <- OR.est1 <- OR.est0 <- list()
IPS.est <- IF <- list()

#################################
# Parametric Model
#################################

OR.fit$glmer.noS <- OR.Estimation(TrueData$Y,TrueData$A,TrueData[,Cov.Pos],TrueData$Cid,type="GLMM",outcome.type=binomial(),SL.hpara=SL.hpara)
OR.fit$gee.noS <- OR.Estimation(TrueData$Y,TrueData$A,TrueData[,Cov.Pos],TrueData$Cid,type="GEE",outcome.type=binomial(),SL.hpara=SL.hpara)
OR.fit$glm.noS <- OR.Estimation(TrueData$Y,TrueData$A,TrueData[,Cov.Pos],TrueData$Cid,type="GLM",outcome.type=binomial(),SL.hpara=SL.hpara)

VarCorr(OR.fit$glmer.noS$output) # var=0.4476^2

OR.est$glmer.noS <- OR.Prediction(OR.fit$glmer.noS,TrueData$A,TrueData[,Cov.Pos],outcome.type=binomial())$output
OR.est$gee.noS <- OR.Prediction(OR.fit$gee.noS,TrueData$A,TrueData[,Cov.Pos],outcome.type=binomial())$output
OR.est$glm.noS <- OR.Prediction(OR.fit$glm.noS,TrueData$A,TrueData[,Cov.Pos],outcome.type=binomial())$output
OR.est1$glmer.noS <- OR.Prediction(OR.fit$glmer.noS,TrueData1$A,TrueData1[,Cov.Pos],outcome.type=binomial())$output
OR.est1$gee.noS <- OR.Prediction(OR.fit$gee.noS,TrueData1$A,TrueData1[,Cov.Pos],outcome.type=binomial())$output
OR.est1$glm.noS <- OR.Prediction(OR.fit$glm.noS,TrueData1$A,TrueData1[,Cov.Pos],outcome.type=binomial())$output
OR.est0$glmer.noS <- OR.Prediction(OR.fit$glmer.noS,TrueData0$A,TrueData0[,Cov.Pos],outcome.type=binomial())$output
OR.est0$gee.noS <- OR.Prediction(OR.fit$gee.noS,TrueData0$A,TrueData0[,Cov.Pos],outcome.type=binomial())$output
OR.est0$glm.noS <- OR.Prediction(OR.fit$glm.noS,TrueData0$A,TrueData0[,Cov.Pos],outcome.type=binomial())$output

IPS.est$glmer.noS <- rep(1/2,sum(M))
IPS.est$gee.noS <-   rep(1/2,sum(M))
IPS.est$glm.noS <-   rep(1/2,sum(M))

IF$IPW <- aggregate(TrueData$Y*TrueData$A*2-TrueData$Y*(1-TrueData$A)*2~TrueData$Cid,FUN=mean)[,2]
IF$glmer1 <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$Cid,OR.est$glmer.noS,OR.est1$glmer.noS,OR.est0$glmer.noS,IPS.est$glmer.noS)
IF$gee1 <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$Cid,OR.est$gee.noS,OR.est1$gee.noS,OR.est0$gee.noS,IPS.est$gee.noS)
IF$glm1 <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$Cid,OR.est$glm.noS,OR.est1$glm.noS,OR.est0$glm.noS,IPS.est$glm.noS)

ResY <- (TrueData$Y - OR.est$gee.noS)
Res.Gp.Y <- aggregate(ResY~TrueData$Cid,FUN="sum")[,2]


MSB <- (sum(Res.Gp.Y^2/M) - sum(Res.Gp.Y)^2/sum(M))/(N-1)
MSW <- (sum(ResY^2)-sum(Res.Gp.Y^2/M))/(sum(M)-N)
N.m <- (sum(M) - sum(M^2)/sum(M))/(N-1)
ICCY <- (MSB-MSW)/(MSB+(N.m-1)*MSW) # 0.02199862 # 0.03583591

ResA <- (TrueData$A - IPS.est$glmer.noS)
Res.Gp.A <- aggregate(ResA~TrueData$Cid,FUN="sum")[,2]

MSB <- (sum(Res.Gp.A^2/M) - sum(Res.Gp.A)^2/sum(M))/(N-1)
MSW <- (sum(ResA^2)-sum(Res.Gp.A^2/M))/(sum(M)-N)
N.m <- (sum(M) - sum(M^2)/sum(M))/(N-1)
ICCA <- (MSB-MSW)/(MSB+(N.m-1)*MSW)

ICCY ; ICCA # ICCA

#################################
# Our Method
#################################

RES <- matrix(0,100,12)

for(NF.num in 1:100){
    
    SSIcsv <- read.csv(sprintf("Outcome/SSlist_NF%0.4d.csv",NF.num))
    colnames(SSIcsv) <- "SS2"
    
    SS <- list()
    SS$MS <- which( unique( cbind(rep(1:N,M),SSIcsv$SS2) )[,2]==1 )
    SS$AS <- which( unique( cbind(rep(1:N,M),SSIcsv$SS2) )[,2]==2 )
    
    SSI <- list()
    SSI$MS <- which( SSIcsv$SS2==1 )
    SSI$AS <- which( SSIcsv$SS2==2 )
    
    
    OR2 <- read.csv(sprintf("Outcome/Dupas_OR2_NF%0.4d.csv",NF.num))
    OR.est$ML2 <- OR2$OR
    OR.est1$ML2 <- OR2$OR1
    OR.est0$ML2 <- OR2$OR0
    
    IPS.est$ML2 <- IPS.est$CPS <- rep(1/2,sum(M))
    
    IF$IPS2.B <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$Cid,OR.est$ML2,OR.est1$ML2,OR.est0$ML2,IPS.est$ML2)
    
    IF$CPS.BnoB <- rep(0,N)
    
    SizeQ <- c(1,quantile(unique(cbind(TrueData$C_size,TrueData$Cid))[,1],(1:50)/50))
    Size <- rep(0,sum(M))
    for(jjj in 1:length(Size)){
        Size[jjj] <- sum(TrueData$C_size[jjj]>SizeQ)
    }
    
    Sized <- rep(0,N)
    for(jjj in 1:length(Sized)){
        Sized[jjj] <- sum(M[jjj]>SizeQ)
    }
    
    table( Sized[SS$MS] )
    table( Sized[SS$AS] )
    
    Clu.index <- Size 
    
    Beta <- findbetaDisc(TrueData$Y[SSI$MS],TrueData$A[SSI$MS],TrueData$Cid[SSI$MS],OR.est$ML2[SSI$MS],IPS.est$CPS[SSI$MS],Clustering=Clu.index[SSI$MS])$beta.nobdd
    IF$CPS.BnoB[SS$MS] <- findIFwithB.Disc(TrueData$Y[SSI$MS],TrueData$A[SSI$MS],TrueData$Cid[SSI$MS],OR.est$ML2[SSI$MS],OR.est1$ML2[SSI$MS],OR.est0$ML2[SSI$MS],IPS.est$CPS[SSI$MS],Beta,Clu.index[SSI$MS])
    
    Beta <- findbetaDisc(TrueData$Y[SSI$AS],TrueData$A[SSI$AS],TrueData$Cid[SSI$AS],OR.est$ML2[SSI$AS],IPS.est$CPS[SSI$AS],Clustering=Clu.index[SSI$AS])$beta.nobdd
    IF$CPS.BnoB[SS$AS] <- findIFwithB.Disc(TrueData$Y[SSI$AS],TrueData$A[SSI$AS],TrueData$Cid[SSI$AS],OR.est$ML2[SSI$AS],OR.est1$ML2[SSI$AS],OR.est0$ML2[SSI$AS],IPS.est$CPS[SSI$AS],Beta,Clu.index[SSI$AS])
    
    RES[NF.num,] <- c(sapply(IF,mean),
                      sqrt(var(IF$IPW)/N),
                      sqrt(var(IF$glmer1)/N),
                      sqrt(var(IF$gee1)/N),
                      sqrt(var(IF$glm1)/N),
                      # 0.027254056182642,   
                      # 0.00956706902624382, 
                      # 0.00957322505393492, 
                      # 0.00957880008829728, 
                      # 0.02551257,
                      # 0.01298269,
                      # 0.01297816,
                      # 0.01302969,
                      sqrt(obtainV(IF$IPS2.B,split=2)/N),
                      sqrt(obtainV(IF$CPS.BnoB,split=2)/N))
    
    print(NF.num)
}


RESULT <- data.frame(RES)
colnames(RESULT) <- c(paste("Est_",c("IPW","GLMM","GEE","GLM","ML","CPS_noB"),sep=""),
                      paste("SE_",c("IPW","GLMM","GEE","GLM","ML","CPS_noB"),sep=""))
write.csv(RESULT,"HIV_Ours.csv",row.names=FALSE)


######################
# R-learner/GRF
######################

library(rlearner)
library(grf)
library(foreach)
library(doParallel)

RRR <- matrix(0,100,10)
RRR <- data.frame(RRR)
colnames(RRR) <- c("Est_R_Lasso",
                   "Est_R_Boost",
                   "Est_U_Lasso",
                   "Est_U_Boost",
                   "Est_GRF",
                   "SE_R_Lasso",
                   "SE_R_Boost",
                   "SE_U_Lasso",
                   "SE_U_Boost",
                   "SE_GRF")

for(ii in 1:100){
    
    index <- list()
    for(jjj in 1:N){
        index[[jjj]] <- which( TrueData$Cid==unique(TrueData$Cid)[jjj] )
    }
    
    Cov.Pos <- 3:(dim(TrueData)[2]-1)
    
    
    set.seed(100*ii + ii)
    
    RL <- rlasso(as.matrix(TrueData[,Cov.Pos]), TrueData[,2], TrueData[,1], p_hat=rep(1/2,sum(M)))
    
    output11 <- ATEcalculate(TrueData, RL, weight.vector=1/rep(M,M), clusters=TrueData[,dim(TrueData)[2]])
    
    RLest <- c( output11[1])
    RLsd <- c( output11[2] )
    
    
    
    RB <- rboost(as.matrix(TrueData[,Cov.Pos]), TrueData[,2], TrueData[,1], p_hat=rep(1/2,sum(M)))
    
    output21 <- ATEcalculate(TrueData, RB, weight.vector=1/rep(M,M), clusters=TrueData[,dim(TrueData)[2]])
    
    RBest <- c( output21[1])
    RBsd <- c( output21[2] )
    
    
    
    UL <- ulasso(as.matrix(TrueData[,Cov.Pos]), TrueData[,2], TrueData[,1], p_hat=rep(1/2,sum(M)))
    
    output31 <- ATEcalculate(TrueData, UL, weight.vector=1/rep(M,M), clusters=TrueData[,dim(TrueData)[2]])
    
    ULest <- c( output31[1])
    ULsd <- c( output31[2] )
    
    
    
    UB <- uboost(as.matrix(TrueData[,Cov.Pos]), TrueData[,2], TrueData[,1], p_hat=rep(1/2,sum(M)))
    
    output41 <- ATEcalculate(TrueData, UB, weight.vector=1/rep(M,M), clusters=TrueData[,dim(TrueData)[2]])
    
    UBest <- c( output41[1])
    UBsd <- c( output41[2] )
    
    CF <- grf::causal_forest(X=TrueData[,Cov.Pos],
                             Y=TrueData[,1],
                             W=TrueData[,2],
                             W.hat=rep(1/2,sum(M)),
                             sample.weights=1/rep(M,M),
                             clusters=TrueData[,dim(TrueData)[2]])
    CFest <- grf::average_treatment_effect(CF)[1]
    CFsd <- grf::average_treatment_effect(CF)[2]
    
    
    RRR[ii,] <- c(RLest,RBest,ULest,UBest, CFest,
                  RLsd, RBsd, ULsd, UBsd , CFsd )
}

write.csv(RRR,"HIV_RUGRF.csv",row.names=FALSE)

####################################################
# Summary
####################################################

R1 <- read.csv("HIV_Ours.csv")
R1 <- RESULT
Est.R1 <- apply(R1,2,median)[1:(dim(R1)[2]/2)]
SE.R1 <- sqrt(apply(R1[,dim(R1)[2]/2+1:(dim(R1)[2]/2)]^2 + 
                        (R1[,1:(dim(R1)[2]/2)] - matrix(Est.R1,dim(R1),(dim(R1)[2]/2),byrow=T))^2,2,median))

R2 <- read.csv("HIV_RUGRF.csv")

Est.R2 <- apply(R2,2,median)[1:5]
SE.R2 <- sqrt(apply(R2[,6:10]^2 + (R2[,1:5] - matrix(Est.R2,dim(R2),5,byrow=T))^2,2,median))

print(data.frame(rbind(sprintf("Estimate & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f \\\\ \\hline",
                               Est.R1[2],Est.R1[3],Est.R1[4],
                               Est.R2[2],Est.R2[1],Est.R2[4],Est.R2[3],Est.R2[5],
                               Est.R1[5],Est.R1[6]),
                       sprintf("SE & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f \\\\ \\hline",
                               SE.R1[2],SE.R1[3],SE.R1[4],
                               SE.R2[2],SE.R2[1],SE.R2[4],SE.R2[3],SE.R2[5],
                               SE.R1[5],SE.R1[6]),
                       sprintf("RE($\\widehat{\\tau}_w,\\cdot$) & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f & %0.4f \\\\ \\hline",
                               (SE.R1[2]/SE.R1[6])^2,(SE.R1[3]/SE.R1[6])^2,(SE.R1[4]/SE.R1[6])^2,
                               (SE.R2[2]/SE.R1[6])^2,(SE.R2[1]/SE.R1[6])^2,(SE.R2[4]/SE.R1[6])^2,(SE.R2[3]/SE.R1[6])^2,(SE.R2[5]/SE.R1[6])^2,
                               (SE.R1[5]/SE.R1[6])^2,(SE.R1[6]/SE.R1[6])^2) )),row.names=FALSE)









