############################################
# Last update : October 09 2021
# Summary of the Results
############################################

PATH <- getwd()

ATE <- 4

Result.Ours <- list()

Result.Ours$NPS.G1   <- Result.Ours$SPS.G1   <- Result.Ours$NPS.G2   <- Result.Ours$SPS.G2   <- list()
TYPE <- c("G1","G2","B1","B2")

TempG1.NPS <- TempG2.NPS <- TempG1.SPS <- TempG2.SPS <- list()

TempG1.NPS[[ 2]] <- read.csv("Ours/N0025_NoPS_G1.csv")
TempG1.NPS[[ 3]] <- read.csv("Ours/N0050_NoPS_G1.csv")
TempG1.NPS[[ 4]] <- read.csv("Ours/N0100_NoPS_G1.csv")
TempG1.NPS[[ 5]] <- read.csv("Ours/N0250_NoPS_G1.csv")
TempG1.NPS[[ 6]] <- read.csv("Ours/N0500_NoPS_G1.csv")
TempG1.SPS[[ 2]] <- read.csv("Ours/N0025_StrongPS_G1.csv")
TempG1.SPS[[ 3]] <- read.csv("Ours/N0050_StrongPS_G1.csv")
TempG1.SPS[[ 4]] <- read.csv("Ours/N0100_StrongPS_G1.csv")
TempG1.SPS[[ 5]] <- read.csv("Ours/N0250_StrongPS_G1.csv")
TempG1.SPS[[ 6]] <- read.csv("Ours/N0500_StrongPS_G1.csv")
TempG2.NPS[[ 2]] <- read.csv("Ours/N0025_NoPS_G2.csv")
TempG2.NPS[[ 3]] <- read.csv("Ours/N0050_NoPS_G2.csv")
TempG2.NPS[[ 4]] <- read.csv("Ours/N0100_NoPS_G2.csv")
TempG2.NPS[[ 5]] <- read.csv("Ours/N0250_NoPS_G2.csv")
TempG2.NPS[[ 6]] <- read.csv("Ours/N0500_NoPS_G2.csv")
TempG2.SPS[[ 2]] <- read.csv("Ours/N0025_StrongPS_G2.csv")
TempG2.SPS[[ 3]] <- read.csv("Ours/N0050_StrongPS_G2.csv")
TempG2.SPS[[ 4]] <- read.csv("Ours/N0100_StrongPS_G2.csv")
TempG2.SPS[[ 5]] <- read.csv("Ours/N0250_StrongPS_G2.csv")
TempG2.SPS[[ 6]] <- read.csv("Ours/N0500_StrongPS_G2.csv")

alpha <- 0.05
LP.NPS.G1 <- LP.NPS.G2 <- LP.SPS.G1 <- LP.SPS.G2 <- list()

for(ii in 2:6){
    
    UDS <- unique( TempG1.NPS[[ii]]$Data.Seed )
    Result.Ours$NPS.G1[[ii]] <- matrix(0,length(UDS),6)
    
    LP.NPS.G1[[ii]] <- rep(0,length(UDS))
    
    for(jj in 1:length(UDS)){
        TD <- TempG1.NPS[[ii]][TempG1.NPS[[ii]]$Data.Seed==UDS[jj],]
        TD <- TD[which(!is.na(apply(TD,1,sum))),]
        if(dim(TD)[1]>=5){
            TD <- TD[1:5,]
        }
        LP.NPS.G1[[ii]][jj] <- dim(TD)[1]
        
        Eff <- apply(TD[,5:6],2,median)
        SE <- sqrt( apply((TD[,5:6]-matrix(Eff,dim(TD)[1],2,byrow=T))^2 +
                              TD[,7:8]^2,2,median) )
        Cover <- as.numeric( Eff-qnorm(1-alpha/2)*SE < ATE & Eff+qnorm(1-alpha/2)*SE > ATE )
        Result.Ours$NPS.G1[[ii]][jj,] <- c(Eff,SE,Cover)
        
    }
    Result.Ours$NPS.G1[[ii]] <- data.frame(Result.Ours$NPS.G1[[ii]])
    colnames(Result.Ours$NPS.G1[[ii]]) <- c("Eff_Est_IPS","Eff_Est_CPS_G",
                                            "SE_Est_IPS","SE_Est_CPS_G",
                                            "Cover_Est_IPS","Cover_Est_CPS_G")
    
    UDS <- unique( TempG2.NPS[[ii]]$Data.Seed )
    Result.Ours$NPS.G2[[ii]] <- matrix(0,length(UDS),6)
    
    LP.NPS.G2[[ii]] <- rep(0,length(UDS))
    
    for(jj in 1:length(UDS)){
        TD <- TempG2.NPS[[ii]][TempG2.NPS[[ii]]$Data.Seed==UDS[jj],]
        TD <- TD[which(!is.na(apply(TD,1,sum))),]
        if(dim(TD)[1]>=5){
            TD <- TD[1:5,]
        }
        LP.NPS.G2[[ii]][jj] <- dim(TD)[1]
        
        Eff <- apply(TD[,5:6],2,median)
        SE <- sqrt( apply((TD[,5:6]-matrix(Eff,dim(TD)[1],2,byrow=T))^2 +
                              TD[,7:8]^2,2,median) )
        Cover <- as.numeric( Eff-qnorm(1-alpha/2)*SE < ATE & Eff+qnorm(1-alpha/2)*SE > ATE )
        Result.Ours$NPS.G2[[ii]][jj,] <- c(Eff,SE,Cover)
        
    }
    Result.Ours$NPS.G2[[ii]] <- data.frame(Result.Ours$NPS.G2[[ii]])
    colnames(Result.Ours$NPS.G2[[ii]]) <- c("Eff_Est_IPS","Eff_Est_CPS_G",
                                            "SE_Est_IPS","SE_Est_CPS_G",
                                            "Cover_Est_IPS","Cover_Est_CPS_G")
    
    UDS <- unique( TempG1.SPS[[ii]]$Data.Seed )
    Result.Ours$SPS.G1[[ii]] <- matrix(0,length(UDS),6)
    
    LP.SPS.G1[[ii]] <- rep(0,length(UDS))
    
    for(jj in 1:length(UDS)){
        TD <- TempG1.SPS[[ii]][TempG1.SPS[[ii]]$Data.Seed==UDS[jj],]
        TD <- TD[which(!is.na(apply(TD,1,sum))),]
        if(dim(TD)[1]>=5){
            TD <- TD[1:5,]
        }
        LP.SPS.G1[[ii]][jj] <- dim(TD)[1]
        
        Eff <- apply(TD[,5:6],2,median)
        SE <- sqrt( apply((TD[,5:6]-matrix(Eff,dim(TD)[1],2,byrow=T))^2 +
                              TD[,7:8]^2,2,median) )
        Cover <- as.numeric( Eff-qnorm(1-alpha/2)*SE < ATE & Eff+qnorm(1-alpha/2)*SE > ATE )
        Result.Ours$SPS.G1[[ii]][jj,] <- c(Eff,SE,Cover)
        
    }
    Result.Ours$SPS.G1[[ii]] <- data.frame(Result.Ours$SPS.G1[[ii]])
    colnames(Result.Ours$SPS.G1[[ii]]) <- c("Eff_Est_IPS","Eff_Est_CPS_G",
                                            "SE_Est_IPS","SE_Est_CPS_G",
                                            "Cover_Est_IPS","Cover_Est_CPS_G")
    
    
    UDS <- unique( TempG2.SPS[[ii]]$Data.Seed )
    Result.Ours$SPS.G2[[ii]] <- matrix(0,length(UDS),6)
    
    LP.SPS.G2[[ii]] <- rep(0,length(UDS))
    
    for(jj in 1:length(UDS)){
        TD <- TempG2.SPS[[ii]][TempG2.SPS[[ii]]$Data.Seed==UDS[jj],]
        TD <- TD[which(!is.na(apply(TD,1,sum))),]
        if(dim(TD)[1]>=5){
            TD <- TD[1:5,]
        }
        LP.SPS.G2[[ii]][jj] <- dim(TD)[1]
        
        Eff <- apply(TD[,5:6],2,median)
        SE <- sqrt( apply((TD[,5:6]-matrix(Eff,dim(TD)[1],2,byrow=T))^2 +
                              TD[,7:8]^2,2,median) )
        Cover <- as.numeric( Eff-qnorm(1-alpha/2)*SE < ATE & Eff+qnorm(1-alpha/2)*SE > ATE )
        Result.Ours$SPS.G2[[ii]][jj,] <- c(Eff,SE,Cover)
        
    }
    Result.Ours$SPS.G2[[ii]] <- data.frame(Result.Ours$SPS.G2[[ii]])
    colnames(Result.Ours$SPS.G2[[ii]]) <- c("Eff_Est_IPS","Eff_Est_CPS_G",
                                            "SE_Est_IPS","SE_Est_CPS_G",
                                            "Cover_Est_IPS","Cover_Est_CPS_G")
    
    
    print(ii)
}


Bias.NPS.G1 <- SE.NPS.G1 <- Cover.NPS.G1 <- matrix(0,2,6)
Bias.NPS.G2 <- SE.NPS.G2 <- Cover.NPS.G2 <- matrix(0,2,6)
Bias.SPS.G1 <- SE.SPS.G1 <- Cover.SPS.G1 <- matrix(0,2,6)
Bias.SPS.G2 <- SE.SPS.G2 <- Cover.SPS.G2 <- matrix(0,2,6)

LP.NPS.G1.F <- LP.NPS.G2.F <- LP.SPS.G1.F <- LP.SPS.G2.F <- list()

for(ii in 2:6){
    LP.NPS.G1.F[[ii]] <- which(LP.NPS.G1[[ii]]==5)[1:200]
    LP.NPS.G2.F[[ii]] <- which(LP.NPS.G2[[ii]]==5)[1:200]
    LP.SPS.G1.F[[ii]] <- which(LP.SPS.G1[[ii]]==5)[1:200]
    LP.SPS.G2.F[[ii]] <- which(LP.SPS.G2[[ii]]==5)[1:200]
}

for(ii in 2:6){
    Bias.NPS.G1[,ii] <- c( apply(Result.Ours$NPS.G1[[ii]][LP.NPS.G1.F[[ii]],],2,mean)[1:2]-ATE )
    SE.NPS.G1[,ii] <- c( apply(Result.Ours$NPS.G1[[ii]][LP.NPS.G1.F[[ii]],],2,sd)[1:2] )
    Cover.NPS.G1[,ii] <- c( apply(Result.Ours$NPS.G1[[ii]][LP.NPS.G1.F[[ii]],],2,mean)[5:6] )
    
    Bias.NPS.G2[,ii] <- c( apply(Result.Ours$NPS.G2[[ii]][LP.NPS.G2.F[[ii]],],2,mean)[1:2]-ATE )
    SE.NPS.G2[,ii] <- c( apply(Result.Ours$NPS.G2[[ii]][LP.NPS.G2.F[[ii]],],2,sd)[1:2] )
    Cover.NPS.G2[,ii] <- c( apply(Result.Ours$NPS.G2[[ii]][LP.NPS.G2.F[[ii]],],2,mean)[5:6] )
    
    Bias.SPS.G1[,ii] <- c( apply(Result.Ours$SPS.G1[[ii]][LP.SPS.G1.F[[ii]],],2,mean)[1:2]-ATE )
    SE.SPS.G1[,ii] <- c( apply(Result.Ours$SPS.G1[[ii]][LP.SPS.G1.F[[ii]],],2,sd)[1:2] )
    Cover.SPS.G1[,ii] <- c( apply(Result.Ours$SPS.G1[[ii]][LP.SPS.G1.F[[ii]],],2,mean)[5:6] )
    
    Bias.SPS.G2[,ii] <- c( apply(Result.Ours$SPS.G2[[ii]][LP.SPS.G2.F[[ii]],],2,mean)[1:2]-ATE )
    SE.SPS.G2[,ii] <- c( apply(Result.Ours$SPS.G2[[ii]][LP.SPS.G2.F[[ii]],],2,sd)[1:2] )
    Cover.SPS.G2[,ii] <- c( apply(Result.Ours$SPS.G2[[ii]][LP.SPS.G2.F[[ii]],],2,mean)[5:6] )
    
}

BBB1 <- rbind(Bias.NPS.G1[1,-1],
              Bias.NPS.G2[1,-1],
              Bias.SPS.G1[1,-1],
              Bias.SPS.G2[1,-1])*100
SSS1 <- rbind(SE.NPS.G1[1,-1],
              SE.NPS.G2[1,-1],
              SE.SPS.G1[1,-1],
              SE.SPS.G2[1,-1])*100
CCC1 <- rbind(Cover.NPS.G1[1,-1],
              Cover.NPS.G2[1,-1],
              Cover.SPS.G1[1,-1],
              Cover.SPS.G2[1,-1])

BBB2 <- rbind(Bias.NPS.G1[2,-1],
              Bias.NPS.G2[2,-1],
              Bias.SPS.G1[2,-1],
              Bias.SPS.G2[2,-1])*100
SSS2 <- rbind(SE.NPS.G1[2,-1],
              SE.NPS.G2[2,-1],
              SE.SPS.G1[2,-1],
              SE.SPS.G2[2,-1])*100
CCC2 <- rbind(Cover.NPS.G1[2,-1],
              Cover.NPS.G2[2,-1],
              Cover.SPS.G1[2,-1],
              Cover.SPS.G2[2,-1])

print(data.frame(cbind(c("\\multirow{4}{*}{$\\overline{\\tau}_w$}","","","",
                         "\\multirow{4}{*}{$\\widehat{\\tau}_w$}","","",""),
                       rep(c("& (0.0,0.5)","& (0.0,1.5)","& (1.5,0.5)","& (1.5,1.5)"),2),
                       rep(c("& (0.05,0.20)","& (0.05,0.69)","& (0.27,0.20)","& (0.27,0.69)"),2),
                       rbind(
                           cbind( sprintf("& %0.2f",BBB1[,1]),
                                  sprintf("& %0.2f",SSS1[,1]),
                                  sprintf("& %0.3f",CCC1[,1]),
                                  sprintf("& %0.2f",BBB1[,2]),
                                  sprintf("& %0.2f",SSS1[,2]),
                                  sprintf("& %0.3f",CCC1[,2]),
                                  sprintf("& %0.2f",BBB1[,3]),
                                  sprintf("& %0.2f",SSS1[,3]),
                                  sprintf("& %0.3f",CCC1[,3]),
                                  sprintf("& %0.2f",BBB1[,4]),
                                  sprintf("& %0.2f",SSS1[,4]),
                                  sprintf("& %0.3f",CCC1[,4]),
                                  sprintf("& %0.2f",BBB1[,5]),
                                  sprintf("& %0.2f",SSS1[,5]),
                                  sprintf("& %0.3f",CCC1[,5]) ),
                           cbind( sprintf("& %0.2f",BBB2[,1]),
                                  sprintf("& %0.2f",SSS2[,1]),
                                  sprintf("& %0.3f",CCC2[,1]),
                                  sprintf("& %0.2f",BBB2[,2]),
                                  sprintf("& %0.2f",SSS2[,2]),
                                  sprintf("& %0.3f",CCC2[,2]),
                                  sprintf("& %0.2f",BBB2[,3]),
                                  sprintf("& %0.2f",SSS2[,3]),
                                  sprintf("& %0.3f",CCC2[,3]),
                                  sprintf("& %0.2f",BBB2[,4]),
                                  sprintf("& %0.2f",SSS2[,4]),
                                  sprintf("& %0.3f",CCC2[,4]),
                                  sprintf("& %0.2f",BBB2[,5]),
                                  sprintf("& %0.2f",SSS2[,5]),
                                  sprintf("& %0.3f",CCC2[,5]) ) ),
                       rep(c("\\\\ \\cline{2-18}","\\\\ \\cline{2-18}","\\\\ \\cline{2-18}","\\\\ \\hline"),2) )),row.names=FALSE)





Result.Others <- list()
Result.Others$NPS.G1 <- Result.Others$SPS.G1 <- Result.Others$NPS.G2 <- Result.Others$SPS.G2 <- list()

Temp.NPS.G1 <- Temp.NPS.G2 <- Temp.SPS.G1 <- Temp.SPS.G2 <- matrix(0,210,21)

for(jj in 1:200){
    Temp.NPS.G1[jj,] <- as.numeric( read.csv(sprintf("Others/Others_N0500_sV_No_B%0.4d.csv",jj))[,(1:21)*4-3] )
    Temp.NPS.G2[jj,] <- as.numeric( read.csv(sprintf("Others/Others_N0500_sV_No_B%0.4d.csv",jj))[,(1:21)*4-2] )
    Temp.SPS.G1[jj,] <- as.numeric( read.csv(sprintf("Others/Others_N0500_sV_Strong_B%0.4d.csv",jj))[,(1:21)*4-3] )
    Temp.SPS.G2[jj,] <- as.numeric( read.csv(sprintf("Others/Others_N0500_sV_Strong_B%0.4d.csv",jj))[,(1:21)*4-2] )
}

colnames(Temp.NPS.G1) <- colnames(Temp.SPS.G1) <- colnames(read.csv(sprintf("Others/Others_N0500_sV_No_B%0.4d.csv",1))[,(1:21)*4-3])
colnames(Temp.NPS.G2) <- colnames(Temp.SPS.G2) <- colnames(read.csv(sprintf("Others/Others_N0500_sV_No_B%0.4d.csv",1))[,(1:21)*4-2])


Bias.NPS.G1.O <- apply( Temp.NPS.G1[LP.NPS.G1.F[[6]],] , 2 , mean )[1:7]-4
Bias.NPS.G2.O <- apply( Temp.NPS.G2[LP.NPS.G2.F[[6]],] , 2 , mean )[1:7]-4
Bias.SPS.G1.O <- apply( Temp.SPS.G1[LP.SPS.G1.F[[6]],] , 2 , mean )[1:7]-4
Bias.SPS.G2.O <- apply( Temp.SPS.G2[LP.SPS.G2.F[[6]],] , 2 , mean )[1:7]-4

SE.NPS.G1.O <- apply( Temp.NPS.G1[LP.NPS.G1.F[[6]],] , 2 , sd )[1:7]
SE.NPS.G2.O <- apply( Temp.NPS.G2[LP.NPS.G2.F[[6]],] , 2 , sd )[1:7]
SE.SPS.G1.O <- apply( Temp.SPS.G1[LP.SPS.G1.F[[6]],] , 2 , sd )[1:7]
SE.SPS.G2.O <- apply( Temp.SPS.G2[LP.SPS.G2.F[[6]],] , 2 , sd )[1:7]

Cover.NPS.G1.O <- apply( Temp.NPS.G1[LP.NPS.G1.F[[6]],] , 2 , mean )[15:21]
Cover.NPS.G2.O <- apply( Temp.NPS.G2[LP.NPS.G2.F[[6]],] , 2 , mean )[15:21]
Cover.SPS.G1.O <- apply( Temp.SPS.G1[LP.SPS.G1.F[[6]],] , 2 , mean )[15:21]
Cover.SPS.G2.O <- apply( Temp.SPS.G2[LP.SPS.G2.F[[6]],] , 2 , mean )[15:21]


print(data.frame(
    cbind( c("GLMM","GEE","R-Learner-B","R-Learner-L","U-Learner-B","U-Learner-L","GRF",
             "$\\overline{\\tau}_w$","$\\widehat{\\tau}_w$") ,
           sprintf(" & %0.2f",c(Bias.NPS.G1.O, Bias.NPS.G1[,6])*100),
           sprintf(" & %0.2f",c(SE.NPS.G1.O,   SE.NPS.G1[,6])*100),
           sprintf(" & %0.3f",c(Cover.NPS.G1.O,Cover.NPS.G1[,6])),
           sprintf(" & %0.2f",c(Bias.NPS.G2.O, Bias.NPS.G2[,6])*100),
           sprintf(" & %0.2f",c(SE.NPS.G2.O,   SE.NPS.G2[,6])*100),
           sprintf(" & %0.3f",c(Cover.NPS.G2.O,Cover.NPS.G2[,6])),
           sprintf(" & %0.2f",c(Bias.SPS.G1.O, Bias.SPS.G1[,6])*100),
           sprintf(" & %0.2f",c(SE.SPS.G1.O,   SE.SPS.G1[,6])*100),
           sprintf(" & %0.3f",c(Cover.SPS.G1.O,Cover.SPS.G1[,6])),
           sprintf(" & %0.2f",c(Bias.SPS.G2.O, Bias.SPS.G2[,6])*100),
           sprintf(" & %0.2f",c(SE.SPS.G2.O,   SE.SPS.G2[,6])*100),
           sprintf(" & %0.3f",c(Cover.SPS.G2.O,Cover.SPS.G2[,6])),
           c(rep("\\\\ \\hline",9)))),row.names=FALSE)

range( c(((c(SE.NPS.G1.O[-c(1,2)],SE.NPS.G1[1,6]))/SE.NPS.G1[2,6]-1)*100,
         ((c(SE.NPS.G2.O[-c(1,2)],SE.NPS.G2[1,6]))/SE.NPS.G2[2,6]-1)*100,
         ((c(SE.SPS.G1.O[-c(1,2)],SE.SPS.G1[1,6]))/SE.SPS.G1[2,6]-1)*100,
         ((c(SE.SPS.G2.O[-c(1,2)],SE.SPS.G2[1,6]))/SE.SPS.G2[2,6]-1)*100) )





