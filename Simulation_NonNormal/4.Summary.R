############################################
# Last update : October 09 2021
# Summary of the Results
############################################

PATH <- getwd()
setwd(PATH)

ATE <- 4

Result.Ours <- list()
Result.Ours$NPS.G1   <- Result.Ours$SPS.G1  <- list()

TempG1.NPS <- read.csv("Ours/N0500_NoPS_G1.csv")
TempG1.SPS <- read.csv("Ours/N0500_StrongPS_G1.csv")

alpha <- 0.05

UDS <- unique( TempG1.NPS$Data.Seed )
LP.NPS.G1 <- rep(0,length(UDS))
Result.Ours$NPS.G1 <- matrix(0,length(UDS),6)

LP.NPS.G1 <- rep(0,length(UDS))

for(jj in 1:length(UDS)){
    TD <- TempG1.NPS[TempG1.NPS$Data.Seed==UDS[jj],]
    TD <- TD[which(!is.na(apply(TD,1,sum))),]
    if(dim(TD)[1]>=5){
        TD <- TD[1:5,]
    }
    LP.NPS.G1[jj] <- dim(TD)[1]
    
    Eff <- apply(TD[,5:6],2,median)
    SE <- sqrt( apply((TD[,5:6]-matrix(Eff,dim(TD)[1],2,byrow=T))^2 +
                          TD[,7:8]^2,2,median) )
    Cover <- as.numeric( Eff-qnorm(1-alpha/2)*SE < ATE & Eff+qnorm(1-alpha/2)*SE > ATE )
    Result.Ours$NPS.G1[jj,] <- c(Eff,SE,Cover)
    
}
Result.Ours$NPS.G1 <- data.frame(Result.Ours$NPS.G1)
colnames(Result.Ours$NPS.G1) <- c("Eff_Est_IPS","Eff_Est_CPS_G",
                                  "SE_Est_IPS","SE_Est_CPS_G",
                                  "Cover_Est_IPS","Cover_Est_CPS_G")


UDS <- unique( TempG1.SPS$Data.Seed )
LP.SPS.G1<- rep(0,length(UDS))
Result.Ours$SPS.G1 <- matrix(0,length(UDS),6)

LP.SPS.G1 <- rep(0,length(UDS))

for(jj in 1:length(UDS)){
    TD <- TempG1.SPS[TempG1.SPS$Data.Seed==UDS[jj],]
    TD <- TD[which(!is.na(apply(TD,1,sum))),]
    if(dim(TD)[1]>=5){
        TD <- TD[1:5,]
    }
    LP.SPS.G1[jj] <- dim(TD)[1]
    
    Eff <- apply(TD[,5:6],2,median)
    SE <- sqrt( apply((TD[,5:6]-matrix(Eff,dim(TD)[1],2,byrow=T))^2 +
                          TD[,7:8]^2,2,median) )
    Cover <- as.numeric( Eff-qnorm(1-alpha/2)*SE < ATE & Eff+qnorm(1-alpha/2)*SE > ATE )
    Result.Ours$SPS.G1[jj,] <- c(Eff,SE,Cover)
    
}
Result.Ours$SPS.G1 <- data.frame(Result.Ours$SPS.G1)
colnames(Result.Ours$SPS.G1) <- c("Eff_Est_IPS","Eff_Est_CPS_G",
                                  "SE_Est_IPS","SE_Est_CPS_G",
                                  "Cover_Est_IPS","Cover_Est_CPS_G")

Bias.NPS.G1 <- SE.NPS.G1 <- Cover.NPS.G1 <- c(0,2)
Bias.SPS.G1 <- SE.SPS.G1 <- Cover.SPS.G1 <- c(0,2)

LP.NPS.G1.F <- which(LP.NPS.G1==5)[1:200]
LP.SPS.G1.F <- which(LP.SPS.G1==5)[1:200]



Bias.NPS.G1 <- c( apply(Result.Ours$NPS.G1[LP.NPS.G1.F,],2,mean)[1:2]-ATE )
SE.NPS.G1 <- c( apply(Result.Ours$NPS.G1[LP.NPS.G1.F,],2,sd)[1:2] )
Cover.NPS.G1 <- c( apply(Result.Ours$NPS.G1[LP.NPS.G1.F,],2,mean)[5:6] )

Bias.SPS.G1 <- c( apply(Result.Ours$SPS.G1[LP.SPS.G1.F,],2,mean)[1:2]-ATE )
SE.SPS.G1 <- c( apply(Result.Ours$SPS.G1[LP.SPS.G1.F,],2,sd)[1:2] )
Cover.SPS.G1 <- c( apply(Result.Ours$SPS.G1[LP.SPS.G1.F,],2,mean)[5:6] )




Result.Others <- list()
Result.Others$NPS.G1 <- Result.Others$SPS.G1 <- list()

TempG1.NPS <- TempG1.SPS <- matrix(0,200,21)

for(jj in 1:200){
    TempG1.NPS[jj,] <- as.numeric( read.csv(sprintf("Others/Others_N0500_sV_No_B%0.4d.csv",jj))[,(1:21)*4-3] )
    TempG1.SPS[jj,] <- as.numeric( read.csv(sprintf("Others/Others_N0500_sV_Strong_B%0.4d.csv",jj))[,(1:21)*4-3] )
    
}

colnames(TempG1.NPS) <- colnames(TempG1.SPS) <- colnames(read.csv(sprintf("Others/Others_N0500_sV_No_B%0.4d.csv",1))[,(1:21)*4-3])


Bias.NPS.G1.O <- apply( TempG1.NPS[LP.NPS.G1.F,] , 2 , mean )[1:7]-4
Bias.SPS.G1.O <- apply( TempG1.SPS[LP.SPS.G1.F,] , 2 , mean )[1:7]-4

SE.NPS.G1.O <- apply( TempG1.NPS[LP.NPS.G1.F,] , 2 , sd )[1:7]
SE.SPS.G1.O <- apply( TempG1.SPS[LP.SPS.G1.F,] , 2 , sd )[1:7]

Cover.NPS.G1.O <- apply( TempG1.NPS[LP.NPS.G1.F,] , 2 , mean )[15:21]
Cover.SPS.G1.O <- apply( TempG1.SPS[LP.SPS.G1.F,] , 2 , mean )[15:21]





print(data.frame(
    cbind( c("GLMM","GEE","R-Learner-B","R-Learner-L","U-Learner-B","U-Learner-L","GRF",
             "$\\overline{\\tau}_w$","$\\widehat{\\tau}_w$") ,
           sprintf(" & %0.2f",c(Bias.NPS.G1.O, Bias.NPS.G1)*100),
           sprintf(" & %0.2f",c(SE.NPS.G1.O,   SE.NPS.G1)*100),
           sprintf(" & %0.3f",c(Cover.NPS.G1.O,Cover.NPS.G1)),
           sprintf(" & %0.2f",c(Bias.SPS.G1.O, Bias.SPS.G1)*100),
           sprintf(" & %0.2f",c(SE.SPS.G1.O,   SE.SPS.G1)*100),
           sprintf(" & %0.3f",c(Cover.SPS.G1.O,Cover.SPS.G1)),
           c(rep("\\\\ \\hline",9)))),row.names=FALSE)

range( (c( c(SE.NPS.G1.O,   SE.NPS.G1)[3:8]/c(SE.NPS.G1.O,   SE.NPS.G1)[9],
           c(SE.SPS.G1.O,   SE.SPS.G1)[3:8]/c(SE.SPS.G1.O,   SE.SPS.G1)[9])-1)*100 )
