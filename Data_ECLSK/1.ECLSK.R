############################################
# Last update : October 09 2021
# Application: ECLS-K
############################################

####################################################
# Cleaning 
# Download data from https://nces.ed.gov/ecls/dataproducts.asp
# Make ECLSK_Kto8_child_STATA.dta
####################################################

library(readstata13)
PATH <- "F:/Dropbox/Chan/Research/2021/ClusterEff_Code_Submit/Data_ECLSK"
setwd(PATH)

# RAWDATA <- read.dta13("ECLSK_Kto8_child_STATA.dta") # <- This takes too long time; recommend to run the following stata .do file
# set maxvar 30000
# use /Users/chanpark/Documents/Dropbox/Chan/Research/2020/ClusterEff/Data/ECLSK/ECLSK_Kto8_child_STATA.dta 
# keep CHILDID PARENTID S1_ID S2_ID S4_ID T1_ID T2_ID T4_ID CREGION KURBAN_R GENDER R1_KAGE C1CMOTOR P1CENTER P1HFAMIL WKPARED WKRACETH WKSESL S2KSCTYP C1R4MTSC C1R4RTSC C2R4MTSC C2R4RTSC C4R4MTSC C4R4RTSC 
# save ${path}/short_ECLSK.dta, replace

RAWDATA <- read.dta13("short_ECLSK.dta") 

ShortData <- RAWDATA[,c("CHILDID","PARENTID","S1_ID","S2_ID","S4_ID","T1_ID","T2_ID","T4_ID",
                        "CREGION","KURBAN_R","GENDER","R1_KAGE","C1CMOTOR","P1CENTER","P1HFAMIL",
                        "WKPARED","WKRACETH","WKSESL","S2KSCTYP","C1R4MTSC","C1R4RTSC","C2R4MTSC","C2R4RTSC","C4R4MTSC","C4R4RTSC")]

Data <- ShortData[!is.na(ShortData$CREGION) & !is.na(ShortData$KURBAN_R) & !is.na(ShortData$S2KSCTYP) ,]
Data <- Data[!is.na(Data$GENDER) & 
               ( !is.na(Data$WKRACETH) & Data$WKRACETH!="NOT ASCERTAINED" & Data$WKRACETH!="NOT APPLICABLE" ) & 
               !is.na(Data$R1_KAGE) & (Data$R1_KAGE>0) & 
               !is.na(Data$P1HFAMIL) & 
               ( !is.na(Data$C1CMOTOR) & Data$C1CMOTOR>0 ) & !is.na(Data$WKSESL) & !is.na(Data$WKPARED) &
               ( !is.na(Data$P1CENTER) & Data$P1CENTER!="NOT ASCERTAINED"& Data$P1CENTER!="NOT APPLICABLE" ) ,]

Data <- Data[ ( !is.na(Data$C1R4RTSC) & Data$C1R4RTSC>0 ) , ]                ## READING score (1st semester)
Data <- Data[ Data$S1_ID!="", ]

My.Data <- Data[,c("CHILDID","PARENTID","S1_ID","S2_ID","S4_ID","T1_ID","T2_ID","T4_ID","CREGION","KURBAN_R","S2KSCTYP","GENDER","R1_KAGE","WKRACETH",
                   "C1CMOTOR","P1HFAMIL","WKPARED","WKSESL","P1CENTER","C1R4MTSC","C1R4RTSC","C2R4MTSC","C2R4RTSC","C4R4MTSC","C4R4RTSC")]
colnames(My.Data) <- c("Child.ID","Parent.ID","Kind.ID1","Kind.ID2","Kind.ID4","Teacher.ID1","Teacher.ID2","Teacher.ID4",
                       "C_region","C_location","C_public","S_sex","S_age","S_race",
                       "S_motor","S_familytype","S_parentaledu","S_sestatus","A","Y1","Y2","Y3","Y4","Y5","Y6")

My.Data$C_region_NE <- as.numeric( My.Data$C_region=="NORTHEAST" )
My.Data$C_region_SOUTH <- as.numeric( My.Data$C_region=="SOUTH" )
My.Data$C_region_WEST <- as.numeric( My.Data$C_region=="WEST" )

My.Data$C_location_city <- as.numeric( My.Data$C_location=="CENTRAL CITY" )
My.Data$C_location_rural <- as.numeric( My.Data$C_location=="SMALL TOWN AND RURAL" )

My.Data$C_public_public <- as.numeric( My.Data$C_public=="PUBLIC" )

My.Data$C_M <- 0
for(jj in 1:dim(My.Data)[1]){
  My.Data$C_M[jj] <- sum(My.Data$Kind.ID1[jj]==My.Data$Kind.ID1)
}

My.Data <- My.Data[My.Data$C_M<=25,]

My.Data$S_sex_male <- as.numeric( My.Data$S_sex=="MALE" )

My.Data$S_race_H <- as.numeric( My.Data$S_race=="HISPANIC, RACE SPECIFIED" | My.Data$S_race=="HISPANIC, RACE NOT SPECIFIED" )
My.Data$S_race_B <- as.numeric( My.Data$S_race=="BLACK OR AFRICAN AMERICAN, NON-HISPANIC" )
My.Data$S_race_A <- as.numeric( My.Data$S_race=="ASIAN" )
My.Data$S_race_W <- as.numeric( My.Data$S_race=="WHITE, NON-HISPANIC" )

My.Data$S_familytype_bothparent <- as.numeric( My.Data$S_familytype=="2 PARENTS PLUS SIBLINGS" | My.Data$S_familytype=="2 PARENTS NO SIBLING" )

My.Data$S_parentaledu_college <- as.numeric( My.Data$S_parentaledu=="DOCTORATE OR PROFESSIONAL DEGREE" | 
                                               My.Data$S_parentaledu=="BACHELOR'S DEGREE" |
                                               My.Data$S_parentaledu=="MASTER'S DEGREE (MA, MS)" |
                                               My.Data$S_parentaledu=="SOME COLLEGE")
My.Data$S_age <- as.numeric(My.Data$S_age)
My.Data$S_motor <- as.numeric(My.Data$S_motor)
My.Data$S_sestatus <- as.numeric(My.Data$S_sestatus)
My.Data$A <- as.numeric(My.Data$A=="YES")
My.Data$Y <- as.numeric(My.Data$Y2)
My.Data$GP <- My.Data$Kind.ID1


indmat <- cbind(1:length(unique( My.Data$Kind.ID )),sort(unique( My.Data$Kind.ID )))
for(ii in 1:dim(My.Data)[1]){
  My.Data$GP[ii] <- indmat[which(indmat[,2]==My.Data$GP[ii]),1]
}
My.Data$Kind.ID <- as.numeric(My.Data$Kind.ID)

Final.Data <- My.Data[,c("Y", "A", "GP", 
                         "C_M",
                         "C_region_NE","C_region_SOUTH","C_region_WEST",
                         "C_location_city","C_location_rural",
                         "C_public_public",
                         "S_sex_male","S_age","S_race_A","S_race_B","S_race_H","S_race_W",
                         "S_motor","S_familytype_bothparent","S_parentaledu_college","S_sestatus")]
write.csv(Final.Data,"Reading1_ECLSK.csv",row.names=FALSE)


####################################################
# Nuisance Function Estimation
####################################################

source("../MySL.R")
source("../ClusterFtSource.R")

Para.Mat <- expand.grid(1:100,1:2)

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
SL.hpara$MTRY <- c(5,10)                  # random forest parameters
SL.hpara$MLPL <- c(2,4)                   # number of nodes in MLP
SL.hpara$NMN <- 50                        # gbm parameter
SL.hpara$MLPdecay <- 10^c(-4,-5)          # MLP decay parameter

for(BATCH in 1:200){
  
  NF.num <- Para.Mat[BATCH,1]
  TYPE <- Para.Mat[BATCH,2]
  DataType <- "Reading"
  
  pos.X <- 3:19
  pos.AX <- 2:29
  CPS.B <- 5
  
  
  Data <- read.csv(sprintf("%s1_ECLSK.csv",DataType))
  Data <- Data[order(Data$GP),]
  TrueData <- Data[,-which(colnames(Data)=="GP")]
  TrueData$GP <- Data[,which(colnames(Data)=="GP")]
  TrueData$C_M <- scale(TrueData$C_M)
  TrueData$S_age <- scale(TrueData$S_age)
  TrueData$S_motor <- scale(TrueData$S_motor)
  TrueData$S_sestatus <- scale(TrueData$S_sestatus)
  
  TrueData0 <- TrueData
  TrueData0[,2] <- 0
  TrueData1 <- TrueData
  TrueData1[,2] <- 1
  
  N <- length( unique(TrueData$GP) )
  M <- as.numeric(table(TrueData$GP))
  
  TrueDataA <- TrueData
  TrueDataA$A.others <- ( rep(aggregate(A~GP,TrueData,FUN="sum")[,2],M)-TrueData$A )/(rep(M,M)-1)
  TrueDataA$S_sex_male.others <- ( rep(aggregate(S_sex_male~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_sex_male )/(rep(M,M)-1)
  TrueDataA$S_age.others <- ( rep(aggregate(S_age~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_age )/(rep(M,M)-1)
  TrueDataA$S_race_A.others <- ( rep(aggregate(S_race_A~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_A )/(rep(M,M)-1)
  TrueDataA$S_race_B.others <- ( rep(aggregate(S_race_B~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_B )/(rep(M,M)-1)
  TrueDataA$S_race_H.others <- ( rep(aggregate(S_race_H~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_H )/(rep(M,M)-1)
  TrueDataA$S_race_W.others <- ( rep(aggregate(S_race_W~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_W )/(rep(M,M)-1)
  TrueDataA$S_motor.others <- ( rep(aggregate(S_motor~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_motor )/(rep(M,M)-1)
  TrueDataA$S_familytype_bothparent.others <- ( rep(aggregate(S_familytype_bothparent~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_familytype_bothparent )/(rep(M,M)-1)
  TrueDataA$S_parentaledu_college.others <- ( rep(aggregate(S_parentaledu_college~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_parentaledu_college )/(rep(M,M)-1)
  TrueDataA$S_sestatus.others <- ( rep(aggregate(S_sestatus~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_sestatus )/(rep(M,M)-1)
  
  TrueDataA <- TrueDataA[,c(colnames(TrueData)[2:19],paste(colnames(TrueData)[c(2,10:19)],".others",sep=""),"GP")]
  
  set.seed(NF.num*1000)
  
  index <- list()
  for(ii in 1:N){
    index[[ii]] <- which( TrueData$GP==unique(TrueData$GP)[ii] )
  }
  
  
  SS <- list()
  SS$SS1 <- SS$SS2 <- SS$SS3 <- NULL
  SS$MS <- SS$AS <- NULL
  Cindex <- list()
  ttt <- list()
  RI <- rep(0,length(unique(M)))
  for(mm in sort(unique(M))){
    Cindex[[mm]] <- which(M==mm)
    ttt[[mm]] <- sample(Cindex[[mm]],length(Cindex[[mm]]))
    tms <- ttt[[mm]][1:round(length(ttt[[mm]])/2)]
    tas <- ttt[[mm]][(round(length(ttt[[mm]])/2)+1):length(ttt[[mm]])]
    
    if(runif(1)<0.5){
      SS$MS <- c(SS$MS,tms) ; SS$AS <- c(SS$AS,tas)
    } else {
      SS$MS <- c(SS$MS,tas) ; SS$AS <- c(SS$AS,tms)
    }
    
    print(mm)
  }
  
  SS$MS  <- sort(SS$MS)
  SS$AS  <- sort(SS$AS)
  
  SSI <- list()
  SSI$MS <- NULL
  for(ii in 1:length(SS$MS)){
    SSI$MS <- c(SSI$MS, index[[ SS$MS[ii] ]])
  }
  SSI$AS <- NULL
  for(ii in 1:length(SS$AS)){
    SSI$AS <- c(SSI$AS, index[[ SS$AS[ii] ]])
  }
  
  
  
  if(TYPE==1){
    
    SSIcsv <- cbind(rep(0,sum(M)))
    SSIcsv[SSI$MS,1] <- 1
    SSIcsv[SSI$AS,1] <- 2
    
    SSIcsv <- data.frame(SSIcsv)
    colnames(SSIcsv) <- c("SS2")
    
    write.csv(SSIcsv,sprintf("%s1_SSlist_NF%0.4d.csv",DataType,NF.num),row.names=FALSE)
    
    OR.fit <- list()
    OR.fit$MS <- OR.Estimation(TrueData$Y[SSI$MS],TrueData$A[SSI$MS],TrueData[SSI$MS,pos.X],TrueData$GP[SSI$MS],type="ML",outcome.type=gaussian(),SL.hpara=SL.hpara)
    OR.fit$AS <- OR.Estimation(TrueData$Y[SSI$AS],TrueData$A[SSI$AS],TrueData[SSI$AS,pos.X],TrueData$GP[SSI$AS],type="ML",outcome.type=gaussian(),SL.hpara=SL.hpara)
    
    OR.est <- list()
    OR.est$ML2 <- rep(0,sum(M))
    OR.est$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,TrueData$A[SSI$AS],TrueData[SSI$AS,pos.X])$output
    OR.est$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,TrueData$A[SSI$MS],TrueData[SSI$MS,pos.X])$output
    
    OR.est1 <- list()
    OR.est1$ML2 <- rep(0,sum(M))
    OR.est1$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,TrueData1$A[SSI$AS],TrueData1[SSI$AS,pos.X])$output
    OR.est1$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,TrueData1$A[SSI$MS],TrueData1[SSI$MS,pos.X])$output
    
    
    OR.est0 <- list()
    OR.est0$ML2 <- rep(0,sum(M))
    OR.est0$ML2[SSI$AS] <- OR.Prediction(OR.fit$MS,TrueData0$A[SSI$AS],TrueData0[SSI$AS,pos.X])$output
    OR.est0$ML2[SSI$MS] <- OR.Prediction(OR.fit$AS,TrueData0$A[SSI$MS],TrueData0[SSI$MS,pos.X])$output
    
    RES <- cbind(rep((1:N),M),(1:sum(M)),rep(BATCH,sum(M)),OR.est$ML2,OR.est1$ML2,OR.est0$ML2)
    RES <- data.frame(RES)
    colnames(RES) <- c("Cluster.ID","Ind.ID","BATCH","OR","OR1","OR0")
    
    write.csv(RES,sprintf("Reading/%s1_OR.ML_NF%0.4d.csv",DataType,NF.num),row.names=FALSE)
    
    PS.fit <- list()
    PS.fit$MS <- IPS.Estimation(TrueData$A[SSI$MS],TrueData[SSI$MS,pos.X],TrueData$GP[SSI$MS],type="ML",SL.hpara=SL.hpara)
    PS.fit$AS <- IPS.Estimation(TrueData$A[SSI$AS],TrueData[SSI$AS,pos.X],TrueData$GP[SSI$AS],type="ML",SL.hpara=SL.hpara)
    
    IPS.est <- list()
    IPS.est$ML2 <- rep(0,sum(M))
    IPS.est$ML2[SSI$AS] <- IPS.Prediction(PS.fit$MS,TrueData$A[SSI$AS],TrueData[SSI$AS,pos.X],TrueData$GP[SSI$AS])$IPS
    IPS.est$ML2[SSI$MS] <- IPS.Prediction(PS.fit$AS,TrueData$A[SSI$MS],TrueData[SSI$MS,pos.X],TrueData$GP[SSI$MS])$IPS
    
    RES <- cbind(rep((1:N),M),(1:sum(M)),rep(BATCH,sum(M)),IPS.est$ML2)
    RES <- data.frame(RES)
    colnames(RES) <- c("Cluster.ID","Ind.ID","BATCH","IPS")
    
    write.csv(RES,sprintf("Reading/%s1_IPS.ML_NF%0.4d.csv",DataType,NF.num),row.names=FALSE)
    
  } else if(TYPE==2) {
    
    PS.fit <- list()
    IPS.est <- list()
    
    IPS.est$CPS <- rep(0,sum(M))
    CPS.temp.mat <- matrix(0,sum(M),CPS.B)
    
    for(cps.iter in 1:CPS.B){
      
      SSI.temp <- list()
      SSI.temp$MS <- SSI.temp$AS <- NULL
      for(ii in 1:N){
        pos <- which(rep(1:N, M)==ii)
        if(length(pos)>1){
          pos <- sample(pos,2)
          if(sum(SS$MS==ii)>0){ 
            SSI.temp$MS <- c(SSI.temp$MS,pos) 
          } else { 
            SSI.temp$AS <- c(SSI.temp$AS,pos) 
          }
        } else {
          if(sum(SS$MS==ii)>0){ 
            SSI.temp$MS <- c(SSI.temp$MS,NA) 
          } else { 
            SSI.temp$AS <- c(SSI.temp$AS,NA) 
          }
        }
      }
      
      SSI.temp$MS <- sort(SSI.temp$MS)
      SSI.temp$AS <- sort(SSI.temp$AS)
      
      DataB <- TrueDataA
      SL.hpara2 <- SL.hpara
      
      CVL.MS <- CVL.AS <- list()
      
      
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
      
      
      PS.fit$CPS.MS <- IPS.Estimation(DataB$A[SSI.temp$MS],DataB[SSI.temp$MS,pos.AX],DataB$GP[SSI.temp$MS],type="ML",SL.hpara=SL.hpara2,CVlist=CVL.MS)
      PS.fit$CPS.AS <- IPS.Estimation(DataB$A[SSI.temp$AS],DataB[SSI.temp$AS,pos.AX],DataB$GP[SSI.temp$AS],type="ML",SL.hpara=SL.hpara2,CVlist=CVL.AS)
      
      SSI$MS.A <- setdiff(SSI$MS,which(rep(M,M)==1))
      SSI$AS.A <- setdiff(SSI$AS,which(rep(M,M)==1))
      
      CPS.temp.mat[SSI$AS.A,cps.iter] <- IPS.Prediction(PS.fit$CPS.MS,TrueDataA$A[SSI$AS.A],TrueDataA[SSI$AS.A,pos.AX],TrueDataA$GP[SSI$AS.A])$IPS
      CPS.temp.mat[SSI$MS.A,cps.iter] <- IPS.Prediction(PS.fit$CPS.AS,TrueDataA$A[SSI$MS.A],TrueDataA[SSI$MS.A,pos.AX],TrueDataA$GP[SSI$MS.A])$IPS
      
      print(cps.iter)
    }
    
    IPS.est$CPS <- apply(CPS.temp.mat,1,median)
    
    IPS.est$CPS[which(rep(M,M)==1)] <- 1000
    
    RES <- cbind(rep((1:N),M),(1:sum(M)),rep(BATCH,sum(M)),IPS.est$CPS)
    RES <- data.frame(RES)
    colnames(RES) <- c("Cluster.ID","Ind.ID","BATCH","CPS")
    
    write.csv(RES,sprintf("Reading/%s1_CPS.ML_NF%0.4d.csv",DataType,NF.num),row.names=FALSE)
    
  } 
  
  
}

####################################################
# Calculate IF
####################################################

DataType <- "Reading"

pos.X <- 3:19
pos.AX <- 2:29
CPS.B <- 10

Data <- read.csv(sprintf("%s1_ECLSK.csv",DataType))
Data <- Data[order(Data$GP),]
TrueData <- Data[,-which(colnames(Data)=="GP")]
TrueData$GP <- Data[,which(colnames(Data)=="GP")]
TrueData$C_M <- scale(TrueData$C_M)
TrueData$S_age <- scale(TrueData$S_age)
TrueData$S_motor <- scale(TrueData$S_motor)
TrueData$S_sestatus <- scale(TrueData$S_sestatus)

TrueData0 <- TrueData
TrueData0[,2] <- 0
TrueData1 <- TrueData
TrueData1[,2] <- 1

N <- length( unique(TrueData$GP) )
M <- as.numeric(table(TrueData$GP))

TrueDataA <- TrueData
TrueDataA$A.others <- ( rep(aggregate(A~GP,TrueData,FUN="sum")[,2],M)-TrueData$A )/(rep(M,M)-1)
TrueDataA$S_sex_male.others <- ( rep(aggregate(S_sex_male~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_sex_male )/(rep(M,M)-1)
TrueDataA$S_age.others <- ( rep(aggregate(S_age~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_age )/(rep(M,M)-1)
TrueDataA$S_race_A.others <- ( rep(aggregate(S_race_A~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_A )/(rep(M,M)-1)
TrueDataA$S_race_B.others <- ( rep(aggregate(S_race_B~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_B )/(rep(M,M)-1)
TrueDataA$S_race_H.others <- ( rep(aggregate(S_race_H~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_H )/(rep(M,M)-1)
TrueDataA$S_race_W.others <- ( rep(aggregate(S_race_W~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_race_W )/(rep(M,M)-1)
TrueDataA$S_motor.others <- ( rep(aggregate(S_motor~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_motor )/(rep(M,M)-1)
TrueDataA$S_familytype_bothparent.others <- ( rep(aggregate(S_familytype_bothparent~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_familytype_bothparent )/(rep(M,M)-1)
TrueDataA$S_parentaledu_college.others <- ( rep(aggregate(S_parentaledu_college~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_parentaledu_college )/(rep(M,M)-1)
TrueDataA$S_sestatus.others <- ( rep(aggregate(S_sestatus~GP,TrueData,FUN="sum")[,2],M)-TrueData$S_sestatus )/(rep(M,M)-1)

TrueDataA <- TrueDataA[,c(colnames(TrueData)[2:19],paste(colnames(TrueData)[c(2,10:19)],".others",sep=""),"GP")]




index <- list()
for(ii in 1:N){
  index[[ii]] <- which( TrueData$GP==unique(TrueData$GP)[ii] )
}

SL.hpara <- list()
SL.hpara$SLL <- 1:11
SL.hpara$MTRY <- c(3,6,9)
SL.hpara$MLPL <- c(1,2,3,4)
SL.hpara$NMN <- round(0.05*N)

Cov.Pos <- 3:19
Cov.Pos.A <- 2:19
Clu.index <- 1:25

OR.fit <- PS.fit <- list()
OR.est <- OR.est1 <- OR.est0 <- list()
IPS.est <- GPS.est <- OPS.est <- list()
IF <- list()

#################################
# Parametric Model
#################################


OR.fit$glmer.noS <- OR.Estimation(TrueData$Y,TrueData$A,TrueData[,Cov.Pos],TrueData$GP,type="GLMM",outcome.type=gaussian(),SL.hpara=SL.hpara)
OR.fit$gee.noS <-   OR.Estimation(TrueData$Y,TrueData$A,TrueData[,Cov.Pos],TrueData$GP,type="GEE",outcome.type=gaussian(),SL.hpara=SL.hpara)
OR.fit$glm.noS <-   OR.Estimation(TrueData$Y,TrueData$A,TrueData[,Cov.Pos],TrueData$GP,type="GLM",outcome.type=gaussian(),SL.hpara=SL.hpara)



OR.est$glmer.noS <- OR.Prediction(OR.fit$glmer.noS,TrueData$A,TrueData[,Cov.Pos])$output
OR.est$gee.noS <- OR.Prediction(OR.fit$gee.noS,TrueData$A,TrueData[,Cov.Pos])$output
OR.est$glm.noS <- OR.Prediction(OR.fit$glm.noS,TrueData$A,TrueData[,Cov.Pos])$output
OR.est1$glmer.noS <- OR.Prediction(OR.fit$glmer.noS,TrueData1$A,TrueData1[,Cov.Pos])$output
OR.est1$gee.noS <- OR.Prediction(OR.fit$gee.noS,TrueData1$A,TrueData1[,Cov.Pos])$output
OR.est1$glm.noS <- OR.Prediction(OR.fit$glm.noS,TrueData1$A,TrueData1[,Cov.Pos])$output
OR.est0$glmer.noS <- OR.Prediction(OR.fit$glmer.noS,TrueData0$A,TrueData0[,Cov.Pos])$output
OR.est0$gee.noS <- OR.Prediction(OR.fit$gee.noS,TrueData0$A,TrueData0[,Cov.Pos])$output
OR.est0$glm.noS <- OR.Prediction(OR.fit$glm.noS,TrueData0$A,TrueData0[,Cov.Pos])$output

PS.fit$glmer.noS <- IPS.Estimation(TrueData$A,TrueData[,Cov.Pos],TrueData$GP,type="GLMM",SL.hpara=SL.hpara)
PS.fit$gee.noS <- IPS.Estimation(TrueData$A,TrueData[,Cov.Pos],TrueData$GP,type="GEE",SL.hpara=SL.hpara)
PS.fit$glm.noS <- IPS.Estimation(TrueData$A,TrueData[,Cov.Pos],TrueData$GP,type="GLM",SL.hpara=SL.hpara)

IPS.est$glmer.noS <- IPS.Prediction(PS.fit$glmer.noS,TrueData$A,TrueData[,Cov.Pos],TrueData$GP)$IPS
IPS.est$gee.noS <- IPS.Prediction(PS.fit$gee.noS,TrueData$A,TrueData[,Cov.Pos],TrueData$GP)$IPS
IPS.est$glm.noS <- IPS.Prediction(PS.fit$glm.noS,TrueData$A,TrueData[,Cov.Pos],TrueData$GP)$IPS

IF$glmer1 <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$GP,OR.est$glmer.noS,OR.est1$glmer.noS,OR.est0$glmer.noS,IPS.est$glmer.noS)
IF$gee1 <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$GP,OR.est$gee.noS,OR.est1$gee.noS,OR.est0$gee.noS,IPS.est$gee.noS)
IF$glm1 <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$GP,OR.est$glm.noS,OR.est1$glm.noS,OR.est0$glm.noS,IPS.est$glm.noS)

ResY <- (TrueData$Y - OR.est$glmer.noS)
Res.Gp.Y <- aggregate(ResY~TrueData$GP,FUN="sum")[,2]

MSB <- (sum(Res.Gp.Y^2/M) - sum(Res.Gp.Y)^2/sum(M))/(N-1)
MSW <- (sum(ResY^2)-sum(Res.Gp.Y^2/M))/(sum(M)-N)
N.m <- (sum(M) - sum(M^2)/sum(M))/(N-1)
ICCY <- (MSB-MSW)/(MSB+(N.m-1)*MSW)

ResA <- (TrueData$A - IPS.est$glmer.noS)
Res.Gp.A <- aggregate(ResA~TrueData$GP,FUN="sum")[,2]

MSB <- (sum(Res.Gp.A^2/M) - sum(Res.Gp.A)^2/sum(M))/(N-1)
MSW <- (sum(ResA^2)-sum(Res.Gp.A^2/M))/(sum(M)-N)
N.m <- (sum(M) - sum(M^2)/sum(M))/(N-1)
ICCA <- (MSB-MSW)/(MSB+(N.m-1)*MSW)

ICCY ; ICCA

#################################
# Our Method
#################################

RES <- matrix(0,100,10)

for(NF.num in 1:100){
  
  SSIcsv <- read.csv(sprintf("%s_u/%s1_SSlist_NF%0.4d.csv",DataType,DataType,NF.num))
  
  SS <- list()
  SS$MS <- which( unique( cbind(rep(1:N,M),SSIcsv$SS2) )[,2]==1 )
  SS$AS <- which( unique( cbind(rep(1:N,M),SSIcsv$SS2) )[,2]==2 )
  
  SSI <- list()
  SSI$MS <- which( SSIcsv$SS2==1 )
  SSI$AS <- which( SSIcsv$SS2==2 )
  
  SSI$MS.A <- setdiff(SSI$MS,which(rep(M,M)==1))
  SSI$AS.A <- setdiff(SSI$AS,which(rep(M,M)==1))
  
  OR2 <- read.csv(sprintf("%s_u/%s1_OR.ML_NF%0.4d.csv",DataType,DataType,NF.num))
  OR.est$ML2 <- OR2$OR
  OR.est1$ML2 <- OR2$OR1
  OR.est0$ML2 <- OR2$OR0
  
  PS2 <- read.csv(sprintf("%s_u/%s1_IPS.ML_NF%0.4d.csv",DataType,DataType,NF.num))
  IPS.est$ML2 <- PS2$IPS
  
  PSC <- read.csv(sprintf("%s_u/%s1_CPS.ML_NF%0.4d.csv",DataType,DataType,NF.num))
  IPS.est$CPS <- PSC$CPS
  
  IPS.est$CPS[rep(M,M)==1] <- IPS.est$ML2[rep(M,M)==1]
  
  IF$IPS2.B <- findIFwithoutB(TrueData$Y,TrueData$A,TrueData$GP,OR.est$ML2,OR.est1$ML2,OR.est0$ML2,IPS.est$ML2)
  
  Clustering <- rep(M,M)-1
  
  Beta.MS <- findbetaDisc(TrueData$Y[SSI$MS],TrueData$A[SSI$MS],TrueData$GP[SSI$MS],OR.est$ML2[SSI$MS],IPS.est$CPS[SSI$MS],Clustering=Clustering[SSI$MS])$beta.nobdd
  Beta.AS <- findbetaDisc(TrueData$Y[SSI$AS],TrueData$A[SSI$AS],TrueData$GP[SSI$AS],OR.est$ML2[SSI$AS],IPS.est$CPS[SSI$AS],Clustering=Clustering[SSI$AS])$beta.nobdd
  
  IF$CPS.BnoB[SS$MS] <- findIFwithB.Disc(Data$Y[SSI$MS],Data$A[SSI$MS],Data$GP[SSI$MS],OR.est$ML2[SSI$MS],OR.est1$ML2[SSI$MS],OR.est0$ML2[SSI$MS],IPS.est$CPS[SSI$MS],Beta.MS,Clustering=Clustering[SSI$MS])
  IF$CPS.BnoB[SS$AS] <- findIFwithB.Disc(Data$Y[SSI$AS],Data$A[SSI$AS],Data$GP[SSI$AS],OR.est$ML2[SSI$AS],OR.est1$ML2[SSI$AS],OR.est0$ML2[SSI$AS],IPS.est$CPS[SSI$AS],Beta.AS,Clustering=Clustering[SSI$AS])
  
  
  RES[NF.num,] <- c(sapply(IF,mean),
                    sqrt(var(IF$glmer1)/N),
                    sqrt(var(IF$gee1)/N),
                    sqrt(var(IF$glm1)/N),
                    sqrt(obtainV(IF$IPS2.B,split=2)/N),
                    sqrt(obtainV(IF$CPS.BnoB,split=2)/N))
  
  print(NF.num)
}

R1 <- RES
Est.R1 <- apply(R1,2,median)[1:5]
SE.R1 <- sqrt(apply(R1[,6:10]^2 + (R1[,1:5] - matrix(Est.R1,dim(R1),5,byrow=T))^2,2,median))
Est.R1  ; SE.R1


RESULT <- data.frame(RES)
colnames(RESULT) <- c(paste("Est_",c("GLMM","GEE","GLM","ML","CPS_noB"),sep=""),
                      paste("SE_",c("GLMM","GEE","GLM","ML","CPS_noB"),sep=""))

write.csv(RESULT,sprintf("%s_Ours.csv",DataType),row.names=FALSE)



######################
# R-learner
######################

library(rlearner)
library(grf)
library(foreach)
library(doParallel)

Data <- read.csv("Reading1_ECLSK.csv")
Data <- Data[order(Data$GP),]
TrueData <- Data[,-which(colnames(Data)=="GP")]
TrueData$GP <- Data[,which(colnames(Data)=="GP")]

TrueData0 <- TrueData
TrueData0[,2] <- 0
TrueData1 <- TrueData
TrueData1[,2] <- 1

N <- length( unique(TrueData$GP) )
M <- as.numeric(table(TrueData$GP))

index <- list()
for(ii in 1:N){
  index[[ii]] <- which(Data$GP==unique(Data$GP)[ii])
}

RRR <- matrix(0,100,8)
RRR <- data.frame(RRR)
colnames(RRR) <- c("Est_R_Lasso",
                   "Est_R_Boost",
                   "Est_U_Lasso",
                   "Est_U_Boost",
                   "SE_R_Lasso",
                   "SE_R_Boost",
                   "SE_U_Lasso",
                   "SE_U_Boost")

for(ii in 1:100){
  
  set.seed(100*ii + ii)
  
  RL <- rlasso(as.matrix(TrueData[,-c(1,2,20)]), TrueData[,2], TrueData[,1])
  
  output11 <- ATEcalculate(TrueData, RL, weight.vector=1/TrueData[,3], clusters=TrueData[,20])
  
  RLest <- c( output11[1])
  RLsd <- c( output11[2] )
  
  RB <- rboost(as.matrix(TrueData[,-c(1,2,20)]), TrueData[,2], TrueData[,1])
  
  output21 <- ATEcalculate(TrueData, RB, weight.vector=1/TrueData[,3], clusters=TrueData[,20])
  
  RBest <- c( output21[1])
  RBsd <- c( output21[2] )
  
  
  
  UL <- ulasso(as.matrix(TrueData[,-c(1,2,20)]), TrueData[,2], TrueData[,1])
  
  output31 <- ATEcalculate(TrueData, UL, weight.vector=1/TrueData[,3], clusters=TrueData[,20])
  
  ULest <- c( output31[1])
  ULsd <- c( output31[2] )
  
  UB <- uboost(as.matrix(TrueData[,-c(1,2,20)]), TrueData[,2], TrueData[,1])
  
  output41 <- ATEcalculate(TrueData, UB, weight.vector=1/TrueData[,3], clusters=TrueData[,20])
  
  UBest <- c( output41[1])
  UBsd <- c( output41[2] )
  
  test <- c(RLest,RBest,ULest,UBest,
            RLsd, RBsd, ULsd, UBsd )
  test <- data.frame(matrix(test,1,8))
  RRR[ii,] <- c(RLest,RBest,ULest,UBest,RLsd, RBsd, ULsd, UBsd )
  
}

write.csv(RRR,"Reading_RU.csv",row.names=FALSE)



######################
# GRF
######################

GRFest <- matrix(0,100,2)
library(grf)
for(ii in 1:100){
    set.seed(100*ii + ii)
    CF <- causal_forest(X=TrueData[,-c(1,2,20)],Y=TrueData[,1],W=TrueData[,2],clusters=TrueData[,20],sample.weights=1/rep(M,M),num.trees=5000)
    GRFest[ii,] <- average_treatment_effect(CF)
    print(ii)
}

write.csv(GRFest,"Reading_GRF.csv",row.names=FALSE)



####################################################
# Summary
####################################################

R1 <- read.csv("Reading_Ours.csv")
Est.R1 <- apply(R1,2,median)[1:5]
SE.R1 <- sqrt(apply(R1[,6:10]^2 + (R1[,1:5] - matrix(Est.R1,dim(R1),5,byrow=T))^2,2,median))

R2 <- read.csv("Reading_RU.csv")
Est.R2 <- apply(R2,2,median)[1:4]
SE.R2 <- sqrt(apply(R2[,5:8]^2 + (R2[,1:4] - matrix(Est.R2,dim(R2),4,byrow=T))^2,2,median))

R3 <- read.csv("Reading_GRF.csv")
Est.R3 <- apply(R3,2,median)[1]
SE.R3 <- matrix(c(sqrt(median(R3[,2]^2 + (R3[,1]-Est.R3)^2))),1,1)

  
  
  
  print(data.frame(cbind(c("Method","Estimate","SE","RE($\\widehat{\\tau}_{\\rm ATE},\\cdot$)"),
                         rbind( c(" & GLMM"," & GEE"," & GLM",
                                  " & R-Learner-B"," & R-Learner-L",
                                  " & U-Learner-B"," & U-Learner-L", "GRF",
                                  " & $\\overline{\\tau}_w$"," & $\\widehat{\\tau}_w$"),
                                sprintf(" & %0.3f", c(Est.R1[c(1,2,3)], Est.R2[c(2,1,4,3)], Est.R3, Est.R1[c(4,5)])),
                                sprintf(" & %0.3f", c(SE.R1[c(1,2,3)],  SE.R2[c(2,1,4,3)],  SE.R3,  SE.R1[c(4,5)]) ),
                                sprintf(" & %0.3f", (c(SE.R1[c(1,2,3)],  SE.R2[c(2,1,4,3)], SE.R3,  SE.R1[c(4,5)])/SE.R1[5])^2)),
                         rep("\\\\ \\hline"))),row.names=FALSE)
  




####################################################
# Balance/Overlap
####################################################

Data <- read.csv("Reading1_ECLSK.csv")
M <- Data$C_M
N <- length(Data$GP)

IPS <- CPS <- list()
for(ii in 1:100){
  IPS[[ii]] <- read.csv(sprintf("Reading_u/Reading1_IPS.ML_NF%0.4d.csv",ii)) # folder is empty in the zip file due to its size
  CPS[[ii]] <- read.csv(sprintf("Reading_u/Reading1_CPS.ML_NF%0.4d.csv",ii)) # folder is empty in the zip file due to its size
}

IPS.Med <- array(unlist(IPS), c( dim(IPS[[1]])[1], dim(IPS[[1]])[2], 100))
CPS.Med <- array(unlist(CPS), c( dim(CPS[[1]])[1], dim(CPS[[1]])[2], 100))
IPS.Med <- data.frame( apply(IPS.Med,c(1:2),median) )
CPS.Med <- data.frame( apply(CPS.Med,c(1:2),median) )
CPS.Med[CPS.Med[,4]>1,4] <- IPS.Med[CPS.Med[,4]>1,4] 

No.Adj.Tstat <- IPS.Adj.Tstat <- CPS.Adj.Tstat <- matrix(0,17,1000)
IND <- cbind(cumsum(M)-M+1,cumsum(M))
V1 <- Data[,4:20]*Data$A/mean(Data$A) - Data[,4:20]*(1-Data$A)/mean(1-Data$A)
V2 <- Data[,4:20]*Data$A/(IPS.Med[,4]) - Data[,4:20]*(1-Data$A)/(IPS.Med[,4])
V3 <- Data[,4:20]*Data$A/(CPS.Med[,4]) - Data[,4:20]*(1-Data$A)/(CPS.Med[,4])

for(ii in 1:1000){
  
  RD <- rep(0,N)
  for(jj in 1:N){
    RD[jj] <- sample( IND[jj,], 1)
  }
  
  
  No.Adj.Tstat[,ii] <- apply( V1[RD,],
                              2,
                              function(x){t.test(x)$statistic})
  
  IPS.Adj.Tstat[,ii] <- apply( V2[RD,],
                               2,
                               function(x){t.test(x)$statistic})
  
  CPS.Adj.Tstat[,ii] <- apply( V3[RD,],
                               2,
                               function(x){t.test(x)$statistic})
  
  
  
  print(ii)
  
}

print(data.frame(cbind(c("Cluster size","Census region=Northeast","Census region=South","Census region=West",
                         "Location=City","Location=Rural","Public kindergarten",
                         "Male","Age","Race=Asian","Race=Black","Race=Hispanic","Race=White",
                         "Motor skill","Intact family","Parental education $\\geq$ College",
                         "Socioeconomic status"),
                       sprintf("& %0.3f",apply(No.Adj.Tstat,1,median)),
                       sprintf("& %0.3f",apply(IPS.Adj.Tstat,1,median)),
                       sprintf("& %0.3f",apply(CPS.Adj.Tstat,1,median)),"\\\\ \\hline")),
      row.names=FALSE)



layout(matrix(c(1,2,3,4,5,3),2,3,byrow=T),widths=c(3,3,1),heights=c(1.5,1))
par(mar=c(0,3,1.5,0.5))

hist(IPS.Med[Data$A==1,4],xlim=c(0.4,1),ylim=c(0,6),breaks=seq(0,1,by=0.01),prob=T,
     col=rgb(1,0,0,0.2),border=NA,xlab="",ylab="Density",xaxt='n',main="")
par(new=T)
hist(1-IPS.Med[Data$A==0,4],xlim=c(0.4,1),ylim=c(0,6),breaks=seq(0,1,by=0.01),prob=T,
     col=rgb(0,0,1,0.2),border=NA,xlab="",ylab="",xaxt='n',yaxt='n',main="")
axis(1,at=c(-10,10))
title(main=expression("Propensity Score (e)"),cex=1.2)

hist(CPS.Med[Data$A==1,4],xlim=c(0.3,1),ylim=c(0,6),breaks=seq(0,1,by=0.01),prob=T,
     col=rgb(1,0,0,0.2),border=NA,xlab="",ylab="",xaxt='n',main="")
par(new=T)
hist(1-CPS.Med[Data$A==0,4],xlim=c(0.3,1),ylim=c(0,6),breaks=seq(0,1,by=0.01),prob=T,
     col=rgb(0,0,1,0.2),border=NA,xlab="",ylab="",xaxt='n',yaxt='n',main="")
axis(1,at=c(-10,10))
title(main=expression("Conditional Propensity Score ("*pi*")"),cex=1.2)

par(mar=c(3,0,1,0))
plot.new()
points(c(0.1,0.1),c(0.7,0.3),pch=c(15,15),cex=c(2,2),col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))
text(c(0.2,0.2),c(0.7,0.3),c("PS Estimate","PS Estimate"),pos=4)
text(c(0.2,0.2),c(0.65,0.25),c("under A=1","under A=0"),pos=4)

par(mar=c(3,3,0,0.5))
boxplot(IPS.Med[,4]*Data$A + (1-IPS.Med[,4])*(1-Data$A)~(1-Data$A),horizontal=TRUE,ylim=c(0.4,1),frame=F,
        col=c(rgb(0,0,1,0.2),rgb(1,0,0,0.2)),yaxt='n',xlim=c(0,3),outpch=20,outcex=0.25)

boxplot(CPS.Med[,4]*Data$A + (1-CPS.Med[,4])*(1-Data$A)~(1-Data$A),horizontal=TRUE,ylim=c(0.3,1),frame=F,
        col=c(rgb(0,0,1,0.2),rgb(1,0,0,0.2)),yaxt='n',xlim=c(0,3),outpch=20,outcex=0.25)

