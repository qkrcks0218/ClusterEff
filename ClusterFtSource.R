############################################
# Last update : October 09 2021
# Collection of Functions
############################################

library(lme4)
library(geepack)

logit <- function(x){
    return( log (x) - log(1-x) )
}

indPS <- function(X.pos, pos, PS.function, v){
    if(length(v)==1){
        indPS <- matrix(0,length(pos),1)
        for(jj in 1:length(pos)){ 
            indPS[jj] <- PS.function(matrix( X.pos[jj,],1,length(X.pos[jj,]) ),v)
        } 
    } else {
        indPS <- matrix(0,length(pos),length(v))
        for(jj in 1:length(pos)){
            for(kk in 1:length(v)){
                indPS[jj,kk] <- PS.function(matrix( X.pos[jj,],1,length(X.pos[jj,]) ),v[kk])
            }
        }
    }
    return(indPS)
}

JointPS.V <- function(A.pos, X.pos, pos, PS.function, v, sV){
    IP <- indPS(X.pos, pos, PS.function, v)
    OIP <- IP*matrix(A.pos,dim(IP)[1],dim(IP)[2]) + (1-IP)*(1-matrix(A.pos,dim(IP)[1],dim(IP)[2]))
    if(sV<10^(-5)){
        apply( OIP , 2, prod )
    } else {
        apply( OIP , 2, prod )*dnorm(v,mean=0,sd=sV)
    }
    
}

JointPS <- function( A,
                     X,
                     PS.function,
                     GP,
                     sV){
    ## Individual level PS
    GP.ID <- unique(GP)
    Joint.PS <- rep(0, length(GP.ID))
    for(ii in 1:length(GP.ID)){
        pos <- which(GP==GP.ID[ii])
        if(!is.null(dim(X))){
            X.pos <- matrix(X[pos,],length(pos),dim(X)[2])
        } else {
            X.pos <- matrix(X[pos],length(pos),1)
        }
        A.pos <- A[pos]
        if(sV<10^(-5)){
            Joint.PS[ii] <- JointPS.V(A.pos,X.pos,pos,PS.function,0,sV)
        } else {
            Joint.PS[ii] <- integrate( function(v){ JointPS.V(A.pos,X.pos,pos,PS.function,v,sV) }, 
                                       lower=-10*(floor(length(A.pos)/10)+1)*sV, upper=10*(floor(length(A.pos)/10)+1)*sV, rel.tol =.Machine$double.eps^0.75)$value
        }
    }
    return(Joint.PS)
}

JointPS.Others <- function( A,
                            X,
                            PS.function,
                            GP,
                            sV){
    ## Individual level PS
    GP.ID <- unique(GP)
    Joint.PS <- rep(0, length(A))
    for(ii in 1:length(GP.ID)){
        pos.original <- which(GP==GP.ID[ii])
        if(length(pos.original)>1){
            for(kkk in 1:length(pos.original)){
                pos <- pos.original[-kkk]
                if(!is.null(dim(X))){
                    X.pos <- matrix(X[pos,],length(pos),dim(X)[2])
                } else {
                    X.pos <- matrix(X[pos],length(pos),1)
                }
                A.pos <- A[pos]
                if(sV<10^(-5)){
                    Joint.PS[pos.original[kkk]] <- JointPS.V(A.pos,X.pos,pos,PS.function,0,sV)
                } else {
                    Joint.PS[pos.original[kkk]] <- integrate( function(v){ JointPS.V(A.pos,X.pos,pos,PS.function,v,sV) }, 
                                                              lower=-10*(floor(length(A.pos)/10)+1)*sV, upper=10*(floor(length(A.pos)/10)+1)*sV, rel.tol =.Machine$double.eps^0.75)$value
                }
            }
        } else {
            Joint.PS[pos.original] <- 1
        }
        
        
    }
    return(Joint.PS)
}

Ind.A <- function(A,weight=NULL,type="IPS"){
    if(is.null(weight)){
        m <- length(A)
        return( matrix(as.numeric(A==1) - as.numeric(A==0),length(A),1)/(m*2^(m-1)) )
    } else {
        if(type=="IPS"){
            III <- matrix(0,length(A),1)
            for(iiiter in 1:length(A)){
                III[iiiter] <- ( 2*A[iiiter]-1 )*prod(weight[-iiiter])/length(A)
            }
            return(III)
        } else if(type=="GPS") {
            III <- matrix(0,length(A),1)
            for(iiiter in 1:length(A)){
                III[iiiter] <- ( 2*A[iiiter]-1 )*weight[iiiter]/length(A)
            }
            return(III)
        } else if(type=="CPS"){
            III <- matrix(0,length(A),1)
            for(iiiter in 1:length(A)){
                III[iiiter] <- ( 2*A[iiiter]-1 )*weight[iiiter]/length(A)
            }
            return(III)
        }
        
    }
    
}

B <- function(M,beta){
    B <- matrix(-beta,M,M)
    diag(B) <- 1
    return(B)
}

expit <- function(x){
    return( (1+exp(-x))^(-1) )
}

obtainB <- function(estb,bound=TRUE,lb=NULL,ub=NULL,name=NULL){
    if(is.null(name)){
        if(bound==TRUE){
            cand <- cbind(as.numeric( gsub("as.factor\\(Mf\\)","",names(estb)) ),apply(cbind(apply(cbind(lb,estb),1,max),ub),1,min))
            return( rbind( cand ) )
        } else {
            cand <- cbind(as.numeric( gsub("as.factor\\(Mf\\)","",names(estb)) ),estb)
            return( rbind( cand ) )
        }
    } else {
        if(bound==TRUE){
            cand <- cbind(name,apply(cbind(apply(cbind(lb,estb),1,max),ub),1,min))
            return( rbind( cand ) )
        } else {
            cand <- cbind(name,estb)
            return( rbind( cand ) )
        }
    }
    
    
}

obtainV <- function(input,split=2){
    if (split==2) {
        return( (var(input[SS$MS])*(length(SS$MS)-1)+
                        var(input[SS$AS])*(length(SS$AS)-1) )/length(input) )
    } else if (split==1) {
        return( var(input)*(length(input)-1)/length(input) )
    }
    
}




IPS.Estimation <- function(A,X,Cid,type="ML",SL.hpara=NULL,CVlist=NULL){
    p <- dim(X)[2]
    tempD <- data.frame(cbind(Cid,A,X))
    colnames(tempD)[1] <- "Cid"
    colnames(tempD)[2] <- "A"
    if(type=="ML"){
        output <- MySL(tempD,locY=2,locX=3:(p+2),Ydist=binomial(),
                       SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN, 
                       CVlist=CVlist)
    } else if (type=="GLMM") {
        output <- glmer(as.formula(paste("A~",paste(colnames(tempD)[3:(p+2)],collapse="+"),"+(1|Cid)",sep="")),data=tempD,
                        family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000) ) )
    } else if (type=="GLM") {
        output <- glm(as.formula(paste("A~",paste(colnames(tempD)[3:(p+2)],collapse="+"),sep="")),data=tempD,family="binomial")
    } else if (type=="GEE"){
        output <- geeglm(as.formula(paste("A~",paste(colnames(tempD)[3:(p+2)],collapse="+"),sep="")),id=Cid,data=tempD,family="binomial",corstr="exchangegable")
    }
    result <- list()
    result$type <- paste("PS_",type,sep="")
    result$output <- output
    return(result)
}

OR.Estimation <- function(Y,A,X,Cid,type="ML",outcome.type=gaussian(),SL.hpara=NULL,CVlist=NULL){
    p <- dim(X)[2]
    tempD <- data.frame(cbind(Cid,Y,A,X))
    colnames(tempD)[1] <- "Cid"
    colnames(tempD)[2] <- "Y"
    colnames(tempD)[3] <- "A"
    
    if(type=="ML"){
        output <- MySL(tempD,locY=2,locX=3:(p+3),Ydist=outcome.type,
                       SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN, CVlist=CVlist)
    } else if (type=="GLMM") {
        if(outcome.type$family=="binomial"){
            output <- glmer(as.formula(paste("Y~A*(",paste(colnames(tempD)[4:(p+3)],collapse="+"),")+(1|Cid)",sep="")),data=tempD,
                            family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000) ) )
        } else if (outcome.type$family=="gaussian"){
            output <- lmer(as.formula(paste("Y~A*(",paste(colnames(tempD)[4:(p+3)],collapse="+"),")+(1|Cid)",sep="")),data=tempD)
        }
    } else if (type=="GLM") {
        if(outcome.type$family=="binomial"){
            output <- glm(as.formula(paste("Y~A*(",paste(colnames(tempD)[4:(p+3)],collapse="+"),")",sep="")),data=tempD,family="binomial")
        } else if (outcome.type$family=="gaussian"){
            output <- lm(as.formula(paste("Y~A*(",paste(colnames(tempD)[4:(p+3)],collapse="+"),")",sep="")),data=tempD)
        }
    } else if (type=="GEE") {
        if(outcome.type$family=="binomial"){
            output <- geeglm(as.formula(paste("Y~A*(",paste(colnames(tempD)[4:(p+3)],collapse="+"),")",sep="")),id=Cid,data=tempD,family="binomial",corstr="exchangegable")
        } else if (outcome.type$family=="gaussian"){
            output <- geeglm(as.formula(paste("Y~A*(",paste(colnames(tempD)[4:(p+3)],collapse="+"),")",sep="")),id=Cid,data=tempD,corstr="exchangegable")
        }
    }
    
    result <- list()
    result$type <- paste("OR_",type,sep="")
    result$output <- output
    return(result)
    
}

IPS.Prediction <- function(IPS.Est,A,X,Cid){
    p <- dim(X)[2]
    tempD <- data.frame(cbind(Cid,A,X))
    colnames(tempD)[1] <- "Cid"
    colnames(tempD)[2] <- "A"
    if(IPS.Est$type=="PS_ML"){
        
        IPS.output.1 <- predict(IPS.Est$output,tempD[,3:(p+2)])$pred
        IPS.output <- IPS.output.1*A + (1-IPS.output.1)*(1-A)
        
    } else if(IPS.Est$type=="PS_GLMM"){
        
        gamma <- fixef( IPS.Est$output )
        sV <- as.data.frame(VarCorr(IPS.Est$output))[1,5]
        
        PS.function.IPS <- function(X,V){
            if(!is.null(dim(X)[1])){
                return( expit( cbind(1,X)%*%gamma + V )  )
            } else {
                return( expit( sum(c(1,X)*gamma)+V ) )
            }
        }
        
        IPS.output <- JointPS(A,as.matrix(X),PS.function.IPS,1:length(A),sV)
        IPS.output.1 <- JointPS(rep(1,length(A)),as.matrix(X),PS.function.IPS,1:length(A),sV)
        
        
    } else if(IPS.Est$type=="PS_GLM"|IPS.Est$type=="PS_GEE"){
        
        gamma <- IPS.Est$output$coefficients
        
        PS.function.IPS <- function(X,V){
            if(!is.null(dim(X)[1])){
                return( expit( cbind(1,X)%*%gamma  )  )
            } else {
                return( expit( sum(c(1,X)*gamma) ) )
            }
        }
        
        IPS.output <- JointPS(A,as.matrix(X),PS.function.IPS,1:length(A),0)
        IPS.output.1 <- JointPS(rep(1,length(A)),as.matrix(X),PS.function.IPS,1:length(A),0)
    } 
    
    result <- list()
    result$type <- IPS.Est$type
    result$IPS <- IPS.output
    result$IPS.1 <- IPS.output.1
    return(result)
}


OR.Prediction <- function(OR.Est,A,X,outcome.type=gaussian()){
    p <- dim(X)[2]
    tempD <- data.frame(cbind(A,X))
    colnames(tempD)[1] <- "A"
    
    if(outcome.type$family=="gaussian"){
        if(OR.Est$type=="OR_ML"){
            output <- predict(OR.Est$output,tempD)$pred
        } else if (OR.Est$type=="OR_GLMM"){
            beta <- fixef(OR.Est$output)
            output <- as.matrix(cbind(1,A,X,A*X))%*%beta
        } else if (OR.Est$type=="OR_GLM"|OR.Est$type=="OR_GEE"){
            beta <- OR.Est$output$coefficients
            output <- as.matrix(cbind(1,A,X,A*X))%*%beta
        }
    } else {
        if(OR.Est$type=="OR_ML"){
            output <- predict(OR.Est$output,tempD)$pred
        } else if (OR.Est$type=="OR_GLMM"){
            beta <- fixef(OR.Est$output)
            sV <- as.data.frame(VarCorr(OR.Est$output))[1,5]
            Predictor <- as.matrix(cbind(1,A,X,A*X))%*%beta
            if(sV>=10^(-5)){
                output <- rep(0,dim(X)[1])
                for(ii in 1:length(output)){
                    output[ii] <- integrate( function(v){ (exp(Predictor[ii]+v)/(1+exp(Predictor[ii]+v)))*dnorm(v,mean=0,sd=sV) }, 
                                             lower=max(-10*sV,-20), upper=min(10*sV,20), rel.tol =.Machine$double.eps^0.75)$value
                }    
            } else {
                output <- exp(Predictor)/(1+exp(Predictor))
            }
        } else if (OR.Est$type=="OR_GLM"|OR.Est$type=="OR_GEE"){
            beta <- OR.Est$output$coefficients
            output <- exp(as.matrix(cbind(1,A,X,A*X))%*%beta)/(1+exp(as.matrix(cbind(1,A,X,A*X))%*%beta))
        }
    }
    
    
    result <- list()
    result$type <- OR.Est$type
    result$output <- output
    return(result)
    
}

findbetaDisc <- function(Y,A,Cid,Yhat,IPShat,Clustering){
    
    N <- length(unique(Cid))
    M <- rep(0,N)
    Clvec <- rep(0,N)
    index <- list()
    for(ii in 1:N){
        index[[ii]] <- which( Cid==unique(Cid)[ii] )
        M[ii] <- length(index[[ii]])
        ck <- unique(Clustering[index[[ii]]])
        if(length(ck)>1){
            print("Clustering values of units in the same cluster should be identical.")
            return()
        } else {
            Clvec[ii] <- ck
        }
    }
    
    Residual <- Y - Yhat
    Sigma <- list()
    SY <- rep(0,N)
    SW <- rep(0,N)
    for(ii in 1:N){
        Sigma[[ii]] <- t(t(Residual[ index[[ii]] ]))%*%t(Residual[ index[[ii]] ])
        S0 <- Sigma[[ ii ]]
        S1 <- matrix( apply(S0,2,sum), dim(S0)[1], dim(S0)[1], byrow=T)
        S2 <- matrix( apply(S0,2,sum), dim(S0)[1], dim(S0)[1])
        S3 <- matrix( sum(S0), dim(S0)[1], dim(S0)[1])
        contrast.vec <- Ind.A(A[ index[[ii]] ], weight=1/IPShat[ index[[ii]] ],type="CPS")
        
        SY[ii] <- (t(contrast.vec)%*%( 2*S0 - S1 - S2 )%*%contrast.vec)/ (-2*t(contrast.vec)%*%( S0-S1-S2+S3 )%*%contrast.vec)
        SW[ii] <- t(contrast.vec)%*%( S0-S1-S2+S3 )%*%contrast.vec
    }
    FD <- data.frame(cbind(SY,Clvec,SW))
    colnames(FD) <- c("SY","Mf","SW")
    Mf.unique <- sort(unique(FD$Mf))
    
    valid.Mat <- cbind(unique(Clvec),0)
    valid.Mat <- valid.Mat[valid.Mat[,1]>=1,]
    
    if(is.null(dim(valid.Mat))){
        valid.Mat <- matrix(valid.Mat,1,2)
        temp.b <- cbind(sort(unique(valid.Mat[,1])),lm(SY~1,weights=SW,data=FD[FD$Mf>0,])$coefficients)
    } else {
        temp.b <- cbind(sort(unique(valid.Mat[,1])),lm(SY~0+as.factor(Mf),weights=SW,data=FD[FD$Mf>0,])$coefficients)
    }
    for(jj in 1:dim(valid.Mat)[1]){
        valid.Mat[jj,2] <- temp.b[which(valid.Mat[jj,1]==temp.b[,1]),2]
    }
    if(dim(valid.Mat)[1]>1){
        valid.Mat <- valid.Mat[order(valid.Mat[,1]),]
    }
    
    
    beta.bdd <- NULL
    beta.nobdd <- obtainB( estb=valid.Mat[,2],bound=FALSE,name=valid.Mat[,1])
    
    output <- list()
    output$beta.bdd <- beta.bdd
    output$beta.nobdd <- beta.nobdd
    return(output)
}

findIFwithoutB <- function(Y,A,Cid,Yhat,Y1hat,Y0hat,IPShat){
    
    N <- length(unique(Cid))
    M <- rep(0,N)
    index <- list()
    for(ii in 1:N){
        index[[ii]] <- which( Cid==unique(Cid)[ii] )
        M[ii] <- length(index[[ii]])
    }
    
    Residual <- Y - Yhat
    
    IF.result <- aggregate( (2*A-1)/IPShat*Residual + Y1hat - Y0hat~Cid, FUN=mean)[,2]
    
    return(IF.result)
}

findIFwithB.Disc <- function(Y,A,Cid,Yhat,Y1hat,Y0hat,IPShat,Beta,Clustering){
    
    N <- length(unique(Cid))
    M <- rep(0,N)
    Clvec <- rep(0,N)
    index <- list()
    for(ii in 1:N){
        index[[ii]] <- which( Cid==unique(Cid)[ii] )
        M[ii] <- length(index[[ii]])
        ck <- unique(Clustering[index[[ii]]])
        if(length(ck)>1){
            print("Clustering values of units in the same cluster should be identical.")
            return()
        } else {
            Clvec[ii] <- ck
        }
    }
    
    Residual <- Y - Yhat
    IF.result <- rep(0,N)
    
    for(iter in 1:N){
        if(M[iter]>1){
            beta <- Beta[Beta[,1]==Clvec[iter],2]
            if(length(beta)==0){ beta <- 0 }
            BLP.mat <- B(M[iter],beta)
            contrast.vec <- Ind.A(A[ index[[iter]] ], weight=1/IPShat[ index[[iter]] ],type="CPS")
            IF.result[iter] <- t(contrast.vec)%*%BLP.mat%*%Residual[ index[[iter]] ] + mean( Y1hat[ index[[iter]] ] - Y0hat[ index[[iter]] ] )
            
        } else {
            contrast.vec <- Ind.A(A[ index[[iter]] ], weight=1/IPShat[ index[[iter]] ],type="CPS")
            IF.result[iter] <- t(contrast.vec)%*%Residual[ index[[iter]] ] + mean( Y1hat[ index[[iter]] ] - Y0hat[ index[[iter]] ] )
            
        }
        
    }
    
    return(IF.result)
    
}

ATEcalculate <- function (Data, learner, weight.vector=NULL, clusters=NULL) {
    
    if(length(weight.vector)==0){ weight.vector <- rep(1,dim(Data)[1]) }
    cluster.se <- length(clusters) > 0
    subset.W.orig <- Data$A
    subset.W.hat <- learner$p_hat
    subset.Y.orig <- Data$Y
    subset.Y.hat <- learner$m_hat
    tau.hat.pointwise <- predict(learner)
    
    control.idx <- which(subset.W.orig == 0)
    treated.idx <- which(subset.W.orig == 1)
    tau.avg.raw <- weighted.mean(tau.hat.pointwise, weight.vector)
    
    
    Y.hat.0 <- subset.Y.hat - subset.W.hat * tau.hat.pointwise
    Y.hat.1 <- subset.Y.hat + (1 - subset.W.hat) * tau.hat.pointwise
    
    gamma.control.raw <- 1/(1 - subset.W.hat[control.idx])
    gamma.treated.raw <- 1/subset.W.hat[treated.idx]
    
    gamma <- rep(0, length(subset.W.orig))
    gamma[control.idx] <- gamma.control.raw/sum(weight.vector[control.idx] * gamma.control.raw) * sum(weight.vector)
    gamma[treated.idx] <- gamma.treated.raw/sum(weight.vector[treated.idx] * gamma.treated.raw) * sum(weight.vector)
    dr.correction.all <- subset.W.orig * gamma * (subset.Y.orig - Y.hat.1) - (1 - subset.W.orig) * gamma * (subset.Y.orig - Y.hat.0)
    dr.correction <- weighted.mean(dr.correction.all, weight.vector)
    if (cluster.se) {
        correction.clust <- Matrix::sparse.model.matrix(~factor(clusters) + 0, transpose = TRUE) %*% (dr.correction.all * weight.vector)
        sigma2.hat <- sum(correction.clust^2)/sum(weight.vector)^2 * length(correction.clust)/(length(correction.clust) - 1)
    } else {
        sigma2.hat <- mean(dr.correction.all^2)/(length(dr.correction.all) - 1)
    }
    
    tau.avg <- tau.avg.raw + dr.correction
    tau.se <- sqrt(sigma2.hat)
    return(c(estimate = tau.avg, std.err = tau.se))
}

