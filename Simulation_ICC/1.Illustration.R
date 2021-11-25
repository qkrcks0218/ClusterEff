############################################
# Last update : October 09 2021
# Illustration of Comparison between \bar{tau} and \hat{\tau}
############################################

PATH <- getwd()

source("../MySL.R")
source("../ClusterFtSource.R")

ATE <- matrix(0,4410,9)
ATE <- data.frame(ATE)
colnames(ATE) <- c("TITER","sV","sU",
                   "M1_CY_Plug","M1_CY_Ours",
                   "M2_CY_Plug","M2_CY_Ours",
                   "OR_ICC","PS_ICC")

for(BATCH in 1:4410){
    
    Vgrid <- seq(0,2,length=21)
    Ugrid <- seq(0,2,length=21)
    Pmat <- expand.grid(Vgrid,Ugrid)
    sV <- Pmat[(BATCH-1)%%441+1,1]
    sU <- Pmat[(BATCH-1)%%441+1,2]
    
    TITER <- 100
    
    ATE1 <- matrix(0,TITER,2)
    ICC <- matrix(0,TITER,4)
    
    for(ii in 1:TITER){
        
        N <- 200
        M <- round( runif(N, min=1000/N*0.8, max=1000/N*1.2) )
        M.start <- cumsum(M)-M+1
        M.end <- cumsum(M)
        GP <- rep(1:N, M)
        
        
        temp <- sample(N,N)
        
        SSMS.C <- sort( temp[1:round(N/2)] )
        SSAS.C <- (1:N)[-SSMS.C]
        
        SSM <- SSA <- NULL
        for(jj in 1:N){
            pos <- which(GP==jj)
            if(sum(SSMS.C==jj)>0){ SSM <- c(SSM,pos) } else if (sum(SSAS.C==jj)) { SSA <- c(SSA,pos) } 
        }
        
        SS <- list()
        SS$MS <- SSMS.C
        SS$AS <- SSAS.C
        
        SSI <- list()
        SSI$MS <- SSM
        SSI$AS <- SSA
        
        
        V.GP <- rnorm(N)*sV
        V <- rep(V.GP,M)
        
        Cl.X <- cbind( rep(rnorm(N),M), rep(rbinom(N,1,0.7),M) )
        X <- cbind( rnorm(sum(M)), rnorm(sum(M)), rbinom(sum(M),1,0.3), Cl.X)
        
        
        
        PS.function <- function(X,V){
            return( expit(-0.5 + 0.5*X[,1] - as.numeric(X[,2]>1) + 0.5*X[,3] - 0.25*X[,4] + X[,5] + V )  )
        }
        
        PS <- PS.function(X,V)
        A <- rbinom(sum(M),1,PS)
        True.GPS <- JointPS(A,X,PS.function,GP,sV) 
        True.IPS <- JointPS(A,X,PS.function,1:sum(M),sV)
        True.OPS <- JointPS.Others(A,X,PS.function,GP,sV)
        True.CPS <- rep(True.GPS,M)/True.OPS
        
        U <- rep(rnorm(N)*sU, M)
        eps <- rnorm(sum(M))
        
        Y <- true.g <- true.g.1 <- true.g.0 <- list()
        
        Y[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*A + 2*X[,1] - X[,4]^2 + X[,2]*X[,5] + U + eps
        true.g[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*A + 2*X[,1] - X[,4]^2 + X[,2]*X[,5]
        true.g.1[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*1 + 2*X[,1] - X[,4]^2 + X[,2]*X[,5]
        true.g.0[[1]] <- 3 + (2.1 + X[,2]^2 + 3*X[,3])*0 + 2*X[,1] - X[,4]^2 + X[,2]*X[,5]
        
        ATE1[ii,1] <- mean( findIFwithoutB(Y[[1]],A,GP,true.g[[1]],true.g.1[[1]],true.g.0[[1]],True.IPS) )
        
        Clustering <- rep(M,M)
        VV <- rep(0,N)
        
        Beta.MS <- findbetaDisc(Y[[1]][SSI$MS],A[SSI$MS],GP[SSI$MS],true.g[[1]][SSI$MS],True.CPS[SSI$MS],Clustering=Clustering[SSI$MS])$beta.nobdd
        Beta.AS <- findbetaDisc(Y[[1]][SSI$AS],A[SSI$AS],GP[SSI$AS],true.g[[1]][SSI$AS],True.CPS[SSI$AS],Clustering=Clustering[SSI$AS])$beta.nobdd
        
        VV[SS$MS] <- findIFwithB.Disc(Y[[1]][SSI$MS],A[SSI$MS],GP[SSI$MS],true.g[[1]][SSI$MS],true.g.1[[1]][SSI$MS],true.g.0[[1]][SSI$MS],True.CPS[SSI$MS],Beta.MS,Clustering=Clustering[SSI$MS])
        VV[SS$AS] <- findIFwithB.Disc(Y[[1]][SSI$AS],A[SSI$AS],GP[SSI$AS],true.g[[1]][SSI$AS],true.g.1[[1]][SSI$AS],true.g.0[[1]][SSI$AS],True.CPS[SSI$AS],Beta.AS,Clustering=Clustering[SSI$AS])
        
        ATE1[ii,2] <- mean( VV )
        
        
        Gp.Sum <- aggregate(Y[[1]]-true.g[[1]]~GP,FUN="sum")[,2]
        MSB <- (sum(Gp.Sum^2/M) - sum(Gp.Sum)^2/sum(M))/(N-1)
        MSW <- (sum((Y[[1]]-true.g[[1]])^2)-sum(Gp.Sum^2/M))/(sum(M)-N)
        N.m <- (sum(M) - sum(M^2)/sum(M))/(N-1)
        ICC1 <- (MSB-MSW)/(MSB+(N.m-1)*MSW)
        
        Gp.Sum <- aggregate(A-True.IPS~GP,FUN="sum")[,2]
        MSB <- (sum(Gp.Sum^2/M) - sum(Gp.Sum)^2/sum(M))/(N-1)
        MSW <- (sum((A-True.IPS)^2)-sum(Gp.Sum^2/M))/(sum(M)-N)
        N.m <- (sum(M) - sum(M^2)/sum(M))/(N-1)
        ICC3 <- (MSB-MSW)/(MSB+(N.m-1)*MSW)
        
        ICC[ii,] <- c(sV,sU,ICC1,ICC3)
        
        
        print(ii)
    }
    
    ATE[BATCH,] <- c(TITER,sV,sU,apply(ATE1,2,sum),apply(ATE1^2,2,sum),apply(ICC[,3:4],2,mean))
    
}

write.csv(ATE,"Illustration_Grid.csv",row.names=FALSE)







R.raw <- read.csv("Illustration_Grid.csv")
RR <- aggregate(.~sV+sU,data=R.raw,FUN="sum")

Summary <- cbind(RR$sV,RR$sU,
                 RR$M2_CY_Plug/RR$TITER -(RR$M1_CY_Plug/RR$TITER)^2,
                 RR$M2_CY_Ours/RR$TITER -(RR$M1_CY_Ours/RR$TITER)^2)
Summary <- data.frame(Summary)
colnames(Summary) <- c("sV","sU",
                       "CY_P","CY_O")


RR2 <- aggregate(.~sV,data=R.raw,FUN="mean")
ICC <- round( RR2$PS_ICC[c(1,6,11,16,21)], 2)


OCY <- matrix(Summary$CY_P/Summary$CY_O,21,21)
ODY <- matrix(Summary$DY_P/Summary$DY_O,21,21)

CY <- DY <- matrix(0,21,21)
ww <- 0
for(ii in 1:21){
    for(jj in 1:21){
        II <- (max(ii-ww,0)):(min(ii+ww,21))
        JJ <- (max(jj-ww,0)):(min(jj+ww,21))
        CY[ii,jj] <- mean(OCY[II,JJ])
        DY[ii,jj] <- mean(ODY[II,JJ])
    }
}




filled.contour3 <-
    function (x = seq(0, 1, length.out = nrow(z)),
              y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes,mar, ...) 
    {
        # modification by Ian Taylor of the filled.contour function
        # to remove the key and facilitate overplotting with contour()
        # further modified by Carey McGilliard and Bridget Ferris
        # to allow multiple plots on one page
        
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                }
                else {
                    z <- x
                    x <- seq.int(0, 1, length.out = nrow(z))
                }
            }
            else stop("no 'z' matrix specified")
        }
        else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
        if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
            stop("increasing 'x' and 'y' values expected")
        # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        # on.exit(par(par.orig))
        # w <- (3 + mar.orig[2]) * par("csi") * 2.54
        # par(las = las)
        # mar <- mar.orig
        plot.new()
        # par(mar=mar)
        plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
        if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
            stop("no proper 'z' matrix specified")
        if (!is.double(z)) 
            storage.mode(z) <- "double"
        .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                        col = col)
        if (missing(plot.axes)) {
            if (axes) {
                title(main = "", xlab = "", ylab = "")
                Axis(x, side = 1)
                Axis(y, side = 2)
            }
        }
        else plot.axes
        if (frame.plot) 
            box()
        if (missing(plot.title)) 
            title(...)
        else plot.title
        invisible()
    }


filled.legend <-
    function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                           length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes, ...) 
    {
        # modification of filled.contour by Carey McGilliard and Bridget Ferris
        # designed to just plot the legend
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                }
                else {
                    z <- x
                    x <- seq.int(0, 1, length.out = nrow(z))
                }
            }
            else stop("no 'z' matrix specified")
        }
        else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
        if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
            stop("increasing 'x' and 'y' values expected")
        #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        #  on.exit(par(par.orig))
        #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
        #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
        #  par(las = las)
        #  mar <- mar.orig
        #  mar[4L] <- mar[2L]
        #  mar[2L] <- 1
        #  par(mar = mar)
        # plot.new()
        plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                    yaxs = "i")
        rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
        if (missing(key.axes)) {
            if (axes) 
                axis(4)
        }
        else key.axes
        box()
    }



C11 <- function(n,alpha=1){
    l1 <- 11
    l2 <- n-11+1
    cols <- c(if (l1 > 0) {hsv(h = 180/360, s = seq.int(0.5,if (TRUE){ 0.5/(l1-1) } else { 0 }, length.out = l1), v = 1, alpha = alpha)} , 
              if (l2 > 1) {hsv(h = 0, s = seq.int(0, 0.75, length.out = l2)[-1L], v = 1, alpha = alpha)} )
    cols
    
}
C2 <- function(n,alpha=1){
    cols <- hsv(h = 240/360, s = seq.int(0.05,0.7, length.out = n), v = 1, alpha = alpha)
    cols
    
}




g1 <- g2 <- seq(0,2,length=21)

layout(matrix(c(1,1,1,1,1,2,
                1,1,1,1,1,2),2,6,byrow=T))

par(mar=c(3.5,3.5,3,0.5),oma=c(2,2,0,0),xpd=TRUE)

range(CY)

LY <- seq(0.8,4.05,by=0.05)

C1 <- function(n,alpha=1){
    l1 <- which(LY==1)-1
    l2 <-  n-l1+1
    cols <- c(if (l1 > 0) {hsv(h = 180/360, s = seq.int(0.5,if (TRUE){ 0.5/(l1-1) } else { 0 }, length.out = l1), v = 1, alpha = alpha)} , 
              if (l2 > 1) {hsv(h = 0, s = seq.int(0, 0.75, length.out = l2)[-1L], v = 1, alpha = alpha)} )
    cols
    
}

filled.contour3( g1, g2, 
                 CY, color = C1, zlim=range(LY), levels=LY,
                 plot.axes = { 
                     points(c(0.0,0.0,1.5,1.5),c(0.5,1.5,0.5,1.5),pch=c(1,2,3,4),cex=2,lwd=2,col=4);
                     contour(g1, g2, CY, levels=c(0.9,1,2,3),add=TRUE,lwd=c(1,2,1,1),col=rgb(0,0,0,0.5))})
title(xlab="",ylab="",line=2.5)
title(main=expression(rho(sigma["V"],sigma["U"])),cex.main=1.5)
axis(2,at=c(0,0.5,1,1.5,2))
axis(2,at=c(0,0.5,1,1.5,2),
     label=sprintf("(%0.2f)",c(0,0.25/1.25,0.5,1.5^2/(1+1.5^2),0.8)), col.axis="green3", col=NA, col.ticks=NA,line=1.25)
axis(1,at=c(0,0.5,1,1.5,2))
axis(1,at=c(0,0.5,1,1.5,2),
     label=sprintf("(%0.2f)",ICC), col.axis="green3", col=NA, col.ticks=NA,line=1.25)




plot.new()
par(mar=c(3.5,3.5,3,3))

filled.legend( c(0,1), c(0,1), matrix(c(-3.5,2.5,0,0),2,2),
               color = C1, zlim=range(LY), levels=c(0.8,0.85,0.9,0.95,seq(1,4,by=0.2),4.05) )


mtext(expression("Treatment Correlation ("*sigma["V"]*")                   "),side=1,line=0.5,outer=TRUE,cex=0.9)
mtext("                                   (ICC)",side=1,line=0.25,outer=TRUE,cex=0.9,col="green3")
mtext(expression("Outcome Correlation ("*sigma["U"]*")                   "),side=2,line=0.5,outer=TRUE,cex=0.9)
mtext("                                   (ICC)",side=2,line=0.75,outer=TRUE,cex=0.9,col="green3")







