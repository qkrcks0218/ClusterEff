############################################
# Last update : October 09 2021
# Collection of Super Learner Functions
############################################

library(SuperLearner)
library(caret)
library(nnet)
library(glmnet)
library(earth)
library(gam)
library(gbm)
library(xgboost)     
library(kernlab)
library(polspline)
library(ranger)

MySL <- function( Data, locY, locX, Ydist=gaussian(), 
                  SL.list=c(1:11), MTRY=c(2,4,6,8), MLPL=c(2,4,6,8), MLPdecay=c(10^(-4),10^(-5)), NMN=c(20), obsWeights=NULL,
                  PS.thr = 10^(-3), CVlist=NULL ){
    
    ## Poisson c(2,4,5)
    if(Ydist$family=="poisson"){
        SL.list <- intersect(c(1,2,4,5),SL.list)
    }
    
    Learners <- list()
    
    ##############
    # Caret Based
    ##############
    
    SL.caret.SLP <- function (Y, X, newX, family, obsWeights, method = "mlpML", 
                              L1,L2,L3,decay,
                              trControl = caret::trainControl(method = "none"), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) 
    {
        if (family$family == "gaussian") {
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      preProc =  c('center', 'scale', 'pca'),
                                      hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                      learnFuncParams=decay,
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      hiddenActFunc = "Act_Identity",
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    SL.caret.MLP <- function (Y, X, newX, family, obsWeights, method = "mlpML", 
                              L1,L2,L3,decay,
                              trControl = caret::trainControl(method = "none"), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) 
    {
        if (family$family == "gaussian") {
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      learnFuncParams=decay,
                                      hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      hiddenActFunc = "Act_Identity",
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    
    SL.caret.gbm <- function (Y, X, newX, family, obsWeights, method = "gbm",
                              ntree,intdepth,sh,nmn,
                              trControl = caret::trainControl(method = "none"),
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...)
    {
        if (family$family == "gaussian") {
            
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights,
                                      metric = metric, method = method,
                                      tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                      metric = metric, method = method,
                                      tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    #############################################
    # SL-based
    #############################################
    
    SL.new.earth <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                              nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...)
    {
        if (family$family == "gaussian") {
            fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights, 
                                      nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                      ncross = ncross, minspan = minspan, endspan = endspan)
        }
        if (family$family == "binomial") {
            fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights,
                                      nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                      ncross = ncross, minspan = minspan, endspan = endspan, 
                                      glm = list(family = binomial))
        }
        pred <- predict(fit.earth, newdata = newX, type = "response")
        fit <- list(object = fit.earth)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.earth")
        return(out)
    }
    
    SL.new.xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                nthread = 1, verbose = 0, save_period = NULL, ...) 
    {
        if (packageVersion("xgboost") < 0.6) 
            stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
        if (!is.matrix(X)) {
            X = model.matrix(~. - 1, X)
        }
        xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
        if (family$family == "gaussian") {
            model = xgboost::xgboost(data = xgmat, objective = "reg:linear", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, nthread = nthread, 
                                     params = params, save_period = save_period)
        }
        if (family$family == "binomial") {
            model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, nthread = nthread, 
                                     params = params, save_period = save_period)
        }
        if (family$family == "multinomial") {
            model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                                     nthread = nthread, params = params, save_period = save_period)
        }
        if (family$family == "poisson") {
            model = xgboost::xgboost(data = xgmat, objective = "count:poisson", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, 
                                     nthread = nthread, params = params, save_period = save_period)
        }
        if (!is.matrix(newX)) {
            newX = model.matrix(~. - 1, newX)
        }
        pred = predict(model, newdata = newX)
        fit = list(object = model)
        class(fit) = c("SL.xgboost")
        out = list(pred = pred, fit = fit)
        return(out)
    }
    
    Learners[[1]] <- create.Learner("SL.glm")
    TOTAL.M <- 1
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.glmnet",tune=list(alpha=c(1,0.5,0),useMin=c(TRUE,FALSE))) 
    # Lasso.min , EN.min , Ridge.min , Lasso.1se , EN.1se , Ridge.1se
    TOTAL.M <- TOTAL.M+1           # 2
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.new.earth",tune=list(degree=c(1,2,3,4,5)))
    # Earth.deg=1 , Earth.deg=2 , Earth.deg=3 , Earth.deg=4 , Earth.deg=5
    TOTAL.M <- TOTAL.M+1           # 3
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.gam",tune=list(deg.gam=c(1,2,3,4,5)))
    # Gam.deg=1 , Gam.deg=2 , Gam.deg=3 , Gam.deg=4 , Gam.deg=5
    TOTAL.M <- TOTAL.M+1           # 4
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.new.xgboost",tune=list(n.trees=c(100,300,500),max_depth=c(1,2,3,4)))
    # xgboost
    TOTAL.M <- TOTAL.M+1           # 5
    
    # Learners[[5]] <- create.Learner("SL.kernelKnn",tune=list(k=c(1,5,10,20)))
    # knn
    
    # Learners[[5]] <- create.Learner("SL.ksvm",tune=list(kernel=c("rbfdot","polydot","tanhdot")))
    # SVM
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.polymars",tune=list(knots=c(2,3,4)))
    # polspline (similar to earth)
    TOTAL.M <- TOTAL.M+1           # 6
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.ranger",tune=list(num.trees=c(500,1000,1500),mtry=MTRY))
    # RF
    TOTAL.M <- TOTAL.M+1           # 7
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.nnet",tune=list(linout=c(TRUE,FALSE), decay=c(0,0.1),size=MTRY))
    # nnet
    TOTAL.M <- TOTAL.M+1           # 8
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.gbm",tune=list(ntree=c(100,300,500),intdepth=c(1,2,3),sh=c(0.1,0.01),nmn=NMN))
    # gbm
    TOTAL.M <- TOTAL.M+1           # 9
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.SLP", tune = list(L1=MLPL,L2=c(0),L3=c(0),decay=MLPdecay))
    # 1layer
    TOTAL.M <- TOTAL.M+1           # 10
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.MLP", tune = list(L1=MLPL,L2=MLPL,L3=MLPL,decay=MLPdecay))
    # 3layer
    TOTAL.M <- TOTAL.M+1           # 11
    
    # Learners[[12]] <- create.Learner("SL.bartMachine",tune=list(num_trees=c(50,100,150),mem_cache_for_speed=FALSE))
    # bart
    
    BaseLearner <- Learners[[ SL.list[1] ]]$names
    
    if( length(SL.list)>1 ){
        for( METHOD in 2:length(SL.list) ){
            BaseLearner <- c(BaseLearner, Learners[[ SL.list[METHOD] ]]$names)
        }
    }
    
    if( length(locX)==1 ){
        dX <- data.frame(matrix( Data[,locX], dim(Data)[1], 1))
        colnames(dX) <- colnames(Data)[locX]
    } else {
        dX <- Data[,locX]
    }
    
    if(is.null(CVlist)){
        CVL <- list(V = 5)
    } else {
        CVL <- list(V = 5, shuffle = FALSE, validRows = CVlist)
    }
    
    
    if(Ydist$family!="binomial"){
        capture.output( Fitted.SL <- SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                                  SL.library=BaseLearner, cvControl = CVL, obsWeights=obsWeights) , file=NULL )
        
    } else {
        capture.output( Fitted.SL <- SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                                   SL.library=BaseLearner, cvControl = CVL, obsWeights=obsWeights) , file=NULL )
        # POS.Z <- which(apply(Fitted.SL2$Z,2,min)>=PS.thr & apply(Fitted.SL2$Z,2,max)<=1-PS.thr)
        # capture.output( Fitted.SL <- PS.Adjust(Fitted.SL2,POS.Z,
        #                                        Y=Data[,locY],X=dX,family=Ydist,
        #                                        SL.library=BaseLearner,cvControl = list(V = 5), obsWeights=obsWeights) , file=NULL )
        
        
    }
    
    
    return(Fitted.SL)
}




PS.Adjust <- function (Fitted.SL2,POS.Z,
                       Y, X, newX = NULL, family = gaussian(), SL.library, 
                       method = "method.NNLS", id = NULL, verbose = FALSE, 
                       control = list(), cvControl = list(), obsWeights = NULL, 
                       env = parent.frame()) {
    
    SL.library <- SL.library[POS.Z]
    
    .createLibrary <- function(SL.library) {
        if (is.character(SL.library)) {
            k <- length(SL.library)
            whichScreen <- matrix(1, nrow = 1, ncol = k)
            screenAlgorithm <- "All"
            library <- data.frame(predAlgorithm = SL.library, rowScreen = 1, stringsAsFactors=FALSE)
        } else if (is.list(SL.library)) {
            predNames <- sapply(SL.library, FUN = "[", 1)
            NumberScreen <- (sapply(SL.library, FUN = length) - 1)
            if (sum(NumberScreen == 0) > 0) {
                for(ii in which(NumberScreen == 0)) {
                    SL.library[[ii]] <- c(SL.library[[ii]], "All")
                    NumberScreen[ii] <- 1
                }
            }
            screenAlgorithmFull <- unlist(lapply(SL.library, FUN="[", -1))
            screenAlgorithm <- unique(screenAlgorithmFull)
            
            library <- data.frame(predAlgorithm = rep(predNames, times=NumberScreen), rowScreen = match(screenAlgorithmFull, screenAlgorithm), stringsAsFactors = FALSE)
        } else {
            stop('format for SL.library is not recognized')
        }
        
        out <- list(library = library, screenAlgorithm = screenAlgorithm)
        return(out)
    }
    library <- .createLibrary(SL.library)
    
    
    
    if (is.character(method)) {
        if (exists(method, mode = "list")) {
            method <- get(method, mode = "list")
        }
        else if (exists(method, mode = "function")) {
            method <- get(method, mode = "function")()
        }
    }
    if (!is.null(method$require)) {
        sapply(method$require, function(x) require(force(x), 
                                                   character.only = TRUE))
    }
    control <- do.call("SuperLearner.control", control)
    cvControl <- do.call("SuperLearner.CV.control", cvControl)
    varNames <- colnames(X)
    N <- dim(X)[1L]
    p <- dim(X)[2L]
    k <- length(POS.Z)
    Z <- Fitted.SL2$Z[,POS.Z]
    if(is.null(dim(Z))){
        Z <- matrix(Z,N,1)
    }
    libraryNames <- paste(library$library$predAlgorithm, library$screenAlgorithm[library$library$rowScreen], 
                          sep = "_")
    
    fitLibEnv <- new.env()
    assign("fitLibrary", vector("list", length = k), 
           envir = fitLibEnv)
    
    if (is.null(obsWeights)) {
        obsWeights <- rep(1, N)
    }
    if (is.null(newX)) {
        newX <- X
    }
    validRows <- CVFolds(N = N, id = id, Y = Y, cvControl = cvControl)
    if (is.null(id)) {
        id <- seq(N)
    }
    errorsInCVLibrary <- Fitted.SL2$errorsInCVLibrary[POS.Z]
    
    getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames,
                                  obsWeights = obsWeights, control = control, verbose = verbose, 
                                  errorsInLibrary = errorsInCVLibrary)
    coef <- getCoef$coef
    names(coef) <- libraryNames
    
    if (!("optimizer" %in% names(getCoef))) {
        getCoef["optimizer"] <- NA
    }
    m <- dim(newX)[1L]
    predY <- matrix(NA, nrow = m, ncol = k)
    .screenFun <- function(fun, list) {
        screen_fn = get(fun, envir = env)
        testScreen <- try(do.call(screen_fn, list))
        if (inherits(testScreen, "try-error")) {
            warning(paste("replacing failed screening algorithm,", 
                          fun, ", with All() in full data", "\n "))
            out <- rep(TRUE, ncol(list$X))
        }
        else {
            out <- testScreen
        }
        return(out)
    }
    
    whichScreen <- sapply(library$screenAlgorithm, FUN = .screenFun, 
                          list = list(Y = Y, X = X, family = family, id = id, obsWeights = obsWeights), 
                          simplify = FALSE)
    whichScreen <- do.call(rbind, whichScreen)
    .predFun <- function(index, lib, Y, dataX, newX, whichScreen, 
                         family, id, obsWeights, verbose, control, libraryNames) {
        pred_fn = get(lib$predAlgorithm[index], envir = env)
        testAlg <- try(do.call(pred_fn, list(Y = Y, X = subset(dataX, 
                                                               select = whichScreen[lib$rowScreen[index], ], drop = FALSE), 
                                             newX = subset(newX, select = whichScreen[lib$rowScreen[index], 
                                                                                      ], drop = FALSE), family = family, id = id, obsWeights = obsWeights)))
        if (inherits(testAlg, "try-error")) {
            warning(paste("Error in algorithm", lib$predAlgorithm[index], 
                          " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
            out <- rep.int(NA, times = nrow(newX))
        }
        else {
            out <- testAlg$pred
            if (control$saveFitLibrary) {
                eval(bquote(fitLibrary[[.(index)]] <- .(testAlg$fit)), 
                     envir = fitLibEnv)
            }
        }
        if (verbose) {
            message(paste("full", libraryNames[index]))
        }
        invisible(out)
    }
    predY <- do.call("cbind", lapply(seq(k), FUN = .predFun, 
                                     lib = library$library, Y = Y, dataX = X, newX = newX, 
                                     whichScreen = whichScreen, family = family, id = id, 
                                     obsWeights = obsWeights, verbose = verbose, control = control, 
                                     libraryNames = libraryNames))
    errorsInLibrary <- apply(predY, 2, function(algorithm) anyNA(algorithm))
    if (sum(errorsInLibrary) > 0) {
        if (sum(coef[as.logical(errorsInLibrary)]) > 0) {
            warning(paste0("Re-running estimation of coefficients removing failed algorithm(s)\n", 
                           "Original coefficients are: \n", paste(coef, 
                                                                  collapse = ", "), "\n"))
            Z[, as.logical(errorsInLibrary)] <- 0
            if (all(Z == 0)) {
                stop("All algorithms dropped from library")
            }
            getCoef <- method$computeCoef(Z = Z, Y = Y, libraryNames = libraryNames, 
                                          obsWeights = obsWeights, control = control, verbose = verbose, 
                                          errorsInLibrary = errorsInLibrary)
            coef <- getCoef$coef
            names(coef) <- libraryNames
        }
        else {
            warning("Coefficients already 0 for all failed algorithm(s)")
        }
    }
    getPred <- method$computePred(predY = predY, coef = coef, 
                                  control = control)
    
    
    
    if (control$saveCVFitLibrary) {
        cvFitLibrary <- lapply(crossValFUN_out, "[[", "model_out")
    } else {
        cvFitLibrary <- NULL
    }
    colnames(predY) <- libraryNames
    if (sum(errorsInCVLibrary) > 0) {
        getCoef$cvRisk[as.logical(errorsInCVLibrary)] <- NA
    }
    
    out <- list(call = call, libraryNames = libraryNames, SL.library = library, 
                SL.predict = getPred, coef = coef, library.predict = predY, 
                Z = Z, cvRisk = getCoef$cvRisk, family = family, fitLibrary = get("fitLibrary", 
                                                                                  envir = fitLibEnv), cvFitLibrary = cvFitLibrary, 
                varNames = varNames, validRows = validRows, method = method, 
                whichScreen = whichScreen, control = control, cvControl = cvControl, 
                errorsInCVLibrary = errorsInCVLibrary, errorsInLibrary = errorsInLibrary, 
                metaOptimizer = getCoef$optimizer, env = env, times = times)
    class(out) <- c("SuperLearner")
    return(out)
}



