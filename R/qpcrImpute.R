# actual function
qpcrImpute <- function(object, dj=NULL, pyfit=NULL, groupVars=NULL, batch=NULL, tol=1,
                       iterMax=100, outform=c("Single","Param","Multy"), 
                       vary_fit=TRUE, vary_model=TRUE, add_noise=TRUE, 
                       formula=NULL, numsam=5, linkglm = c("logit", "probit", "cloglog"))
{
  linkglm <- match.arg(linkglm)
  outform <- match.arg(outform)
  ## check input
  if(class(object)!="qPCRset") stop("Input data must be of class qPCRset.")
  if(sum(featureCategory(object) == "Undetermined") == 0){
    stop("No values in featureCategory(object) are Undetermined.
         There is no need to impute missing values.")
  }
  if(!identical(colnames(exprs(object)),
                as.character(pData(object)$sampleName))){
    stop("colnames(exprs(object)) and pData(object)$sampleName
         must be identical.")
  }

  if(outform=="Multy" & ((vary_fit+vary_model+add_noise)==0)){
    warning("vary_fit, vary_model, and add_noise are all FALSE. Performing single imputation.")
    outform <- "Imput"
  }

  ## select normalization genes
  if(is.null(dj)){
    i.hk <- grep("control", featureType(object), ignore.case=TRUE)
    if(length(i.hk) == 0){
      message("No control genes found in featureType(object).
              Data will not be normalized.")
      dj <- rep(0, ncol(object))
    } else {
      if(length(i.hk) == 1){
        dj <- exprs(object)[i.hk,]
        dj <- dj-mean(dj)
      } else{
        dj <- colMeans(exprs(object)[i.hk,])
        dj <- dj-mean(dj)
      }
    }
  }
  ## accounting for batch effect 
  if (!is.null(batch)) 
     {
      B <- batch
     } else { 
      B <- matrix(0, nrow=nrow(exprs(object)), ncol=ncol(exprs(object)))
  }
  
  
  ## select target genes
  ind <- grep("target", featureType(object), ignore.case=TRUE)
  if(length(ind) == 0) stop("No target genes found in featureType(object).")

  ## sample grouping variable(s)
  if(is.null(formula)){
    if(is.null(groupVars)){
      nrep <- apply(pData(object), 1, paste, collapse=":")
    } else {
      if(length(groupVars)>1){
        nrep <- apply(pData(object)[,groupVars], 1, paste, collapse=":")
      } else nrep <- pData(object)[,groupVars]
    }
    if(class(nrep)!="character"){
    nrep <- droplevels(nrep)
    }
    if(length(levels(nrep))==1) formula=as.formula("~1") else formula=as.formula("~0+nrep")
    DesLM=model.matrix(formula,as.data.frame(nrep))
  } else{
    formula <- as.formula(formula)
    DesLM=model.matrix(formula, data=pData(object))
  }
  print(formula)
  
  ## take out batch effect
  Ct.nobatch <- exprs(object)-B
  Y          <- sweep(Ct.nobatch[ind,], 2, dj)
  tst        =  lmFit(Y, design=DesLM)

  repName <- nrep
  nrep <- as.numeric(as.factor(nrep))

  ## initial fit
  if(is.null(pyfit)){
    if(length(nrep)>1){
    
      p.nd <- as.vector(
        apply(featureCategory(object)[ind,], 1,
              function(x) by(x == "Undetermined", nrep, mean)))
      zs <- as.vector(
        apply(exprs(object)[ind,], 1, function(x) by(x, nrep, median)))
      ws <- as.vector(
        apply(featureCategory(object)[ind,], 1,
              function(x) by(x, nrep, length)))
    } else {
        p.nd <- as.vector(
          apply(featureCategory(object)[ind,], 1,
                function(x) mean(x == "Undetermined")))
        zs <- as.vector(
          apply(exprs(object)[ind,], 1, median))
        ws <- as.vector(
          apply(featureCategory(object)[ind,], 1, length))
      }
   # pyfit <- glm(p.nd~zs, family=quasibinomial(link=linkglm), weights=ws) # GLM
    pyfit <- bayesglm(p.nd~zs, family=binomial(link=linkglm), weights=ws)# Bayes GLM arm package
  }
  
  ## format data for EM algorithm
  Ct    <- as.vector(exprs(object)[ind,])
  B.vec <- as.vector(B[ind,])
  ngene <- as.numeric(as.factor(rep((1:nrow(object))[ind], ncol(object))))
  i.nd  <- as.vector(featureCategory(object)[ind,] == "Undetermined")
  dj    <- rep(dj, each=length(ind))
  ntype <- rep(nrep, each=length(ind))

  ## initial means
  thetaVec <- vector(length=length(Ct))
  thetaMat <- matrix(nrow=max(ngene), ncol=max(ntype))

  thetaVec.tst <- vector(length=length(Ct))
  thetaMat.tst <- matrix(nrow=max(ngene), ncol=max(ntype))
  #####
  for(i in 1:max(ngene)){
    for(j in 1:max(ntype)){
      ind2 <- intersect(which(ngene == i), which(ntype == j))
      ###
      thetaVec[ind2] <- thetaMat[i,j] <- mean(Ct[ind2]-dj[ind2]-B.vec[ind2])
      thetaVec.tst[ind2]<-thetaMat.tst[i,j]<-tst$coefficients[i,j]
    }
  }

  ## initial variance
  s2Vec <- vector(length=length(Ct))
  s2Mat <- matrix(nrow=max(ngene), ncol=max(ntype))

  s2Vec.tst <- vector(length=length(Ct))
  Name.s2<-s2Mat.tst <- matrix(nrow=max(ngene), ncol=max(ntype))
  #####
  for(i in 1:max(ngene)){
    for(j in 1:max(ntype)){
      ind2            <- which(ngene == i)
  ###
      s2Vec[ind2]     <- s2Mat[i,j] <- sum((Ct[ind2]-thetaVec[ind2]-dj[ind2]-B.vec[ind2])^2)/(length(ind2)-max(ntype))
      s2Vec.tst[ind2] <- s2Mat.tst[i,j] <- tst$sigma[i]^2
      Name.s2[ind2]   <- names(tst$sigma[i])
    }
  }
  names(s2Vec.tst) <- Name.s2
  
  ## begin EM algorithm
  params.new <- params.old <- list(
                                    s2Vec=s2Vec.tst,
                                    s2Mat=s2Mat.tst,
                                    thetaVec=thetaVec.tst,
                                    thetaMat=thetaMat.tst,
                                    sigma=tst$sigma,
                                    cov.matrix=tst$cov.coefficients
                                  )
  ll <- vector(length=iterMax)
  iter <- 1
  cond <- TRUE

  ## EM algorithm
  while(cond){

    print(paste(iter, "/", iterMax))
    params.old <- params.new

    ## E-Step
    ## missing values
    ez <- rep(NA, length=length(Ct))
    for(i in 1:length(Ct)){
      if(i.nd[i]){
        #####
        xi <- params.new$thetaVec[i]+dj[i]+B.vec[i]
        ez[i] <- EZ(xi, params.new$s2Vec[i], pyfit)
      }
    }

    ## M-Step -- theta
    thetas <- updateTheta(Y, ez, dj, B.vec, i.nd, ngene, ntype, DesLM)
    params.new$thetaVec   <- thetas[[1]]
    params.new$thetaMat   <- thetas[[2]]
    params.new$sigma      <- thetas[[3]]
    params.new$cov.matrix <- thetas[[4]]

    ## E-Step
    ## missing values
    ez <- ez2 <- rep(NA, length=length(Ct))
    for(i in 1:length(Ct)){
      if(i.nd[i]){
        #####
        xi <- params.new$thetaVec[i]+dj[i]+B.vec[i]
        ez[i] <- EZ(xi, params.new$s2Vec[i], pyfit)
        ez2[i] <- EZ2(xi, params.new$s2Vec[i], pyfit, ez[i])
      }
    }

    ## M-Step -- sigma2
    sigmas <- updateS2(Ct, params.new$thetaVec, dj, B.vec, ez, ez2,
                       i.nd, ngene, ntype)
    params.new$s2Vec <- sigmas[[1]]
    params.new$s2Mat <- sigmas[[2]]

    #####
    ll[iter] <- logLik(Ct, ez, ez2, params.new$s2Vec,
                       params.new$thetaVec, dj, B.vec, i.nd)
    message(ll[iter])

    if(iter>1) cond <- (abs(ll[iter]-ll[iter-1]) > tol) & (iter < iterMax)
    iter <- iter+1

    ## update fit
    tmp <- Ct
    tmp[which(as.logical(i.nd))] <- ez[which(as.logical(i.nd))]
    zs <- tmp
    #gavg <- as.vector(by(tmp, paste(ntype, ngene, sep=":"), median))
    #p.nd <- as.vector(by(i.nd, paste(ntype, ngene, sep=":"), mean))
    #ws <- as.vector(by(i.nd, paste(ntype, ngene, sep=":"), length))
    #pyfit <- glm(p.nd~gavg, family=binomial(link=logit), weights=ws)
    
    #pyfit <- glm(as.vector(i.nd)~zs, family=binomial(link=logit))
    pyfit <- bayesglm(as.vector(i.nd)~zs, family=binomial(link=linkglm))# Bayes GLM arm package
  }

  ############# if statement about returning parameters or values
  print(outform)
  #####
  if (outform == "Multy") {
    multylist<-multy(object, pyfit, numsam, params.new, Ct, Y, dj, B.vec, ez, ez2,
                     i.nd, ngene, ntype, DesLM, iterMax, tol,
                     vary_fit, vary_model, add_noise)
  } 
  else {
    if (outform == "Single") {
      Ct[which(as.logical(i.nd))] <- ez[which(as.logical(i.nd))]
      ind <- grep("target", featureType(object), ignore.case=TRUE)
      exprs(object)[ind,] <- Ct

      fc <- as.matrix(featureCategory(object))
      fc[which(fc == "Undetermined", arr.ind=TRUE)] <- "Imputed"
      featureCategory(object) <- as.data.frame(fc)

      if (nrow(getCtHistory(object)) == 0){
        setCtHistory(object) <- data.frame(
          history = "Manually created qPCRset object.",
          stringsAsFactors = FALSE)
      }
      setCtHistory(object) <- rbind(getCtHistory(object),
                                    capture.output(match.call(qpcrImpute)))

      object
    }
    else {
      if (outform == "Param") {
        sigma2<-params.new$s2Mat.tst
        Parameters <- list("sigma2"=params.new$s2Mat,
                           "theta"=params.new$thetaMat,
                           "formula"=formula, "tst"=tst,
                           "fit"=pyfit) # added fit to the output
        colnames(Parameters$sigma2) <- colnames(Parameters$theta) <- substring(colnames(Parameters$tst$design),5)
        rownames(Parameters$sigma2) <- rownames(Parameters$theta) <- names(Parameters$tst$Amean)
        message("Parameters are saved in the list")
        return(Parameters)}
    }

  }
  }
