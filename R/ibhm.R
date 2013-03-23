require(DEoptim)

source('R/scals.R')
source('R/activs.R')
source('R/optim.R')
source('R/stopCriteria.R')
source('R/config.R')

# IBHM method -------------------------------------------------------------

TrainIBHM <- function(x,y, config = ConfigureIBHM()  ){
  if(config$jit){
    require(compiler)
    compiler::enableJIT(3)
  }
  
  ctx <- CreateContext(x,y,config)
  
  ctx$m$w0 <- mean(y)  
  
  while(ctx$continue){
    ctx$iteration <- ctx$iteration + 1
    
    ctx <- ctx$sc$eval( DoIteration(ctx) )    
  }  
    
  
  if(config$final.estimation == 'weights'){
    ctx <- OptimizeAllWeights(ctx)
  }else{
    ctx <- OptimizeAllParams(ctx)
  }
  
  ctx$log('\nFinal RMSE: ',CalcRMSE(ctx))
  return(ctx$m)
}


DoIteration <- function(ctx){
  ctx$log('\nIteration ',ctx$iteration, ' RMSE: ', CalcRMSE(ctx))    
  
  ctx <- FindActiv( FindScal(ctx) ) 
  ctx$m$w <- append(ctx$m$w,1)  
  OptimizeWeights(ctx)
}

CalcRMSE <- function(ctx){
  sqrt(mean((ctx$y-predict(ctx$m,ctx$x))^2))
}




# Method runtime context --------------------------------------------------

CreateContext <- function(x,y, params){
  x <- as.matrix(x)
  y <- as.matrix(y)
  stopifnot(nrow(x)==nrow(y), nrow(x)>20)  
  
  if(params$verbose){
    log <- function(...) cat(...,'\n')
  }else{
    log <- function(...) {}
  }
  
  
  append(params,
         list( x = x,
               y = y,               
               m=CreateEmptyModel(x,y), 
               scal.y = list(),
               continue=TRUE, 
               iteration=0,
               log = log
         )
  )
}


# IBHM model structure ----------------------------------------------------
CreateEmptyModel <- function(x=NULL,y=NULL){
  m <- list(w0=0,
            w=vector(),
            scals=list(),
            d=NULL, 
            activs=list(),
            a=vector(),
            b=vector(),
            train=list(x=x,y=y))
  class(m) <- 'IBHM'
  return(m)
}



# Standard S3 methods -----------------------------------------------------
predict.IBHM <- function(object,x = NULL,...){  
  m<-object
  if(is.null(x)){
    x <- m$train$x
  }else{
    x <- as.matrix(x)  
  }
  
  
  d <- dim(x)
  
  y <- matrix(m$w0,d[[1]],1)
  
  len <- length(m$scals)
  if(len > 0){
    for(i in 1:len){    
      act <- m$activs[[i]]
      a <- m$a[[i]]
      b <- m$b[[i]]
      
      scal <- m$scals[[i]]
      d <- m$d[,i]
      
      w <- m$w[[i]]    
      y <- y + w * act( a*scal(x,d)+b)
    }
  }
  
  y
}


length.IBHM <- function(x){
  length(x$w)
}

ToString <- function(m){  
  f <- function(v){
    sprintf('%.2e',v)
  }
  
  str <- f(m$w0)
  len <- length(m$w)
  
  if(len>0){  
    for(i in 1:length(m)){
      w <- m$w[[i]]
      if(w>0) s<-'+' else s<-''
      
      act <- expr(m$activs[[i]])
      a <- f(m$a[[i]])
      b <- f(m$b[[i]])
      
      scal <- expr(m$scals[[i]])
      d <- paste(f(m$d[,i]),collapse=' ')
      
      
      str <- paste(str,s,f(w),act,' (',a,'*',scal,'(x,[',d,'])','+',b,') ')
    }
  }
  
  str
  
}

summary.IBHM <- function(object,...){  
  y <- object$train$y
  x <- object$train$x
  
  se <- c((y - predict(object))^2)
  
  res <- list(model = ToString(object),
              model.size = length(object),
              TrainSize = nrow(y),
              TrainDim = ncol(x),
              mse = mean(se),
              sd.se = sd(se),
              rmse = sqrt(mean(se)),
              cor=cor(predict(object),y))
  class(res) <- 'summary.IBHM'
  res
}

print.summary.IBHM <- function(x,...){
  cat('Model equation: ', x$model,'\n',
      'Model size: ',x$model.size,'\n\n',
      'Train set dim: ', x$TrainDim, ' Train set size: ', x$TrainSize,'\n',      
      'MSE:  ', x$mse, ' Std. dev:', x$se.sd,'\n',
      'RMSE: ', x$rmse, '\n',
      'Pearson correlation coefficient: ', x$cor,'\n')
}


# Finding a scalarization function ----------------------------------------

FindScal <- function(ctx){
  ctx$log(' Finding next scal')
  
  nScals <- length(ctx$scal$candidates)  
  candidates <- vector("list",length=nScals)
  
  for( i in  1:nScals){
    scal <- ctx$scal$candidates[[i]]    
    candidates[[i]] <- OptimizeScal(scal,ctx)
    ctx$log('  ',expr(scal), ' eval : ', candidates[[i]]$eval)
  }
  
  bestIdx <- which.min( sapply(candidates, function(candidate){return(candidate$eval)}))  
  best <- candidates[[bestIdx]]
  ctx$log(' Best scal: ',expr(best$scal))
  
  k <- ctx$iteration
  ctx$m$scals[[k]]  <- best$scal
  ctx$m$d <- cbind(ctx$m$d,best$d)
  ctx$w <- ctx$wf(best$scal(ctx$x, best$d))
  ctx$scal.y[[k]] <- best$scal(ctx$x, best$d)
  
  ctx
}


OptimizeScal <- function(scal, ctx){
  
  goal <- function(par){              
    yh <- scal(ctx$x, par)
    
    if(var(yh)>0){
      val<-(1-abs(ctx$scal$eval(ctx$y,yh,ctx)))
    }else{
      val<-(Inf)
    }
    
    return(val)
  }
  n <- dim(ctx$x)[[2]]  
  res <- ctx$scal$optim(goal,n,ctx$scal$optim.params)    
  
  append(res, list(scal = scal))
}

EvaluateScal <- function(y,yh,ctx){
  w <- ctx$wf(yh)
  if(weighted.var(y,w) > 1e-4 && weighted.var(yh,w) > 1e-4){
    cor.y <- weighted.rho(y, yh,w)    
    indep <- 1.0
    for(ys in ctx$scal.y){    
      indep <- indep * (1-abs(weighted.r(yh,ys,w)))
    }
    
    if(is.nan(cor.y)){
      cor.y<-0
    }
    if(is.nan(indep)){
      indep <- 1.0  
    }
    
    return(abs(cor.y) * indep)
  }else{
    return(0.0)
  }
}

# Finding an activation function #####################


FindActiv <- function(ctx){
  ctx$log(' Finding next activ')
  nActivs <- length(ctx$activ$candidates)  
  candidates <- vector("list",length=nActivs)
  
  yh <- RunLastScal(ctx)
  
  for( i in  1:nActivs){
    activ <- ctx$activ$candidates[[i]]
    
    candidates[[i]] <- OptimizeActiv(activ, yh, ctx)
    ctx$log('  ',expr(activ),' eval:',candidates[[i]]$eval)
  }
  
  bestIdx <- which.min( sapply(candidates, function(candidate){return(candidate$eval)}))  
  best <- candidates[[bestIdx]]
  ctx$log(' Best activ: ', expr(best$activ))
  
  k <- ctx$iteration
  ctx$m$activs[[k]] <- best$activ
  ctx$m$a[[k]]      <- best$a
  ctx$m$b[[k]]      <- best$b
  
  ctx
}


OptimizeActiv <- function(activ, yh, ctx){  
  goal <- function(par){          
    a <- par[[1]]
    b <- par[[2]]
    
    ya <- activ(a*yh+b)    
    if(var(ya) > 0){
      val <- 1-abs(ctx$activ$eval(ctx$y,ya,ctx$w))
    }else{
      val <- Inf
    }
    
    if(is.nan(val)){
      val <- Inf
    }
    return(val)
  }  
  
  res <- ctx$activ$optim(goal, ctx$activ$optim.params)  
  
  append(res, list(activ=activ))
}

EvaluateActiv <- function(y,ya,w){
  val<-weighted.r(y, ya,w)  
  if(is.nan(val)){
    val<-0
  }
  return(val)
}

RunLastScal <-function(ctx){
  lastScalIdx <- length(ctx$m$scals)
  scal <- ctx$m$scals[[lastScalIdx]]  
  d <- ctx$m$d[,lastScalIdx]
  
  scal(ctx$x,d)  
}


# Estimating regression weights -------------------------------------------

OptimizeWeights <- function(ctx){  
  m <- ctx$m
  goal <- function(w){        
    m$w <- w
    
    return(weighted.mean((ctx$y - predict(m,ctx$x))^2,ctx$w))    
  }
  
  res <- optim(m$w,goal, method="BFGS")  
  m$w <- res$par
  
  ctx$m <-m
  
  ctx
}

OptimizeAllWeights <- function(ctx){
  ctx$log('Final weight optimization ')
  
  m <- ctx$m
  goal <- function(par){    
    m$w0 <- par[[1]]
    m$w <- par[2:length(par)]
    
    return(mean((ctx$y - predict(m,ctx$x))^2))
  }
  
  res <- optim(c(m$w0,m$w),goal, method="BFGS")
  m$w0 <- res$par[[1]]
  m$w <- res$par[2:length(res$par)]
  
  ctx$m <-m
  
  ctx
}

OptimizeAllParams <- function(ctx){
  ctx$log('Final parameter reoptimization')
  m <- ctx$m
  goal <- function(par){    
    m <- extract(par,m)
    
    return(mean((ctx$y - predict(m,ctx$x))^2))
  }
  
  mask <- m
  mask$w0 <- 1
  mask$w[] <- 2
  mask$a[] <- 3
  mask$b[] <- 4
  mask$d[,] <- 5
  mask <- c(mask$w0, mask$w, mask$a, mask$b, as.vector(mask$d))
  
  
  extract <- function(par,m){    
    s <- split(par, mask)
    m$w0 <- s[[1]]    
    m$w <- s[[2]]
    m$a <- s[[3]]
    m$b <- s[[4]]
    m$d <- matrix(s[[5]], nrow=nrow(m$d))
    m
  }
  
  res <- optim(c(m$w0,m$w,m$a,m$b,as.vector(m$d)),goal, method="BFGS")  
  ctx$m <- extract(res$par,m)
  
  ctx
}






