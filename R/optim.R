MultiOptim <- function(retries, optim.fun, ...){
  best <- list(eval=Inf)
  for(i in 1:retries){
    res <- optim.fun(...)
    if(res$eval < best$eval){
      best <- res
    }
  }
  
  res
}

# Scalarizing function optimization ---------------------------------------

OptimizeScalMultiCMAES <- function(goal, n, optim.params){
  MultiOptim(optim.params$retries, OptimizeScalCMAES, 
                 goal,n, optim.params$inner)
}

OptimizeScalCMAES <- function(goal,n, optim.params){
  require(cmaes)
  res <- cmaes::cma_es(runif(n,min=-10,max=10),goal,lower=rep(-1e+1,n),upper=rep(1e+1,n), control=optim.params)                  
  
  list(eval=res$value, d=res$par)
}

OptimizeScalMultiDE <- function(goal, n, optim.params){
  MultiOptim(optim.params$retries, OptimizeScalDE, 
                 goal,n, optim.params$inner)
}

OptimizeScalDE <- function(goal,n, optim.params){
  require(DEoptim)
  res <- DEoptim::DEoptim(goal,lower=rep(-1e+1,n),upper=rep(1e+1,n), optim.params)                  
  
  list(eval=res$optim$bestval, d=res$optim$bestmem)
}



# Activation function optimization ----------------------------------------

OptimizeActivMultiCMAES <- function(goal, optim.params){
  MultiOptim(optim.params$retries, OptimizeActivCMAES, goal,optim.params$inner)
}

OptimizeActivCMAES<- function(goal, optim.params){
  require(cmaes)
  res <- cmaes::cma_es(runif(2,min=-10,max=10),goal,lower=rep(-1e+1,2),upper=rep(1e+1,2), control=optim.params)                  
  
  list( eval=res$value, 
        a=res$par[[1]],
        b=res$par[[2]]
  )
}
OptimizeActivDE<- function(goal, optim.params){
  require(DEoptim)
  res <- DEoptim::DEoptim(goal,lower=rep(-1.0e+1,2),upper=rep(1.0e+1,2), optim.params)                 
  
  list( eval=res$optim$bestval, 
        a=res$optim$bestmem[[1]],
        b=res$optim$bestmem[[2]]
  )
}

OptimizeActivMultiDE <- function(goal, optim.params){
  MultiOptim(optim.params$retries, OptimizeActivDE, goal,optim.params$inner)
}