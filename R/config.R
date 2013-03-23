
# IBHM configuration ------------------------------------------------------

ConfigureIBHM  <- function(stop.criterion = IterationSC(3),
                           weighting.function = function(y){0.01+dcauchy(y)},                                                      
                           scal.optim = 'multi.CMAES',
                           scal.optim.params = list(retries=3, inner=list(maxit=50, stopfitness=-1)),
                           scal.candidates = c('dot.pr','radial','root.radial'),                                                      
                           activ.optim = 'multi.CMAES',
                           activ.optim.params = list(retries=3, inner=list(maxit=100, stopfitness=-1)),
                           activ.candidates = c('tanh','logsig','lin'),
                           jit=TRUE,
                           verbose=FALSE
){  
  require(DEoptim)
  list( sc = stop.criterion,
        wf = weighting.function,
        final.estimation = 'weights',
        jit=jit,
        verbose=verbose,
        scal = list(
          eval = EvaluateScal,          
          optim = switch(scal.optim, 
                         CMAES = OptimizeScalCMAES,
                         multi.CMAES = OptimizeScalMultiCMAES,
                         DE = OptimizeScalDE,
                         multi.DE = OptimizeScalMultiDE,
                         stop('Unknown scal optimization method :',activ.optim)),
          optim.params = scal.optim.params,
          candidates = ScalFunctions(scal.candidates)
        ),
        activ = list(
          eval = EvaluateActiv,
          optim = switch(activ.optim, 
                         CMAES = OptimizeActivCMAES,
                         multi.CMAES = OptimizeActivMultiCMAES,
                         DE = OptimizeActivDE,
                         multi.DE = OptimizeActivMultiDE,
                         stop('Unknown activ optimization method :',activ.optim)),
          optim.params = activ.optim.params,
          candidates = ActivationCandidates(activ.candidates)
        )
  )
}