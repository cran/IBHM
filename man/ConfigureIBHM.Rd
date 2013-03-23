\name{ConfigureIBHM}
\alias{ConfigureIBHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ConfigureIBHM
}
\description{
Creates a configuration object for \code{\link{TrainIBHM}} function which creates IBHM approximation
models.
}
\usage{

ConfigureIBHM(stop.criterion = IterationSC(3),
              weighting.function = function(y){0.01+dcauchy(y)},                                                      
              
              scal.optim = 'multi.CMAES',
              scal.optim.params = list(retries=3, inner=list(maxit=50, stopfitness=-1)),
              scal.candidates = c('dot.pr','radial','root.radial'),                                                      
              
              activ.optim = 'multi.CMAES',
              activ.optim.params = list(retries=3, inner=list(maxit=100, stopfitness=-1)),
              activ.candidates = c('tanh','logsig','lin'),
              
              jit=TRUE,
              verbose=FALSE
            )          
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stop.criterion}{
  The stop criterion for the model construction process. Possible values include objects created using \code{\link{IterationSC}} (the default) or \code{\link{ValidationSC}}. 
  
  For simplicity, the defalut value causes the method to construct a model with a fixed number of components (3), however it's actually best to   use \code{\link{ValidationSC}} criterion with a separate validation set.
}
  \item{weighting.function}{
  Definition of the weighting function used during model construction. This function puts emphasis on local features of the approximated function  (see details).
}
  \item{scal.optim}{
  The optimization method used to estimate scalarization functions' parameters . Possible values are: \code{CMAES}, \code{multi.CMAES}(default), \code{DE}, \code{multi.DE}. 
  
  The parameter values set the optimization methods as follows: 
  \code{CMAES} uses covariance matrix evolution strategy implemented in the \code{cmaes} package, 
  \code{multi.CMAES} is multistart \code{cmaes}, 
  \code{DE} uses differential evolution  implemented in \code{DEoptim} package,
  \code{multi.DE} is multistart \code{DE}.
}
  \item{scal.optim.params}{
  The parameters passed to the optimization method used to estimate scalarization functions' parameters.
  In case of \code{CMAES} this should be a list, in case of \code{DE} this should be a \code{DEoptim::DEoptim.control} object. 
  
  In case of multistart versions of optimization methods this parameter is a \code{list} with two fields:
  \code{retries} denoting the number of restarts and \code{inner} containing the parameters passed
  during each restart to the underlying optimisation method.
}
  \item{scal.candidates}{
  Candidate scalarization functions (see details).
}
  \item{activ.optim}{
  The optimization method used to estimate activation functions' parameters - see description of \code{scal.optim}.
}
  \item{activ.optim.params}{
  The parameters passed to the optimization method used to estimate activation functions' parameters - see descripton of the \code{scal.optim.params}.
  
}
  \item{activ.candidates}{
  Candidate activation functions  (see details).
}
  \item{jit}{
  Enables the just-in-time compilation feature provided by the \code{compiler} package.
  }
  
  \item{verbose}{
  Enables verbose output (disabled by default).
  }
}

\details{
The model constructed by IBHM has the following form:
 \deqn{ \hat{f}(x) = w_0 + \sum_{i=1}^m w_i g(a_i h(x,d_i)+b_i)\ ,}{%
 f(x) = w_0 + \sum w_i g(a_i h(x,d_i)+b_i) ,
 }
 where \eqn{h:R^n\rightarrow R}{h:R^n->R} is a scalarization function, \eqn{g:R\rightarrow R}{g:R->R} is an activation function,
 \eqn{d_i} is a parameter vector, and \eqn{a_i, b_i, w_i} are scalar parameters.
 
 The parameter estimation is based on optimizaing weighted correlation 
 measures between the model output and the approximation residual. This allows
 for an iterative model construction process which estimates both model structure
 and parameter values. For more details see [Zawistowski and Arabas].
}

\value{
A configuration object for \code{\link{TrainIBHM}}.
}

\references{
Zawistowski, P. and Arabas, J.: "Benchmarking IBHM method using NN3 competition dataset." In \emph{Proc. 6th int. conf.  on Hybrid artificial intelligent systems} - Vol. 1, HAIS'11, pp 263--270, 2011. Springer-Verlag.
}


\seealso{
 \code{\link{TrainIBHM}}, \code{\link{ValidationSC}},\code{\link{IterationSC}}
}


\examples{
  x <- seq(-3,3,length.out=400)
  y <- tanh(x)
  
  x.val <- runif(50,min=-6,max=6)
  y.val <- tanh(x.val)
  
  m <- TrainIBHM(x,y, ConfigureIBHM(  scal.candidates = 'dot.pr', 
                                      activ.candidates = 'tanh',
                                      stop.criterion = ValidationSC(x.val, y.val)))
                                      
  summary(m)
  plot(y.val,predict(m,x.val),asp=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~regression }
\keyword{ ~nonlinear }

