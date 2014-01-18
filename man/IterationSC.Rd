\name{IterationSC}
\alias{IterationSC}

\title{
IterationSC
}
\description{
An iteration counting stop criterion for \code{\link{TrainIBHM}} function creating IBHM approximation models. Should be passed to \code{\link{ConfigureIBHM}}while creating a configuration object.
}
\usage{
IterationSC(max.iterations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{max.iterations}{
  The number of iterations after which model construction should be stopped.
}
}

\seealso{
\code{\link{ValidationSC}},\code{\link{ConfigureIBHM}}, \code{\link{TrainIBHM}}
}
\examples{
  x.train <-seq(-2,2,length.out=100)
  y.train <-tanh(x.train)
  
  m <- IBHM::TrainIBHM(x.train,y.train, ConfigureIBHM(stop.criterion=IterationSC(2)))  
  length(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~regression }
\keyword{ ~nonlinear }


