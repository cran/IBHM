\name{length.IBHM}
\alias{length.IBHM}
\alias{length}
\title{
length.IBHM
}
\description{
Returns the number of components in a IBHM model.
}
\usage{
\method{length}{IBHM}(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  An object inheriting from class "IBHM". 
}
}
\value{
The number of components in the given IBHM model.
}

\seealso{
\code{\link{TrainIBHM}}
}
\examples{
  x.train <-seq(-2,2,length.out=100)
  y.train <-tanh(x.train)
  
  m <- IBHM::TrainIBHM(x.train,y.train)  
  length(m)
}
\keyword{ ~models }
\keyword{ ~regression }
\keyword{ ~nonlinear }


