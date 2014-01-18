\name{predict.IBHM}
\alias{predict.IBHM}
\alias{predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict.IBHM
}
\description{
Predicted values based on IBHM model object.
}
\usage{
\method{predict}{IBHM}(object, x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of class inheriting from "ibhm"
}
  \item{x}{
  An object for which the predictions will be calculated - should be convertible to a \code{matrix}. If 
  not given the values from the training data are used.
}
  \item{\dots}{
    Further arguments.
}
}

\value{
\code{predict.IBHM} returns a vector of predicted values with length corresponding to the number of rows of the \code{x} parameter.
}

\seealso{
\code{\link{TrainIBHM}}
}
\examples{
  x.train <-seq(-2,2,length.out=100)
  y.train <-tanh(x.train)
  
  m <- IBHM::TrainIBHM(x.train,y.train)  
  cat('Train MSE:',mean((y.train-predict(m))^2),'\n')
  
  x.test <- runif(100, min=-4, max=4)
  yh.test <- predict(m, x.test)
    
  plot(x.test, yh.test)
}
\keyword{ ~models }
\keyword{ ~regression }
\keyword{ ~nonlinear }


