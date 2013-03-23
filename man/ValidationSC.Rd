\name{ValidationSC}
\alias{ValidationSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ValidationSC
}
\description{
A validation set based stop criterion for \code{\link{TrainIBHM}} function creating IBHM approximation models. Should be passed to \code{\link{ConfigureIBHM}}while creating a configuration object.

}
\usage{
ValidationSC(x, y)
}
\arguments{
  \item{x}{
  Validation set input arguments, should be convertible to a \code{matrix}.
}
  \item{y}{
  Validation set predicted argument, should be convertible to a single column \code{matrix}.
}
}

\details{
The cirterion is checked after each iteration and the current model is used to predict values on the validation data set. When the error increases in comparison to the previous iteration,
the construction process is stopped, and the changes in the model from the last iteration are undone.
}

\seealso{
\code{\link{IterationSC}},\code{\link{ConfigureIBHM}}, \code{\link{TrainIBHM}}
}
\examples{
# Training data
  x <- seq(-3,3,length.out=400)
  y <- tanh(x)
  
  # A held out validation set for the stop criterion
  x.val <- runif(50,min=-6,max=6)
  y.val <- tanh(x.val)
  
  
  # Training the model using the validation set to prevent overfitting
  m <- TrainIBHM(x,y, 
            ConfigureIBHM(stop.criterion = ValidationSC(x.val, y.val))
       )
                      
  summary(m)
  plot(y.val,predict(m,x.val),asp=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~regression }
\keyword{ ~nonlinear }


