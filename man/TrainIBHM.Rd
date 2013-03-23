\name{TrainIBHM}
\alias{TrainIBHM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
TrainIBHM
}
\description{
Trains an IBHM model usinggiven training data with input variables \code{x} and a single output variable \code{y}.
}
\usage{
TrainIBHM(x, y, config = ConfigureIBHM())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Training data input variables. Should be convertible to a \code{matrix} with each row corresponding to a single data point.
}
  \item{y}{
  Training data output variable. Should be convertible to a \code{matrix} with a single column and the number of rows equal to the number of rows of \code{x}. 
}
  \item{config}{
  A configuration object created using \code{\link{ConfigureIBHM}}. Contains settings such as the stop criterion, optimization method parameters etc.
}
}

\value{
The created model which is an object of class "IBHM".
}

\references{
Zawistowski, P. and Arabas, J.: "Benchmarking IBHM method using NN3 competition dataset." In \emph{Proc. 6th int. conf.  on Hybrid artificial intelligent systems} - Vol. 1, HAIS'11, pp 263--270, 2011. Springer-Verlag.
}

\seealso{
\code{\link{ConfigureIBHM}},\code{\link{summary.IBHM}},\code{\link{predict.IBHM}}
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

