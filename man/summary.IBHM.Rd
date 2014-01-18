\name{summary.IBHM}
\alias{summary.IBHM}
\alias{print.summary.IBHM}
\alias{summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
summary.IBHM
}
\description{
summary method for class "IBHM".
}
\usage{
\method{summary}{IBHM}(object, ...)

\method{print}{summary.IBHM}(x, ...)
}


\arguments{
  \item{object}{
    An object of class "IBHM", usually created using \code{\link{TrainIBHM}}.
}
  \item{x}{
    An object of class "summary.IBHM" created using \code{summary}.
  }
  \item{\dots}{
  Further arguments.
}
}

\value{
  An object of class 'summary.IBHM' containing the following fields:  
  \item{model}{Equation stating the obtained model in human readable form.}
  \item{model.size}{Number of components in the model.}
  \item{TrainSize}{Size of the training set used to construct the model.}
  \item{TrainDim}{Number of input attributes in the training set.}
  \item{mse}{Mean squared error.}
  \item{se.sd}{Standard deviation of the squared error.}
  \item{rmse}{Root mean squared error.}
  \item{cor}{Linear correlation between the actual and predicted values on the train data set.}
}

\seealso{
\code{\link{TrainIBHM}}
}
\examples{
  x.train <-seq(-2,2,length.out=100)
  y.train <-tanh(x.train)
  
  m <- TrainIBHM(x.train,y.train)  
  summary(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~models }
\keyword{ ~regression }
\keyword{ ~nonlinear }

