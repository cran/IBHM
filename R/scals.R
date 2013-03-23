# Scalarization functions ##############################

ScalFunctions <- function(filter = NULL){
  # Dot product scalarization
  DotPr <- function(x,d){  
    return(x%*%d)
  }
  expr(DotPr) <- 'dot.pr'



  # Radial scalarization
  Radial <- function(x,d){
    n <- dim(x)[[2]]
    m <- dim(x)[[1]]
    
    d_m <- matrix(rep(d,m),nrow=m,ncol=n)
    
    rowSums((x-d_m)^2)
  }
  expr(Radial) <- 'radial'

  # Radial scalarization - euclidean distance
  RootRadial <- function(x,d){
    n <- dim(x)[[2]]
    m <- dim(x)[[1]]
    
    d_m <- matrix(rep(d,m),nrow=m,ncol=n)
    
    rowSums(sqrt((x-d_m)^2))
  }
  expr(RootRadial) <- 'root.radial'
  
  res <- list(DotPr, Radial, RootRadial)  
  if(!is.null(filter)){
    res <- Filter(function(v) any(expr(v)==filter), res)
  }
  if(length(res)==0){
    stop('Invalid scal list given : ', paste(filter,sep=', '))
  }
  
  res
}