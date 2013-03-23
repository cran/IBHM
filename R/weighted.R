# Weighted mean
weighted.mean<-function(y,w,unif=FALSE){  
  s <- sum(w)
  if(unif || s < 1e-10){
    return(sum(y*w))
  }else{
    return(sum(y*w)/s)      
  }
}

# Weighted variance
weighted.var<-function(y,w, y.mean=weighted.mean(y,w,unif), unif=TRUE){
  weighted.mean((y-y.mean)^2,w, unif)
}

# Weighted Spearman correlation coefficient
weighted.rho<-function(y1,y2,w){
  r1 <- rank(y1,'average')
  r2 <- rank(y2,'average')
  
  weighted.r(r1,r2,w)
}

# Weighted Pearson correlation coefficient
weighted.r<-function(y1,y2,w){
  w <- w/sum(w)
  wm.y1 <- weighted.mean(y1,w,unif=TRUE)
  wm.y2 <- weighted.mean(y2,w,unif=TRUE)
  
  return(weighted.mean((y1-wm.y1)*(y2-wm.y2),w,unif=TRUE)/sqrt(weighted.var(y1,w,wm.y1,unif=TRUE)*weighted.var(y2,w,wm.y2,unif=TRUE)))
}

# Weighted cov  coefficient
weighted.cov<-function(y1,y2,w){
  w <- w/sum(w)
  wm.y1 <- weighted.mean(y1,w,unif=TRUE)
  wm.y2 <- weighted.mean(y2,w,unif=TRUE)
  
  return(weighted.mean((y1-wm.y1)*(y2-wm.y2),w,unif=TRUE))
}

