#' @useDynLib SmeiDMFMC
#' @importFrom Rcpp sourceCpp


CUM <- function(X,...){
  if(ncol(X) == 0){
    cum <- M4.MM(X)
  }else{
    n = length(X[1,])
    t = length(X[,1])
    m2 <- (t-1)/t*cov(X)
    cum<-M4.MM(X) - t(c(m2))%x%m2 - matrix((m2)%x%(m2),n,n^3) - m2%x%t(c(m2))}
  return(cum)
}

mm<-function(h,u,sigma,lambda,p,q){
  mm=vector()

  v=q^(-1/p)*((3*lambda^2+1)*(beta(3/p,q-2/p)/beta(1/p,q))-4*lambda^2*(beta(2/p,q-1/p)/beta(1/p,q))^2)^(-1/2)

  m=(2*v*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)

  for (r in 0:h) {
    mm[r+1]=choose(h,r)*((1+lambda)^(r+1)+(-1)^(r)*(1-lambda)^(r+1))*(-lambda)^(h-r)*((v*sigma)^h*q^(h/p)/(2^(r-h+1)))*
      (beta((r+1)/p,q-r/p)*beta(2/p,q-1/p)^(h-r)/(beta(1/p,q)^(h-r+1)))
  }
  return(sum(mm))
}


fnorm <- function(a){
  sqrt(sum((a)^2))
}

msst = function(h,u,sigma,l,q){
  mh = gamma((q-h)/2)*gamma((h+1)/2)*(q-2)^(h/2+0.5)/(sqrt(pi*(q-2))*gamma(q/2))*((l^(h+1)+(-1)^h/l^(h+1))/(l+1/l))
  return(mh)
}

msged = function(h,g,v){
  c = sqrt(2^(-2/v)*gamma(1/v)/gamma(3/v))
  mh = (2^(1/v)*c)^h*gamma((h+1)/v)/gamma(1/v)
  return(mh)
}

calmoments = function(skew, shape, dist = "sstd"){
  ###  the first four moments
  n = max(length(shape),length(skew))
  mi = matrix(NA, nrow = n, ncol = 4)
  r = 1:4
  ### compute the absolute of the distribution
  if(dist == "sstd") {
    sstd.m = function(shape){
      r = 1:4
      num1 = gamma((shape - r)/2)*gamma((r + 1)/2)*(shape - 2)^((r + 1)/2)
      num2 = sqrt(pi*(shape - 2))*gamma(shape/2)
      return( num1 / num2)
    }
    mi = t(mapply(FUN = sstd.m , shape = shape))
  }
  if(dist == "sged"){
    sged.m = function(shape){
      r = 1:4
      lambda = sqrt(2^(-2/shape) * gamma(1/shape) / gamma(3/shape))
      num = (2^(1/shape)*lambda)^r *gamma((r + 1)/shape) / gamma(1/shape)
      return(num)
    }
    mi = t(mapply(FUN = sged.m, shape = shape))
  }
  ### compute the scale factor
  scale = function(skew){
    r = 1:4
    num = (skew^(r + 1) + (-1)^r / skew^(r + 1)) / (skew + 1/skew)
    return(num)
  }
  scales = t(mapply(FUN = scale, skew = skew))
  Mi = scales*mi
  ### compute the skewness and kurtosis
  m1 = Mi[,1]
  m2 = Mi[,2]
  m3 = Mi[,3] - 3*Mi[,1]*Mi[,2] + 2*Mi[,1]^3
  m4 = Mi[,4] - 4*Mi[,1]*Mi[,3] + 6*Mi[,2]*Mi[,1]^2 - 3*Mi[,1]^4
  skewness = m3 / (m2 - m1^2)^(3/2)
  kurtosis = m4 / (m2 - m1^2)^2
  value = data.frame(skewness = skewness, kurtosis = kurtosis)
  return(value)
}

VAR = function(w,mu,m2,m3,m4,mv = F,MU = F,...){
  w = w/sum(w)
  if(mv == F){
    sig = sqrt(t(w)%*%m2%*%w)
    ske = t(w)%*%m3%*%(w%x%w)/sig^3
    kurt = t(w)%*%m4%*%(w%x%w%x%w)/(sig^4)
    z = qnorm(alpha)
    VAR = sig*(z + ske/6*(z^2-1)+kurt/24*(z^3-3*z)-ske^2/36*(2*z^3-5*z))
    if(MU == T){
      VAR = t(mu)%*%w + VAR
    }
  }
  if(mv == T){
    sig = sqrt(t(w)%*%m2%*%w)
    z = qnorm(alpha)
    VAR = t(mu)%*%w + sig*z
  }
  return(-VAR)
}


Uc=function(w,mu,m2,m3,m4,mv = F,gamma = gamma,MU = F,...){
  w = w/sum(w)
  n=length(w)
  
  if(mv == F){
    
    U = gamma/2*(t(w)%*%m2%*%w) + gamma*(gamma+1)/6*t(w)%*%m3%*%(w%x%w) - gamma*(gamma+1)*(gamma+2)/24*(t(w)%*%m4%*%(w%x%w%x%w))
    
    if(MU == T){
      U = t(w)%*%mu + U  
    }
    
  }
  if(mv == T){
    U = t(mu)%*%w - gamma/2*(t(w)%*%m2%*%w)
  }
  return(-U)
}

ES = function(w,mu,m2,m3,m4,mv = F,MU = F,...){
  w = w/sum(w)

  I = function(q,...){

    if(q %in% c(2,4,6,8)){
      ss = 0
      for (i in 1:(q/2)) {
        ss = ss +  (2^(q/2)*factorial(q/2))/(2^i*factorial(i))*g^(2*i)*phig
      }
      ss = ss + (2^(q/2)*factorial(q/2))*phig
    }
    if(q %in% c(1,3,5,7)){
      qs = (q-1)/2
      ss = 0
      for (i in 0:(qs)) {
        ss = ss + (factorial(2*qs+1)/((2^(qs)*factorial(qs))))/(factorial(2*i+1)/((2^(i)*factorial(i))))*g^(2*i+1)*phig
      }
      ss = ss - (factorial(2*qs+1)/((2^(qs)*factorial(qs))))*Phig
    }

    return(ss)
  }

  if(mv == F){
    sig = sqrt(t(w)%*%m2%*%w)
    ske = t(w)%*%m3%*%(w%x%w)/sig^3
    kurt = t(w)%*%m4%*%(w%x%w%x%w)/(sig^4)

    z = qnorm(alpha)
    P1 = (z^2-1)*ske/6
    P2 = (z^3-3*z)*kurt - (2*z^3-5*z)/36*ske^2
    G2 = z + P1 + P2

    g = G2
    phig = dnorm(g)
    Phig = pnorm(g)

    I1 = I(1,g,phig,Phig)
    I2 = I(2,g,phig,Phig)
    I3 = I(3,g,phig,Phig)
    I4 = I(4,g,phig,Phig)
    I6 = I(6,g,phig,Phig)

    EG = -1/alpha*(phig + kurt/24*(I4-6*I2+3*phig) + ske/6*(I3-3*I1) + (ske^2/72)*(I6-15*I4+45*I2-15*phig))
    MES = t(mu)%*%w + sig*EG

  }
  if(mv == T){
    sig = sqrt(t(w)%*%m2%*%w)
    ske = 0
    kurt = 0

    z = qnorm(alpha)
    phig = dnorm(z)

    MES = t(mu)%*%w - sig*phig/alpha
  }
  return(-MES)

}


INVCUM <- function(c4,m2,...){
  mm = c4 + t(c(m2))%x%m2 + matrix((m2)%x%(m2),n,n^3) + m2%x%t(c(m2))
  return(mm)
}

INVM <- function(m4,m2,...){
  cum = m4 - t(c(m2))%x%m2 - matrix((m2)%x%(m2),n,n^3) - m2%x%t(c(m2))
  return(cum)
}


Portfolio.Cumulants = function (w, mm_factor, mm_eps, A) 
{
  m2f = mm_factor[[1]]
  m3f = mm_factor[[2]]
  c4f = mm_factor[[3]]
  m2e = mm_eps[[1]]
  m3e = mm_eps[[2]]
  c4e = mm_eps[[3]]
  B = t(w) %*% A
  m2P = sum((B^2) * m2f) + sum((w^2) * m2e)
  m3P = sum((B^3) * m3f) + sum((w^3) * m3e)
  m4P = sum((B^4) * c4f) + sum((w^4) * c4e) + 3 * m2P^2
  mmP = c(m2P, m3P, m4P)
  return(mmP)
}

Portfolio.Cumulants.Mat = function (w, mm_factor, mm_eps, A){
    m2f = mm_factor[[1]]
    m3f = mm_factor[[2]]
    c4f = mm_factor[[3]]
    m2e = mm_eps[[1]]
    m3e = mm_eps[[2]]
    c4e = mm_eps[[3]]
    B = t(w) %*% A
    m2P = B %*% m2f %*% t(B)                     + sum((w^2) * m2e)
    m3P = B %*% m3f %*% (t(B) %x% t(B))          + sum((w^3) * m3e)
    m4P = B %*% c4f %*% (t(B) %x% t(B) %x% t(B)) + sum((w^4) * c4e) + 3 * m2P^2
    mmP = c(m2P, m3P, m4P)
    return(mmP)
  }

Portfolio_Weight = function(MM,opt.type = c("VAR","Uc","MV"),gamma = 10, 
                            alpha = 0.05, shortselling = T,lb = -1,ub = 2,MU = F){
  M1 = MM[[1]];M2 = MM[[2]];M3 = MM[[3]];C4 = MM[[4]];M4 = MM[[5]];
  if(shortselling == T){
    lb = lb
    ub = ub
  }else{
    lb = 0
    ub = 1
  }
  if(opt.type == "VAR"){
    ww =  Rsolnp::solnp(rep(1/n,n),fun = 
                          function(s,...){VAR(s,M1,M2,M3,C4,alpha = alpha,mv = F,MU = MU)},
                        eqfun = function(s){sum(s)},eqB = 1,
                        LB = rep(lb,n),UB = rep(ub,n))$pars
  }
  if(opt.type == "Uc"){
    ww =  Rsolnp::solnp(rep(1/n,n),fun = 
                          function(s,...){Uc(s,M1,M2,M3,M4,gamma = gamma,MU = MU)},
                        eqfun = function(s){sum(s)},eqB = 1,
                        LB = rep(lb,n),UB = rep(ub,n))$pars
  }
  if(opt.type == "MV"){
    ww =  Rsolnp::solnp(rep(1/n,n),fun = 
                          function(s,...){VAR(s,M1,M2,M3,C4,alpha = alpha,mv = T)},
                        eqfun = function(s){sum(s)},eqB = 1,
                        LB = rep(lb,n),UB = rep(ub,n))$pars
  }
  return(ww)
}

Mat.k = function(A,k,eps = 10^-6){
  ev = eigen(A)$values
  mark = which(ev > eps)
  ev = ev[mark]
  evc = as.matrix(eigen(A)$vectors)[,mark]
  
  Matk = evc%*%diag(ev^k)%*%t(evc)
  
  return(Matk)
}


SuperDiag = function(x,k){
  n = length(x)
  if(k==2){
    return(diag(x,n,n))
  }
  if(k == 3){
    mat = matrix(0,n,n^2)
    for (i in 1:n) {
      mat[i,i + (i-1)*n] <- x[i]
    }
    return(mat)
  }
  if(k == 4){
    mat = matrix(0,n,n^3)
    for (i in 1:n) {
      mat[i,i + (i-1)*n + (i-1)*n^2] <- x[i]
    }
    return(mat)
  }
}


Vm = function(X,k){
  X = scale(X,scale = F)
  return(colMeans(X^k))
}


calcSecs             = function(tt){round(as.numeric(Sys.time())-tt,2)}


rescale              = function(v){return(v/norm(as.matrix(v),"F"))}

