require(mixtools)
require(Matrix)


#'
#' @export
MDSP_gamma = function (x, m.iter, eps, center = 2){
  #x: numeric vector
  #m.iter: maximum number of iterations
  #eps: tolerance level for stopping rule
  #g.int=3*mean(x)

  if(length(unique(x))<2)
  {
    class.g=rep(1,length(x)); g.b=mean(x)
    # print("homogeneous!!! in MDSP Gamma")
  }else{
    k.int=kmeans(x, centers = min(center, length(unique(x))))
    g.int=k.int$centers[which.max(abs(k.int$centers))]*0.95
    iter=1
    g.a=g.int; g.b=g.int+1
    while (abs(g.b-g.a)>eps & iter < m.iter){
      class.g=as.numeric(abs(x-g.a)< abs(x))
      g.b=sum(class.g*x)/sum(class.g)
      iter=iter+1
      g.a=g.b
    }
    if(iter==m.iter) {print("WARNING: maxium iteration achieves")}
  }
  return(list(l.group=class.g, gamma=g.b))
}

#l.group: group indicator (non-zero subgroup =1)
# gamma: sub-homogeneous effect of non-zero subgroup



### Modified BIC function


#########################################
#'
#' @importFrom Matrix KhatriRao
#' @import mixtools
#' @importFrom magrittr %>%
#' @export
MDSP_l <- function(Y, X, Z, m_vec, l1, kp, beta.int=NULL, centers = NULL,
                   max.iter = 300, eps=0.01){
  #############################################
  #############  ADMMM Algorithm ##############
  #############################################

  ## Input data X (heterogeneous, nm by p.b),
  ## Z (homogeneous including intercept, nm by p.a),
  ## and Y (nm by 1): by default, the longitudinal id is rep(1:n, each=m)
  ## m: number of individual measurements
  ## l1: shrinkage parameter
  ## kp: ADMM step
  ## beta.int: initial of individualized parameters (np.b by 1)

  # n=length(Y)/m
  n <- length(m_vec)
  p.b = dim(X)[2]
  p.a = ifelse(is.null(dim(Z)[2]), 0, dim(Z)[2])
  # M <- diag(n)%x%matrix(1,1,m)
  M <- lapply(1:n, function(i){matrix(1, 1, m_vec[i])}) %>%
    do.call(bdiag, .) %>% as.matrix()
  # print(class(M))
  # print(class(t(X)))

  X.td <- cbind(Z, t(KhatriRao(M, t(X))))
  # print(class(X.td))

  ### set tunning hyperparameters and initial values
  ## hyper para
  kpp = kp #ADMM step
  lm.s = l1 #L1 regularization

  ## initials
  if (is.null(beta.int)){nu.int=rep(0,n*p.b)}else{nu.int=beta.int}
  lm.int=rep(0,n*p.b)

  lm.h.new=lm.int
  nu.h.new=nu.int

  a.h.new=rep(0,p.a)
  b.h.new=rep(0,n*p.b)

  theta.h=rep(0,n*p.b+p.a)
  theta.h.new=rep(1,n*p.b+p.a)

  iter=1
  error=rep(1,max.iter-1)

  if(is.null(centers)){ centers <- rep(2, p.b)}

  #NOTE: It's better to have a warm start

  ptm <- proc.time()
  dif <- 1
  ##### Main loop
  while(iter < max.iter & dif > eps){
    error[iter]=sqrt(mean((theta.h-theta.h.new)^2))

    ###
    lm.h=lm.h.new
    nu.h=nu.h.new

    a.h=a.h.new
    b.h=b.h.new
    theta.h=theta.h.new

    ############## theta updates#################
    theta.h.new=solve(t(X.td)%*%X.td +
    bdiag(list(matrix(0, p.a, p.a), kpp*diag(n*p.b)))) %*%
      (t(X.td)%*%Y + matrix(c(rep(0,p.a), as.numeric(kpp*nu.h-lm.h)), ncol=1))
    if(p.a>0){
      a.h.new=theta.h.new[1:p.a]
      b.h.new=theta.h.new[-(1:p.a)]
    } else {
      a.h.new <- numeric()
      b.h.new <- theta.h.new
    }

    ######### nu,gamma updates ###############
    B.h.new=matrix(b.h.new,p.b,n)
    Nu.h.new=matrix(nu.h,p.b,n)
    g.h=rep(0,p.b)
    g.h.new=rep(0,p.b)
    cl.h=matrix(0,p.b,n)
    Lm=matrix(lm.h,p.b,n)

    for (k in 1:p.b){
      sp.0 = MDSP_gamma(B.h.new[k,], max.iter, 0.001, centers[k])
      g.h[k] = sp.0$gamma - 0.01
      g.h.new[k] = g.h[k] + 0.01
      cl.h[k,] = sp.0$l.group
    }

    iter.nu=1
    ### find gamma, nu
    for (k in 1:p.b){
      b.t.new=B.h.new[k,]+kpp^{-1}*Lm[k,]
      while (iter.nu < max.iter & abs(g.h[k]-g.h.new[k])>0.001){
        g.h[k]=g.h.new[k]
        #lm.s.ad=lm.s*min(iter/5,1)
        lm.s.ad = lm.s * mean(m_vec)
        #lm.s.ad=lm.s
        nu.0=sign(b.t.new)*as.numeric(abs(b.t.new)>lm.s.ad/kpp)*(abs(b.t.new)-lm.s.ad/kpp)
        nu.g=g.h[k]+sign(b.t.new-g.h[k])*as.numeric(abs(b.t.new-g.h[k])>lm.s.ad/kpp)*(abs(b.t.new-g.h[k])-lm.s.ad/kpp)
        Nu.h.new[k,]=nu.g*cl.h[k,]+nu.0*(1-cl.h[k,])

        sp.h=MDSP_gamma(Nu.h.new[k,], max.iter, 0.001, centers[k])
        g.h.new[k] = sp.h$gamma
        cl.h[k,] = sp.h$l.group
        iter.nu = iter.nu+1
        if(iter.nu == max.iter) {print("WARNING: maxium iteration (nu) achieves")}
      }
    }

    ### lambda updates
    nu.h.new=as.vector(Nu.h.new)
    lm.h.new=lm.h+kpp*(b.h.new-nu.h.new)

    dif <- sqrt(mean((theta.h-theta.h.new)^2))
    # cat("diff=", round(dif, 4), "\n")

    iter=iter+1
    if(iter==max.iter) {print("WARNING: maxium iteration achieves")}
  }

  return(list(alpha=a.h.new, Beta=B.h.new, Nu=Nu.h.new,
              gamma.h=g.h.new, group=cl.h, it.errors=error[1:(iter-1)]))

  print(proc.time()-ptm)

}
