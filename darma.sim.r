#Exponentiated generalized exponential Dagum distribution
#-------------------------------------------------------------------------------
#Gera Números Aleatórios
reged<-function(n,alpha,lambda,sigma,delta,eta,gama)
{
  u<-runif(n)
  A<-(1-u)^(1/lambda)
  B<-(1-A)^(1/eta)
  C<-(1-B)^(1/gama)
  D<-(1-C)^(-1/sigma)
  result<-((1/alpha)*(D-1))^(-1/delta)
  return(result)
}
#-------------------------------------------------------------------------------
#theta<-c(alpha,delta,lambda,sigma,delta,eta,gama)
#Podemos reparametizar por alpha
reged.rep<-function(n,delta,lambda,sigma,eta,gama,mut,tau)
{
  #mut são os quantís condicionais
  #Parametrizando alpha
  A<-(1-tau)^(1/lambda)
  B<-(1-A)^(1/eta)
  C<-(1-B)^(1/gama)
  D<-(1-C)^(-1/sigma)
  alpha<-(D-1)*exp(delta*log(mut))
  xt<-reged(n=n,alpha=alpha,lambda=lambda,sigma=sigma,delta=delta,eta=eta,gama=gama)
  return(xt)
}
#-------------------------------------------------------------------------------
#Geração de uma a.a. EGED-ARMA
#-------------------------------------------------------------------------------
darma.sim <- function(n,kapa=0,phi=NA,theta=NA,delta,lambda,sigma,eta,gama,beta,
                      X,burn.in=50,tau=0.5,freq=1,link="log")
{
  #kapa<-constante
  m<-burn.in
  ar <- NA
  ma <- NA
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("log", "identity")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"log\"  and \"identity\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv
  )
  )
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  #X <- cbind(sin(2*pi*(1:(n+m))/50))
  ###### ARMA model
  if(any(is.na(phi)==F) && any(is.na(theta)==F))
  {
    #print("ARMA model")
    
    p <- max(ar)
    q <- max(ma)
    #m <- 50
    maxx=max(p,q)
    
    ynew <-rep(kapa,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) # E(error)=0 
    etat<- y <- NULL
    for(i in (maxx+1):(n+m))
    {
      etat[i] <- kapa + X[i,]%*%as.matrix(beta)+(phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta)))+(theta%*%error[i-ma]) 
      mu[i] <- linkinv(etat[i])
      y[i]<-reged.rep(n=1,delta=delta,lambda=lambda,sigma=sigma,eta=eta,gama=gama,mut=mu[i],tau=tau)
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-etat[i]   
    }
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # ARMA model
  
  # AR model
  if(any(is.na(phi)==F) && any(is.na(theta)==T))
  {
    #print("AR model")    
    
    p <- max(ar)
    #m <- 50
    
    ynew <-rep(kapa,(n+m))
    mu <- linkinv(ynew)
    
    etat <- y <- NULL
   
    
    etat <- kapa + X[(p+1),]%*%as.matrix(beta) + (phi%*%(ynew[(p+1)-ar]-X[(p+1)-ar,]%*%as.matrix(beta) ))
    muu <- linkinv(etat)
    
    for(i in (p+1):(n+m))
    {
      etat[i]<-kapa + X[i,]%*%as.matrix(beta)+(phi%*%(ynew[i-ar]-X[i-ar,]%*%as.matrix(beta)))
      mu[i]<-linkinv(etat[i])
      y[i]<-y[i]<-reged.rep(n=1,delta=delta,lambda=lambda,sigma=sigma,eta=eta,gama=gama,mut=mu[i],tau=tau)
      ynew[i]<-linkfun(y[i])
      
    }
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # AR model
  
  # MA model
  if(any(is.na(phi)==T) && any(is.na(theta)==F))
  {
    #print("MA model")    
    
    q <- max(ma)
    #m <- 50
    
    ynew <-rep(kapa,(n+m))
    mu <- linkinv(ynew)
    
    etat <- y <- error <- rep(0,n+m) 

    for(i in (q+1):(n+m))
    {
      #### Lais: note que está subtraindo o theta%*%error, então na verossimilhança
      # precisa ser assim tbm. Estava somando.
      etat[i] <- kapa + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma]) 
      mu[i] <- linkinv(etat[i])
      y[i]<-reged.rep(n=1,delta=delta,lambda=lambda,sigma=sigma,eta=eta,gama=gama,mut=mu[i],tau=tau)
      ynew[i] <- linkfun(y[i])  ## Laís incluiu
      error[i]<- ynew[i]-etat[i]   
    }
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  } # fim MA model  
  
}
