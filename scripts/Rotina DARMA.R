
verossimilhanca_condicional_DARMA = function(x,y,p=1,q=1,mut=0.5){
  delta = x[p+q+1]
  lambda = x[p+q+2]
  sigma = x[p+q+3]
  eta = x[p+q+4]
  gama = x[p+q+5]
  beta = x[p+q+6:length(x)]
  A<-(1-tau)^(1/lambda)
  B<-(1-A)^(1/eta)
  C<-(1-B)^(1/gama)
  D<-(1-C)^(-1/sigma)
  alpha<-(D-1)*exp(delta*log(mut))
  n<-length(y)
  e<-rep(0,n) 
  for(t in (p+1):n){
    tmp<-sum(c(x[1:p])*y[t-c(1:p)],-x[(p+1):(p+q)]*e[t-c(1:q)])
    e[t]<-y[t]-tmp
  }
  y = e
  fdp<-suppressWarnings(lambda*alpha*sigma*delta*eta*gama*
                          (y^(-delta-1))*A*(B^(gama-1))*((1-(B^gama))^(eta-1))*
                          ((1-(1-(B^gama))^eta)^(lambda-1)))
  fdp = fdp[!(is.infinite(fdp) | is.na(drop(fdp)))]
  loglike<-suppressWarnings(sum(log(fdp)))
  return(loglike)
}


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
#Podemos reparametizar por alpha ou delta
reged.rep<-function(n,alpha=NULL,delta=NULL,lambda,sigma,eta,gama,beta,mut,tau=0.5)
{
  #X matriz de covariáveis
  #tau<-quantil utilizado na Parametrização
  A<-(1-tau)^(1/lambda)
  B<-(1-A)^(1/eta)
  C<-(1-B)^(1/gama)
  D<-(1-C)^(-1/sigma)
  if(is.null(alpha)){
    #Parametrizando alpha
    alpha<-(D-1)*exp(delta*log(mut))
  }
  if(is.null(delta)){
    #Parametrizando delta
    delta<-(-1/log(mut))*log((1/alpha)*(D-1))
  }
  xt<-reged(n=n,alpha=alpha,lambda=lambda,sigma=sigma,delta=delta,eta=eta,gama=gama)
  return(xt)
}
#-------------------------------------------------------------------------------
#Geração de uma a.a. EGED-ARMA
#-------------------------------------------------------------------------------
egedarma.sim<-function(model=list(phi=NA,theta=NA),n,const,alpha=NULL,delta=NULL,lambda,sigma,eta,
                       gama,beta,X,tau=0.5,freq,link="log")
{
  #model=list(ar=c(phi_1,...,phi_p),ma=c(theta_1,...,theta_q))
  if (!is.list(model)) 
    stop("'model' must be list")
  if (n <= 0L) 
    stop("'n' must be strictly positive")
  
  p<-NA
  q<-NA
  if(any(is.na(model$phi)==F)){p<-1:length(model$phi)}
  if(length(p)>1){
    minroots<-min(Mod(polyroot(c(1, -model$phi))))
    if(minroots<=1) 
      stop("'ar' part of model is not stationary")
  }
  if(any(is.na(model$theta)==F)){q<-1:length(model$theta)}
  linktemp<-substitute(link)
  if(!is.character(linktemp))
  {
    linktemp<-deparse(linktemp)
    if (linktemp=="link")
      linktemp<-eval(link)
  }
  if(any(linktemp == c("log", "identity")))
  {  
    stats<-make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"log\"  and \"identity\""))
  } 
  link<-structure(list(link=linktemp,
                       linkfun=stats$linkfun,
                       linkinv=stats$linkinv))
  linkfun<-link$linkfun
  linkinv<-link$linkinv
  #-------BEGIN ARMA MODEL------------------------------------------------------
  if(any(is.na(phi)==F) && any(is.na(theta)==F))
  {
    #print("ARMA model")
    p<-max(phi)
    q<-max(theta)
    m<-50 #burn-in
    maxx<-max(p,q)
    ynew<-rep(const,(n+m))
    mu<-linkinv(ynew)
    error<-rep(0,n+m) # E(error)=0 
    etat<-y<-NULL
    #X<-cbind(sin(2*pi*(1:(n+m))/12))
    for(i in (maxx+1):(n+m))
    {
      etat[i]<-const+X[i,]%*%as.matrix(beta)+(phi%*%(ynew[i-phi]-X[i-phi,]%*%as.matrix(beta)))+(theta%*%error[i-theta]) 
      mu[i]<-linkinv(etat[i])
      y[i]<-reged.rep(n=1,alpha=alpha,lambda=lambda,sigma=sigma,delta=delta,eta=eta,gama=gama,mut=mu[i],tau=tau)
      ynew[i]<-linkfun(y[i]) ## Laís incluiu
      error[i]<-ynew[i]-etat[i]   
    }
    return(ts(y[(m+1):(n+m)],frequency=freq))
  }
}
#-------------------------------------------------------------------------------
alpha<-NULL
delta<-1.5
sigma<-2.0
gama<-3.0
eta<-5.0
lambda<-2.0
n<-1000
beta<-c(1,0.5,2.3)
tau<-0.5
freq<-1
const<-2
phi<-c(0.5)
theta<-c(0.8)
k<-length(beta)-1
burn.in<-50
X<-cbind(rep(1,n+burn.in),matrix(runif((n+burn.in)*k),nrow=n+burn.in,ncol=k))
mut<-exp(X%*%beta)

y<-reged.rep(n=1050,alpha=alpha,lambda=lambda,sigma=sigma,delta=delta,eta=eta,gama=gama,mut=mut,tau=tau)
plot.ts(y)



ts.sim<-SIMU.EGEDARMA(n=n,const=const,phi=phi,theta=theta,alpha=alpha,lambda=lambda,
                      sigma=sigma,delta=delta,eta=eta,gama=gama,beta=beta,X=X,tau=tau,freq=freq)

ts.sim <- egedarma.sim(list(phi, theta), n = n, const = const, alpha = alpha,
                       delta = delta, lambda = lambda, sigma = sigma, eta  = eta, 
                       gama = gama, beta = beta, X = X, freq = 1)

n=101
alpha<-NULL
delta<-3.5
sigma<-11.2
gama<-3.0
eta<-5
lambda<-3
n<-1000
beta<-c(0.5,0.2,1.5)
tau<-0.5
freq<-1
const<-2
phi<-c(0.5)
theta<-c(0.8)
k<-length(beta)-1
burn.in<-50
X<-cbind(rep(1,n+burn.in),matrix(runif((n+burn.in)*k),nrow=n+burn.in,ncol=k))
mut<-exp(X%*%beta)
R=10
contador_nao_convergencias=contador_falhas=0
estimativa_parametros_condicional=matrix(rep(NA, R*3), ncol=3)
i=1
while(i<=R)
{
  dado = egedarma.sim(list(phi, theta), n = n, const = const, alpha = alpha,
                      delta = delta, lambda = lambda, sigma = sigma, eta  = eta, 
                      gama = gama, beta = beta, X = X, freq = 1) #ação numeros aleatorios
  dado = dado[-length(dado)]
  theta0 <- runif(8,0,1) # Chute inicial
  
  fit=try(optim(par = theta0,fn=verossimilhanca_condicional_DARMA, y=dado,method = "CG",hessian = T,p=1,q=1,
                control = list(fnscale=-1)))
  
  if(class(fit) == "try-error"){contador_falhas=contador_falhas+1} 
  if(class(fit) != "try-error")
  {
    
    if(fit$conv!=0){contador_nao_convergencias=contador_nao_convergencias+1}
    estimativa_parametros_condicional[i,]=c(fit$par)
    
    
    if((100*i/R)%%10 == 0)
      print(c(100*i/R, "%"), quote = F)
    i=i+1
    
  }
}


alpha<-NULL
delta<-1.5
sigma<-2.0
gama<-3.0
eta<-5.0
lambda<-2.0
n<-1000
beta<-c(1,0.5,2.3)
tau<-0.5
freq<-1
const<-2
phi<-c(0.5)
theta<-c(0.8)
k<-length(beta)-1
burn.in<-50
X<-cbind(rep(1,n+burn.in),matrix(runif((n+burn.in)*k),nrow=n+burn.in,ncol=k))
mut<-exp(X%*%beta)

y<-reged.rep(n=1050,alpha=alpha,lambda=lambda,sigma=sigma,delta=delta,eta=eta,gama=gama,mut=mut,tau=tau)
plot.ts(y)

ts.sim <- egedarma.sim(list(phi, theta), n = n, const = const, alpha = alpha,
                       delta = delta, lambda = lambda, sigma = sigma, eta  = eta, 
                       gama = gama, beta = beta, X = X, freq = 1)

nsp=5
mint=20
inf<-rep(0.01,10)
sup<-rep(5,10)
y = ts.sim[-length(ts.sim)]
x = c(phi,theta, delta, sigma, gama, eta,lambda, beta)
out<-gosolnp(fun=verossimilhanca_condicional_DARMA,n.restarts=nsp,LB=inf,UB=sup,n.sim=mint,x=x)
optim(par =x, fn =verossimilhanca_condicional_DARMA, method = 'CG', y = y,
      control = list(fnscale=-1))
