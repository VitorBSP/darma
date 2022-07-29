library(tibble)
library(dplyr)

#Funções:
reg.eged<-function(n,alpha=NULL,delta=NULL,lambda,sigma,eta,gama,beta,X,tau=0.5)
{
  #X matriz de covariáveis
  #tau<-quantil utilizado na Parametrização
  mut<-exp(X%*%beta)
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
# Estimação usando a Log Verossimilhança Reparametrizada pelo quantil tau
LQV.Reg.EGEDD<-function(theta,dados,X,tau,q=1,repa="alpha0")
{
  beta<-theta[6:length(theta)]
  y<-dados
  mut<-exp(X%*%beta)
  if(repa=="alpha0"){
    #Parametrizando alpha
    delta<-theta[1]
    lambda<-theta[2]
    sigma<-theta[3]
    eta<-theta[4]
    gama<-theta[5]
    A<-(1-tau)^(1/lambda)
    B<-(1-A)^(1/eta)
    C<-(1-B)^(1/gama)
    D<-(1-C)^(-1/sigma)
    alpha<-(D-1)*exp(delta*log(mut))
  }
  if(repa=="delta0"){
    #Parametrizando delta
    alpha<-theta[1]
    lambda<-theta[2]
    sigma<-theta[3]
    eta<-theta[4]
    gama<-theta[5]
    A<-(1-tau)^(1/lambda)
    B<-(1-A)^(1/eta)
    C<-(1-B)^(1/gama)
    D<-(1-C)^(-1/sigma)
    delta<-(-1/log(mut))*log((1/alpha)*(D-1))
  }
  A<-(1+alpha*(y^(-delta)))^(-sigma-1)
  B<-1-(1+alpha*(y^(-delta)))^(-sigma)
  fdp<-suppressWarnings(lambda*alpha*sigma*delta*eta*gama*
                          (y^(-delta-1))*A*(B^(gama-1))*((1-(B^gama))^(eta-1))*
                          ((1-(1-(B^gama))^eta)^(lambda-1)))
  if(q<=0 | q>1){paste("Parameter q out of range 0<q<1")}
  if(q>0 & q<1){#fdp função densidade probabilidade
    LQ<-((fdp^(1-q))-1)/(1-q)
  }
  if(q==1){
    LQ<-log(fdp)
  }
  LQV<-sum(LQ)
  return(-LQV)
}

# Definindo os valores dos parâmetros e condições da simulação
n_ = c(50,100,200)
R = 1000
alpha = NULL
delta = 3.5
lambda = 11.2
sigma = 3
eta = 5
gama = 3
beta = c(0.5, 0.2, 1.5)
tau = 0.5
theta0<-c(alpha,delta,lambda,sigma,eta,gama,beta)
k<-length(beta)-1
l = vector(mode='list', 3)
contador_falhas = 0


#Simulação utilizando OPTIM
for(n in n_){
  m = matrix(nrow=R, ncol=8)
  i = 1
  while(i<=R){
    X<-cbind(rep(1,n),matrix(runif(n*k),nrow=n,ncol=k))
    y<-as.numeric(reg.eged(n=n,alpha=alpha,delta=delta,lambda=lambda,sigma=sigma,eta=eta,
                           gama=gama,beta=beta,X=X,tau=tau))
    ynew <-log(y)
    ajuste<-try(lm(ynew~X+0),T)
    mqo <- c(ajuste$coef)
    inf<-c(rep(0.03,5),rep(min(mqo),length(mqo)))
    sup<-c(rep(10,5),rep(max(abs(mqo)),length(mqo)))
    tmp<-try(optim(theta0,fn=LQV.Reg.EGEDD,lower=inf,upper=sup,
                   method="CG",dados=y,X=X,tau=tau,q=0.95,
                   hessian=TRUE),T)
    
    if(class(tmp) == "try-error"){
      contador_falhas=contador_falhas+1
    }else{
      m[i,]= tmp$par
      i = i + 1
    } 
  }
  l[[match(n,n_)]] =  m
}

# Avaliando a simulação
resultados = data.frame(matrix(ncol = 8, nrow  = 0))
colnames(resultados) = c('mean', 'sd' ,'parametros', 'Vies', 'EQM',
                         'Assimetria', 'Curtose', 'N')
j = 1
for(df in l){
 df_ = df |> 
    colMeans() %>%
    cbind(apply(df, 2, sd)) %>%
    as.data.frame() |> 
    rename(mean = `.`, sd = V2) |>
    mutate(parametros = theta0) |>
    mutate(Vies = mean - parametros) |>
    mutate(EQM = Vies^2 + sd^2,
           Assimetria = moments::skewness(df),
           Curtose = moments::kurtosis(df),
           N = rep(n_[j], length(theta0))) 
 resultados = resultados %>%
   add_row(df_)
 j = j +1
}

# Exportando simulação:
write.csv(resultados_delta, 'simu_delta.csv')