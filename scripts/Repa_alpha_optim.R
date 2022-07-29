library(tibble)
library(ggplot2)
library(dplyr)
library(Rsolnp)

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
LQV.Reg.EGEDD<-function(theta,dados,X,tau,q=1,repa="delta0")
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
]
# Definindo os valores dos parâmetros e condições da simulação
# ATENÇÃO AO TAU
n_ = c(50,200,500)
R = 1000
alpha = 3.5
delta = NULL
lambda = 11.2
sigma = 3
eta = 5
gama = 3
beta = c(0.5, 0.2, 1.5)
tau = 0.9 # Devemos varia o tau
theta0<-c(alpha,delta,lambda,sigma,eta,gama,beta)
k<-length(beta)-1
replicas_alpha = vector(mode='list', 3)
contador_falhas = 0
X<-cbind(rep(1,n),matrix(runif(n*k),nrow=n,ncol=k))

#Simulação utilizando OPTIM
for(n in n_){
  estimativa_parametros = matrix(nrow=R, ncol=8)
  i = 1
  while(i<=R){
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
      estimativa_parametros[i,]= tmp$par
      i = i + 1
    } 
  }
  replicas_alpha[[match(n,n_)]] =  estimativa_parametros
}

# Avaliando a simulação
resultados_alpha = data.frame(matrix(ncol = 9, nrow  = 0))
colnames(resultados_alpha) = c('mean', 'sd' ,'parametros', 'Vies', 'EQM',
                         'Assimetria', 'Curtose', 'N', 'tau')
j = 1
for(df in replicas_alpha){
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
           N = rep(n_[j], length(theta0)),
           tau = rep(tau, length(theta0))) 
  resultados_alpha = resultados_alpha %>%
    add_row(df_)
  j = j +1
}

resultados_alpha = resultados_alpha %>% rename(Média = mean, Viés = Vies, 
                                               DP = sd, EQM = EQM, 
                                               CA = Assimetria, CC = Curtose, 
                                               N = N)


#Construindo gráfico
resultados_alpha %>% filter(parametros == 3.5) %>%
ggplot() +
  geom_line(aes(x=N, y=Viés, group= as.factor(tau), colour = as.factor(tau)), 
            size = 2) +
  geom_point(aes(x=N, y=Viés, colour = as.factor(tau)), size = 3) + 
  theme_minimal() +
  labs(x = 'Tamanho amostral', y= "Viés", colour = "Quantis")

#Construindo latex
resultados_alpha %>% filter(N == 50, tau = 0.2) %>% 
  select(-parametros, -N) %>% 
  mypdf1::pdf1_tbl(format='latex')

#Exportando simulação
resultados_alpha = read.csv('simu.csv')
