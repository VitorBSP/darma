#------Pacotes Necessários------------------------------------------------------
library(bigstatsr)
library(plyr)
library(parallel)
library(doParallel)
library(foreach)
#------Carrega as Funções Necessárias-------------------------------------------
source("SIM.DARMA.R")
#-------------------------------------------------------------------------------

cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
#stopCluster(cl)

parallel.sim<-function(n){
  re<-100
  delta<-3.5
  sigma<-11.2
  gama<-3.0
  eta<-5
  lambda<-3
  beta<-c(0.5,0.2,1.5)
  tau<-0.5
  freq<-1
  const<-0
  phi<-c(0.5)
  theta<-c(0.8)
  theta0<-c(phi,theta,delta,sigma,gama,eta,lambda)
  burn.in<-100
  k<-length(beta)-1
  #Com Covariáveis
  #X<-cbind(rep(1,n+burn.in),matrix(runif((n+burn.in)*k),nrow=n+burn.in,ncol=k))
  #Sem Covariáveis
  X<-cbind(rep(0,n+burn.in),matrix(rep(0,(n+burn.in)*k),nrow=n+burn.in,ncol=k))
  amostra<-matrix(ncol=re,nrow=n)
  for(j in 1:re){
    amostra[,j]<-darma.sim(n=n,kapa=0,phi=phi,theta=theta,delta=delta,lambda=lambda,sigma=sigma,
                           eta=eta,gama=gama,beta=beta,X=X,burn.in=burn.in,tau=tau,freq=1,link="log")
  }
  estimate <- FBM(re,7)
  time<-system.time(
    foreach(j=1:re, .packages=c("foreach")) %dopar%{ 
      source("SIM.DARMA.R")
      estimate[j,]<-SIM.DARMA(start.theta=theta0,y=amostra[,j],method=c("BFGS"))
    }
  )
  output<-as.data.frame(estimate[])
  save(output,file=paste0("simun",n,"CG",".RData"))
  #save(output,file=paste0("simun",n,"CG",".RData"))
  return(time)
}

parallel.sim(100)




parallel.sim<-function(n){
  re<-5100
  delta<-3.5
  sigma<-11.2
  gama<-3.0
  eta<-5
  lambda<-3
  beta<-c(0.5,0.2,1.5)
  tau<-0.5
  freq<-1
  const<-0
  phi<-c(0.5)
  theta<-c(0.8)
  theta0<-c(phi,theta,delta,sigma,gama,eta,lambda)
  burn.in<-100
  k<-length(beta)-1
  #Com Covariáveis
  #X<-cbind(rep(1,n+burn.in),matrix(runif((n+burn.in)*k),nrow=n+burn.in,ncol=k))
  #Sem Covariáveis
  X<-cbind(rep(0,n+burn.in),matrix(rep(0,(n+burn.in)*k),nrow=n+burn.in,ncol=k))
  amostra<-matrix(ncol=re,nrow=n)
  for(j in 1:re){
    amostra[,j]<-darma.sim(n=n,kapa=0,phi=phi,theta=theta,delta=delta,lambda=lambda,sigma=sigma,
                           eta=eta,gama=gama,beta=beta,X=X,burn.in=burn.in,tau=tau,freq=1,link="log")
  }
  #estimate <- FBM(re,7)
  estimate <- matrix(nrow=re,ncol=7)
  time<-system.time(
    for(j in 1:re){ 
      estimate[j,]<-SIM.DARMA(start.theta=theta0,y=amostra[,j],method=c("BFGS"))
    }
  )
  #output<-as.data.frame(estimate[])
  output<-as.data.frame(estimate)
  save(output,file=paste0("simun",n,"CG",".RData"))
  #save(output,file=paste0("simun",n,"CG",".RData"))
  return(time)
}






parallel.sim(100)


parallel.sim(60)
