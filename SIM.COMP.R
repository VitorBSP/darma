SIM.COMP<-function(coef,par,names.pars=NULL,arquivo=NULL)
{
  library(knitr)
  library(xtable)
  #-----------------------------------------------------------------------------
  #nsign<-nível de significância do Teste de Wald
  theta0<-par
  re<-length(coef[,1])
  coeff<-coef
  #-------Calculation of accuracy measures--------------------------------------
  vmean<-colMeans(coeff)
  vbias<-vmean-theta0
  vsd<-apply(coeff,2,sd)
  vmse<-vbias^2+vsd^2
  vac<-moments::skewness(coeff)
  vk<-moments::kurtosis(coeff)
  #-----------IC e TH-----------------------------------------------------------
  results<-matrix(0,nrow=6,ncol=length(theta0))
  results<-round(rbind(vmean,vbias,vsd,vmse,vac,vk),4)
  if(!is.null(names.pars)){colnames(results)<-names.pars}
  P<-c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11")
  if(is.null(names.pars)){colnames(results)<-P[1:ne]}
  rownames(results)<-c("Mean","Bias","SD","MSE","AC","K")
  if(!is.null(arquivo)){
    #salva a Tabela de Resultados no arquivo------------------------------------
    write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    print(xtable(results,type="latex",digits = 4),file=arquivo,append=TRUE)
    write(c("#---------------------------"),file=arquivo,append=TRUE,ncolumns=1)
    print(xtable(t(results),type="latex",digits = 4),file=arquivo,append=TRUE)
    #---------------------------------------------------------------------------
  }else{
    return(kable(results))
  }
}

which(is.na(output))
SIM.COMP(coef=output,par=theta0,names.pars=c("phi1","phi2","theta","delta","sigma","gama","eta","lambda"),arquivo="teste.tex")
SIM.COMP(coef=output,par=theta0,names.pars=c("phi1","phi2","theta","delta","sigma","gama","eta","lambda"))
