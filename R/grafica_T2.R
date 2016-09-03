#' Grafica T2 de Hotelling
#'
#' Esta funcion permite realizar la grafica de control T2 de Hotelling y obtener sus ARL's bajo y fuera de control.
#' @param mu: El vector de las medias de las variables (Parametros del proceso cuando esta bajo control).
#' @param sigma: La matriz de varianzas y covarianzas (Parametros del proceso cuando esta bajo control).
#' @param X: es la matriz que contiene las variables
#' @param LCS: El límite de control superior.
#' @param d: La distancia de mahalanobis.
#' @param k: El tamaño de la muestra.
#' @return Grafico de control multivariante T2 de Hotelling;
#' @return ARL0: Promedio de longitud de rachas cuando el proceso esta bajo control;
#' @return ARL1: Promedio de longitud de rachas cuando el proceso esta fuera de control.
#' @keywords T2 de Hotelling
#' @references [1] Aparisi, F.; Epprecht, E.; Ruiz, O.; Veiga A., "Reducing Sampling Costs of Multivariate SPC with a Double-Dimension T2 Control Chart",” International Journal of Production Research, vol. 1, no. 1, pag. 1-15, 2013;
#' @references [2] Aparisi, F.; Epprecht, E.; Ruiz, O., "T2 Control Charts with Variable Dimension", Journal of Quality Technology, vol. 44, no. 4, pp. 375-393, 2012;
#' @references [3] Ruiz, O., “Gráficos de control de calidad multivariantes con dimension variable” Tesis doctoral, Dept. de Estadistica e Investigación Operativa, Univ. Pol. de Valencia, España, 2013.
#' @examples
#' ###Simulación de datos normales multivariados###
#'mu=c(4.5,7,8.45)
#'sigma=matrix(c(2,1.5,2.4,1.5,3,3.1,2.4,3.1,4),ncol=3)
#'#la matriz de covarianzas debe ser definida positiva
#'#comprobación por medio de los valores propios de la matriz (>0)
#'eigen(sigma)
#'library(MASS)
#'#datos multivariantes
#'dmulti=mvrnorm(250,mu,sigma)
#'colnames(dmulti)=c("va1","va2","va3")
#'#Datos a utilizar
#'#n1----->50 dmulti[1:50,]
#'#n2----->100 dmulti[1:100,]
#'#n5----->250 dmulti
#'dmulti_n1=dmulti[1:50,]
#'dmulti_n2=dmulti[1:100,]

#'#vector de medias bajo control
#'mu0=c(5.4,6.8,8.5)

#'#uso de la funcion T2 de Hotelling 
#'grafica_T2(mu0,sigma,dmulti_n1,k=1,LCS=14.321,d=1.2)#n=1
#'grafica_T2(mu0,sigma,dmulti,k=5,LCS=14.321,d=2)#n=5
grafica_T2<-function(mu,sigma,X,LCS,d,k){
  requireNamespace("calibrate","stats","graphics","MASS")
    if (k>1){  
    covarianzas<-sigma
    medias=matrix(c(rep(0,ncol(X)*I(nrow(X)/k))),ncol=ncol(X))
    for(i in 1:I(nrow(X)/k)){
      medias[i,]=colMeans(X[I(abs(2*i-k+1)):I(i*k),])
    }
    T2<-c(rep(0,nrow(medias)))
    muestras<-1:nrow(medias)
    for(i in 1:(nrow(medias))){
      T2[i]<-k*(t(medias[i,]-mu)%*%solve(covarianzas)%*%(medias[i,]-mu))
    }
    T2
    LCI<-0
    plot(muestras,T2,type='l',col="blue",xaxt='n',xlim=c(0,(length(T2)+18)),ylim=c(min(LCI,min(T2))-0.5,max(max(T2),LCS)+0.5))
    axis(1,seq(from=0,to=length(T2),by=2),cex.axis=0.8,las=2)
    
    abline(h=LCS,lty=2,col="brown")
    text((max(muestras)+5),LCS,paste('LCS=',round(LCS,digits = 6)),pos=3,font=2,cex=0.8)
    abline(h=LCI,lty=2,col="brown")
    text((max(muestras)+5),LCI,paste('LCI=',LCI),pos=3,font=2,cex=0.8)
    for(i in 1:length(muestras)){
      if((T2[i]<LCI)|(T2[i]>LCS)){points(muestras[i],T2[i],pch=4,col="red")}
      else{points(muestras[i],T2[i],pch=20,col="green")}
    }
    mtext('Grafico de control T2 de Hotelling',side=3,font=2)}
  if(k==1){
    covarianzas<-sigma
    T2<-vector(length=nrow(X))
    muestras<-1:nrow(X)
    for(i in 1:(nrow(X))){T2[i]<-t(X[i,]-mu)%*%solve(covarianzas)%*%(X[i,]-mu)}
    T2
    LCI<-0
    plot(muestras,T2,type='l',col="blue",xaxt='n',xlim=c(0,(length(T2)+10)),
         ylim=c(min(LCI,min(T2))-0.5,max(max(T2),LCS)+1))
    axis(1,seq(from=0,to=length(T2),by=2),cex.axis=0.8,las=2)
    
    abline(h=LCS,lty=2,col="brown")
    text((max(muestras)+5),LCS,paste('LCS=',round(LCS,digits = 6)),pos=3,font=2,cex=0.8)
    abline(h=LCI,lty=2,col="brown") 
    text((max(muestras)+5),LCI,paste('LCI=',LCI),pos=3,font=2,cex=0.8)
    for(i in 1:length(muestras)){        
      if((T2[i]<LCI)|(T2[i]>LCS)){ points(muestras[i],T2[i],pch=4,col="red")}
      else {points(muestras[i],T2[i],pch=20,col="green")  }
    }
    mtext('Grafico de control T2 de Hotelling',side=3,font=2)
  }
  #CALCULANDO LOS ARL
  ARL_T2=function(p,d,n,LCS){
    beta=pchisq(LCS,df=p,ncp=n*d*d,lower.tail=FALSE)
    ARL=1/beta
    return(ARL)
  }
  ARL0=ARL_T2(p=ncol(X),d=0,n=k,LCS)
  ARL1=ARL_T2(p=ncol(X),d=d,n=k,LCS)
  return(list(ARL0=ARL0,ARL1=ARL1))
}