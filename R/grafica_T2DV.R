#' Grafico de dimension variable VDT2
#'
#' Esta funcion permite realizar graficos de control de dimension variable y obtener sus ARL's bajo y fuera de control.
#' @param p1: Numero de variables p1 (Variables fáciles de medir o baratas).
#' @param p2: Numero de variables p2 (Variables complejas de medir o costosas).
#' @param p: Numero total de variables.
#' @param w: Limite de alerta o de advertencia.
#' @param CLP1: Limite de control p1.
#' @param CLP: Limite de control p.
#' @param disp1: Distancia de mahalanobis de p1.
#' @param disp: Distancia de mahalanobis de p.
#' @param mu: El vector de las medias de las variables (Parametros del proceso cuando esta bajo control).
#' @param sigma: La matriz de varianzas y covarianzas (Parametros del proceso cuando esta bajo control).
#' @param X: es la matriz que contiene las variables.
#' @param k: El tamaño de la muestra.
#' @return Grafico de control multivariante de dimension variable.
#' @return ARL0: Promedio de longitud de rachas cuando el proceso esta bajo control.
#' @return ARL1: Promedio de longitud de rachas cuando el proceso esta fuera de control.
#' @keywords grafico de dimension variable
#' @references [1] Aparisi, F.; Epprecht, E.; Ruiz, O.; Veiga A., "Reducing Sampling Costs of Multivariate SPC with a Double-Dimension T2 Control Chart",” International Journal of Production Research, vol. 1, no. 1, pag. 1-15, 2013.
#' @references [2] Aparisi, F.; Epprecht, E.; Ruiz, O., "T2 Control Charts with Variable Dimension", Journal of Quality Technology, vol. 44, no. 4, pp. 375-393, 2012.
#' @references [3] Ruiz, O., “Gráficos de control de calidad multivariantes con dimensión variable” Tesis doctoral, Dept. de Estadística e Investigacion Operativa, Univ. Pol. de Valencia, Espana, 2013.
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

#'uso de la función T2 con dimensión variable
#'grafica_T2DV(p1=2,p=3,w=3.83,CLp1=17.09,CLp=10.62,disp1=1,disp=1.2,mu=mu0,sigma=sigma,X=dmulti_n1,k=1)
#'grafica_T2DV(p1=2,p=3,w=3.59,CLp1=18.94,CLp=10.82,disp1=1.5,disp=2,mu=mu0,sigma=sigma,X=dmulti,k=5)
#' #La muestra X debe estar ordenada [p1,p] primero las variables p1: faciles de medir ; y luego las p2: complejas de medir
grafica_T2DV<-function(p1,p,w,CLp1,CLp,disp1,disp,mu,sigma,X,k,ProbInic_p1=1){
  requireNamespace("calibrate","stats","graphics","MASS")
  if(k>1){
    ###########Para las p1 variables############################
    muestra_p1variables=X[,1:p1]
    mu_p1variables=mu[1:p1]
    sigma_p1variables=sigma[1:p1,1:p1]
    medias_p1variables=matrix(c(rep(0,p1*I(nrow(X)/k))),ncol=p1)
    for(i in 1:I(nrow(muestra_p1variables)/k)){
      medias_p1variables[i,]=colMeans(muestra_p1variables[I(abs(2*i-k+1)):I(i*k),])
    }
    T2_p1variables<-c(rep(0,nrow(medias_p1variables)))
    muestras_p1variables<-1:nrow(medias_p1variables)
    for(i in 1:(nrow(medias_p1variables))){
      T2_p1variables[i]<-k*(t(medias_p1variables[i,]-mu_p1variables)%*%solve(sigma_p1variables)%*%(medias_p1variables[i,]-mu_p1variables))
    }
    #####################################para las p variables
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
    #################T2 para graficar################################33
    t2p1=T2_p1variables
    t2p=T2
    t2=rep(0,length(t2p));lb=NULL
    i=1
    t2[1]=t2p1[1];lb[1]="Tp1"
    repeat{
      if(t2[i]>=w){
        t2[i+1]=t2p[i+1]
        lb[i+1]="Tp"}
      else{t2[i+1]=t2p1[i+1]
           lb[i+1]="Tp1"}
      i=i+1
      if(i >=length(t2))break}
    ##########################Graficando para k>1#############
    T2=t2
    LCS<-CLp
    LCI<-CLp1
    W<-w
    plot(muestras,T2,type='l',col="blue",xaxt='n',xlim=c(0,(length(T2)+7)),ylim=c(min(W,min(T2))-0.5,max(max(T2),LCS,LCI)+0.5))
    axis(1,seq(from=0,to=length(T2),by=2),cex.axis=0.8,las=2)
    abline(h=LCS,lty=2,col="brown")
    text((max(muestras)+2),LCS,paste('LCP=',LCS),pos=3,font=2,cex=0.6)
    abline(h=LCI,lty=2,col="brown")
    text((max(muestras)+6),LCI,paste('LCP1=',LCI),pos=3,font=2,cex=0.6)
    abline(h=W,lty=2,col="brown")
    text((max(muestras)+6),W,paste('W=',W),pos=3,font=2,cex=0.6)
    
    for(i in 1:length(lb)){
      if(lb[i]=="Tp1"){ points(muestras[i],T2[i],pch=20,col="green")}
      else{points(muestras[i],T2[i],pch=20,col="yellow")}
    }
    for(i in 1:length(muestras)){        
      if(T2[i]>=LCS &T2[i]>=LCI){points(muestras[i],T2[i],pch=4,col="red")}
    }
    
    #library calibrate
    library(calibrate)
    textxy(muestras,T2,lb)
    mtext('Gráfico de control T2 con Dimensión Variable',side=3,font=2)}
  
  if(k==1){
    covarianzas<-sigma
    T2<-vector(length=nrow(X))
    muestras<-1:nrow(X)
    
    ###########Para las p1 variables############################
    muestra_p1variables=X[,1:p1]
    mu_p1variables=mu[1:p1]
    sigma_p1variables=sigma[1:p1,1:p1]
    T2_p1variables<-c(rep(0,nrow(muestra_p1variables)))
    muestras_p1variables<-1:nrow(muestra_p1variables)
    for(i in 1:(nrow(muestra_p1variables))){
      T2_p1variables[i]<-k*(t(muestra_p1variables[i,]-mu_p1variables)%*%solve(sigma_p1variables)%*%(muestra_p1variables[i,]-mu_p1variables))
    }
    #####################################para las p variables
    T2<-c(rep(0,nrow(X)))
    for(i in 1:(nrow(X))){
      T2[i]<-k*(t(X[i,]-mu)%*%solve(covarianzas)%*%(X[i,]-mu))}
    #################T2 para graficar################################
    t2p1=T2_p1variables
    t2p=T2
    t2=rep(0,length(t2p));lb=NULL
    i=1
    t2[1]=t2p1[1];lb[1]="Tp1"
    repeat{
      if(t2[i]>=w){
        t2[i+1]=t2p[i+1]
        lb[i+1]="Tp"}
      else{t2[i+1]=t2p1[i+1]
           lb[i+1]="Tp1"}
      i=i+1
      if(i >=length(t2))break}
    ##########Graficando k=1
    T2=t2
    LCS<-CLp
    LCI<-CLp1
    W<-w
    plot(muestras,T2,type='l',col="blue",xaxt='n',xlim=c(0,(length(T2)+7)),ylim=c(min(LCI,min(T2))-0.5,max(max(T2),LCS,LCI)+0.5))
    axis(1,seq(from=0,to=length(T2),by=2),cex.axis=0.8,las=2)
    abline(h=LCS,lty=2,col="brown")
    text((max(muestras)+2),LCS,paste('LCP=',LCS),pos=3,font=2,cex=0.6)
    abline(h=LCI,lty=2,col="brown") 
    text((max(muestras)+6),LCI,paste('LCP1=',LCI),pos=3,font=2,cex=0.6)
    abline(h=W,lty=2,col="brown")
    text((max(muestras)+6),W,paste('W=',W),pos=3,font=2,cex=0.6)
    
    for(i in 1:length(lb)){
      if(lb[i]=="Tp1"){ points(muestras[i],T2[i],pch=20,col="green")}
      else{points(muestras[i],T2[i],pch=20,col="yellow")}
    }
    for(i in 1:length(muestras)){        
      if(T2[i]>=LCS&T2[i]>=LCI){points(muestras[i],T2[i],pch=4,col="red")}
    }
    
    #library calibrate
    library(calibrate)
    textxy(muestras,T2,lb)
    mtext('Gráfico de control T2 con Dimensión Variable',side=3,font=2)
  }
  
  ##Calculo del ARL
  n=k
  lamb1 = n*disp1*disp1; lamb2 = n*disp*disp
  a= pchisq( w, p1, lamb1)
  b= pchisq( CLp1, p1, lamb1) - a
  c= pchisq( w, p, lamb2)
  d= pchisq( CLp, p, lamb2) - c
  if(ProbInic_p1==0){
    b1=1;b2=0}
  else{b1=0;b2=1-b1}
  
  ##Calculo de Probabilidades para Qo
  a0= pchisq(w, p1, 0)
  b0= pchisq(CLp1,p1,0)-a0
  c0= pchisq(w, p, 0)
  d0= pchisq(CLp,p,0)-c0
  ##multiplicacion del vector B por la Inversa de I-Qo
  k=1/(1-d0-a0+(a0*d0)-(b0*c0))    #determinante de I-Qo
  s11=(b1*(1-d0)+b2*c)*k
  s12=(b1*b0+b2*(1-a0))*k
  k=s11+s12      #calculo ARLo utilizando la informacion anterior
  s11=s11/k   #calculo el vector S
  s12=s12/k
  #calculo ARLss(d)
  deno = a+d-(a*d)+(b*c)-1
  num = s12*(a-1) + (s11*(d-1)) - (b*s11) - (s12*c)
  ARL = num / deno
  return(ARL)
}