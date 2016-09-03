#' Grafico de doble dimension DDT2
#'
#' Esta funcion permite realizar la grafica de control de doble dimension y obtener sus ARL's bajo y fuera de control.
#' @param p1: Numero de variables p1 (Variables faciles de medir o baratas).
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
#' @param k: El tamano de la muestra.
#' @return Grafico de control multivariante de doble dimension.
#' @return  ARL0: Promedio de longitud de rachas cuando el proceso esta bajo control.
#' @return  ARL1: Promedio de longitud de rachas cuando el proceso esta fuera de control.
#' @keywords grafico de doble dimension
#' @references [1] Aparisi, F.; Epprecht, E.; Ruiz, O.; Veiga A., "Reducing Sampling Costs of Multivariate SPC with a Double-Dimension T2 Control Chart",” International Journal of Production Research, vol. 1, no. 1, pag. 1-15, 2013;
#' @references [2] Aparisi, F.; Epprecht, E.; Ruiz, O., "T2 Control Charts with Variable Dimension", Journal of Quality Technology, vol. 44, no. 4, pp. 375-393, 2012;
#' @references [3] Ruiz, O., “Graficos de control de calidad multivariantes con dimension variable” Tesis doctoral, Dept. de Estadistica e Investigacion Operativa, Univ. Pol. de Valencia, Espana, 2013.
#' @examples
#' ###Simulacion de datos normales multivariados###
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

#'#Uso de la función T2 con doble dimensión
#'grafica_T2DD(p1=2,p=3,w=2.89,CLP1=14.07,CLP=14.11,disp1=1,disp=1.2,mu=mu0,sigma=sigma,X=dmulti_n1,k=1)
#'grafica_T2DD(p1=2,p=3,w=2.6,CLP1=13.89,CLP=14.16,disp1=1.5,disp=2,mu=mu0,sigma=sigma,X=dmulti,k=5)
#'#La muestra X debe estar ordenada [p1,p] primero las variables p1: faciles de medir ; y luego las p2: complejas de medir

grafica_T2DD<-function(p1,p,w,CLP1,CLP,disp1,disp,mu,sigma,X,k,inter=0.05)
{
  requireNamespace("calibrate","stats","graphics","MASS")
  
  ###########Para las p1 variables############################
  if(k>1){
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
    t2=NULL;lb=NULL
    cuenta=0;po=NULL
    for (i in 1:length(t2p1)){
      if(t2p1[i]>w){
        cuenta=cuenta+1
        po=c(po,i)}}
    
    tem=NULL
    for ( i in po){
      tem=c(tem,t2p[i])
    }
    v1t2=c(1:length(t2p1),po)
    v2t2=c(t2p1,tem)
    mt2=matrix(c(v1t2,v2t2),ncol=2)
    mt2=mt2[order(mt2[,1]),]
    muestras=mt2[,1]
    T2=mt2[,2]
    
    lb[1]="Tp1"
    for(i in 1:I(length(muestras)-1)){
      if(muestras[i]==muestras[i+1]){
        lb[i+1]="Tp"}
      else{lb[i+1]="Tp1"}      
    }
    #################### Graficando
    LCS<-CLP
    LCI<-CLP1
    W<-w
    plot(muestras,T2,type='l',col="blue",xaxt='n',xlim=c(0,I(nrow(X)/k)+7),ylim=c(min(W,min(T2))-0.5,max(max(T2),LCS)+0.5))
    axis(1,seq(from=0,to=length(T2),by=2),cex.axis=0.8,las=2)
    abline(h=LCS,lty=2,col="brown")
    text((max(muestras)+5),LCS,paste('LCP=',LCS),pos=3,font=2,cex=0.6)
    abline(h=LCI,lty=2,col="brown")
    text(1,LCI,paste('LCP1=',LCI),pos=3,font=2,cex=0.6)
    abline(h=W,lty=2,col="brown")
    text((max(muestras)+5),W,paste('W=',W),pos=3,font=2,cex=0.6)
    for(i in 1:length(lb)){
      if(lb[i]=="Tp1"){ points(muestras[i],T2[i],pch=20,col="green")}
      else{points(muestras[i],T2[i],pch=20,col="yellow")}
    }
    for(i in 1:length(muestras)){        
      if(T2[i]>=LCS){points(muestras[i],T2[i],pch=4,col="red")}
    }
    
    #library calibrate
    library(calibrate)
    textxy(muestras,T2,lb)
    mtext('Gráfico de control T2 con Doble Dimensión',side=3,font=2)}
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
    #################T2 para graficar################################33
    t2p1=T2_p1variables
    t2p=T2
    t2=NULL;lb=NULL
    cuenta=0;po=NULL
    for (i in 1:length(t2p1)){
      if(t2p1[i]>w){
        cuenta=cuenta+1
        po=c(po,i)}}
    
    tem=NULL
    for ( i in po){
      tem=c(tem,t2p[i])
    }
    v1t2=c(1:length(t2p1),po)
    v2t2=c(t2p1,tem)
    mt2=matrix(c(v1t2,v2t2),ncol=2)
    mt2=mt2[order(mt2[,1]),]
    muestras=mt2[,1]
    T2=mt2[,2]
    
    lb[1]="Tp1"
    for(i in 1:I(length(muestras)-1)){
      if(muestras[i]==muestras[i+1]){
        lb[i+1]="Tp"}
      else{lb[i+1]="Tp1"}      
    }
    ##########Graficando k=1
    LCS<-CLP
    LCI<-CLP1
    W<-w
    plot(muestras,T2,type='l',col="blue",xaxt='n',xlim=c(0,I(nrow(X)/k)+7),ylim=c(min(LCI,min(T2))-0.5,max(max(T2),LCS)+5))
    axis(1,seq(from=0,to=length(T2),by=2),cex.axis=0.8,las=2)
    abline(h=LCS,lty=2,col="brown")
    text((max(muestras)+5),LCS,paste('LCP=',LCS),pos=3,font=2,cex=0.6)
    abline(h=LCI,lty=2,col="brown") 
    text(1,LCI,paste('LCP1=',LCI),pos=3,font=2,cex=0.6)
    abline(h=W,lty=2,col="brown")
    text((max(muestras)+5),W,paste('W=',W),pos=3,font=2,cex=0.6)
    for(i in 1:length(lb)){
      if(lb[i]=="Tp1"){ points(muestras[i],T2[i],pch=20,col="green")}
      else{points(muestras[i],T2[i],pch=20,col="yellow")}
    }
    for(i in 1:length(muestras)){        
      if(T2[i]>=LCS){points(muestras[i],T2[i],pch=4,col="red")}
    }
    #library calibrate
    library(calibrate)
    textxy(muestras,T2,lb)
    mtext('Gráfico de control T2 con Doble Dimensión',side=3,font=2)}
  
  ################## Calculando ARL
  #(p1,p,w,CLP1,CLP,disp1,disp,mu,sigma,X,k,inter=0.05)
  ARLDDT2=function(p1, p, n, iw, iCLP1, iCLP, ides1, ides2){
    #----------ARL DDT-----
    INT1=function(v, iCLP,lamb1,lamb2,p,p1){
      arriba = (iCLP-v)
      if(arriba<0){
        resul=0}
      else{
        if(lamb2==0){
          resul=pchisq(arriba,p-p1)}
        else{
          if(I(p-p1)==1 |p1==1){
            resul=pnorm(sqrt(arriba)-sqrt(lamb2-lamb1))- pnorm(-sqrt(arriba)-sqrt(lamb2-lamb1))}
          else{dchisq(arriba, p-p1, lamb2-lamb1) - 0}
        }
      }
      return(resul)
    }
    
    #--------------------
    interv= function(iCLP1, iW){
      if(lamb2==0){
        inter=0.0001}
      else{
        if (iCLP1<5.5) {inter=0.2}
        if ((iCLP1>=5.5) & (iCLP1<11)){inter= 0.0129*iCLP1*iCLP1 - 0.2378*iCLP1 + 1.1241}
        if ((iCLP1>=11) & (iCLP1<18)) {inter= -0.0009*iCLP1*iCLP1*iCLP1 + 0.0393*iCLP1*iCLP1 - 0.5987*iCLP1 + 3.0314}
        if (iCLP1>=18){ inter=0.005}
        if (p1==1){
          if ((ides1<1) & (ides2<1)){ inter=0.0001}}  
        else{ inter= 0.0005}
      }  
      if (inter<=0){ inter=0.001}
      res2=inter
      return(res2)
    }
    
    #-------------------------
    lamb1 = n*ides1*ides1; lamb2 = n*ides2*ides2;
    pnosignal=0; part1=0
    if (p1==1){
      part1=pnorm(sqrt(iw)-sqrt(lamb1))-pnorm(-sqrt(iw)-sqrt(lamb1))}
    else{
      if (lamb2>0){ part1 = pchisq(iw, p1, lamb1)}
      else{part1= pchisq(iw,p1)}
      inter=interv(iCLP,iw)}
    xi=iw
    repeat{
      if (lamb2>0) {
        alto = INT1(xi, iCLP, lamb1,lamb2,p,p1)* dchisq(xi, p1, lamb1)}
      else{
        alto = INT1(xi, iCLP, lamb1,lamb2,p,p1)* dchisq(xi,p1)}
      pnosignal = pnosignal + alto*inter
      xi = xi + inter
      if(xi>iCLP1)break} #(iCLP1 - xi) < inter
    
    pnosignal = pnosignal + part1
    ARL = 1 / (1- pnosignal)
    return(ARL)
  }
  ARL_DDT2=ARLDDT2(p1=p1,p=p,n=k,iw=w,iCLP1=CLP1,iCLP=CLP,ides1=disp1,ides2=disp)
  return(ARL_DDT2)
}
