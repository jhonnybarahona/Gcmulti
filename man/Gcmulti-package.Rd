\name{Gcmulti-package}
\alias{Gcmulti-package}
\alias{Gcmulti}
\docType{package}
\title{
Graficos de Control Multivariante
}
\description{
Permite realizar los graficos de control multivariante T2 de Hotelling, de doble dimension y de dimension variable. Ademas calcula los ARL  bajo y fuera de control.
}
\details{
\tabular{ll}{
Package: \tab Gcmulti\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-08-29\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Jhonny Barahona\cr
Omar Ruiz\cr
Maintainer: Jhonny Barahona <jfbaraho@espol.edu.ec>
}
\references{
[1] Aparisi, F.; Epprecht, E.; Ruiz, O.; Veiga A., "Reducing Sampling Costs of Multivariate SPC with a Double-Dimension T^2 Control Chart",” International Journal of Production Research, vol. 1, no. 1, pag. 1-15, 2013.\cr
[2] Aparisi, F.; Epprecht, E.; Ruiz, O., "T2 Control Charts with Variable Dimension", Journal of Quality Technology, vol. 44, no. 4, pp. 375-393, 2012.\cr
[3] Ruiz, O., “Gráficos de control de calidad multivariantes con dimensión variable” Tesis doctoral, Dept. de Estadística e Investigación Operativa, Univ. Pol. de Valencia, España, 2013.\cr
}
\examples{
###Simulaci?n de datos normales multivariados###
mu=c(4.5,7,8.45)
sigma=matrix(c(2,1.5,2.4,1.5,3,3.1,2.4,3.1,4),ncol=3)
#la matriz de covarianzas debe ser definida positiva
#comprobaci?n por medio de los valores propios de la matriz (>0)
eigen(sigma)
library(MASS)
#datos multivariantes
dmulti=mvrnorm(250,mu,sigma)
colnames(dmulti)=c("va1","va2","va3")
View(dmulti)
#Datos a utilizar
#n1----->50 dmulti[1:50,]
#n2----->100 dmulti[1:100,]
#n5----->250 dmulti
dmulti_n1=dmulti[1:50,]
dmulti_n2=dmulti[1:100,]

#vector de medias bajo control
mu0=c(5.4,6.8,8.5)

#uso de la funci?n T2 de Hotelling 
#para p1
grafica_T2(mu0[1:2],sigma[1:2,1:2],dmulti_n1[,1:2],k=1,LCS=11.983,d=1) #n=1
#para p
grafica_T2(mu0,sigma,dmulti_n1,k=1,LCS=14.321,d=1.2)#n=1

#uso de la funci?n T2 con dimensi?n variable
grafica_T2DV(p1=2,p=3,w=3.83,CLp1=17.09,CLp=10.62,disp1=1,disp=1.2,mu=mu0,sigma=sigma,X=dmulti_n1,k=1)

#Uso de la funci?n T2 con doble dimensi?n
grafica_T2DD(p1=2,p=3,w=2.89,CLP1=14.07,CLP=14.11,disp1=1,disp=1.2,mu=mu0,sigma=sigma,X=dmulti_n1,k=1)
}
