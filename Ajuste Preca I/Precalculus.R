################################################################################
#                                   Packages
################################################################################
#----
library(dplyr)
library(pscl)
library(writexl)
library(maxLik)
library(VGAM)
library(readr)
library(caret)
library(ggplot2)
library(scales)
################################################################################
#                                     Fit
################################################################################
#----
Zipc.maxlik<-function(X,Z,Y,status,beta){
  P<-ncol(X);s<-ncol(Z);n<-nrow(X)
  f<-function(C,mu){dpois(C,mu)} #Funcion de densidad Poisson
  I<-as.numeric(Y==0) # Indicadora de y_i=0
  I.<-1-I # Indicadora de y_i>0
  di<-status # Datos censurados, en status 0 representa censura
  C<-Y*(1-status) # Constantes de censura 
  # Vector de derivadas 
  grad<-function(beta){
    mu<-exp(X%*%beta[1:P]) # Generando mui
    w<- exp(Z%*%beta[(P+1):(P+s)]) # Generando wi 
    b<-1/(w+1) 
    t<-1/(w+f(0,mu))
    pi<-f(C,mu)/(1-ppois(C-1,mu)) # Phi teórico
    U1<- t(Z)%*%(di*w*t*I - w*b) 
    U2<- t(X)%*%(-di*f(0,mu)*mu*t*I+di*(Y-mu)*I. + (1-di)*C*pi)
    as.vector(rbind(U2,U1))} # Vector de derivadas total
  
  # Matriz Hessiana
  Hess<-function(beta){
    mu<-exp(X%*%beta[1:P]) # Generando mui
    w<- exp(Z%*%beta[(P+1):(P+s)]) # Generando wi 
    b<-1/(w+1) 
    t<-1/(w+f(0,mu))
    pi<-f(C,mu)/(1-ppois(C-1,mu)) # Phi teorico
    R<-diag(as.vector(di*f(0,mu)*w*t^2*I-w*b^2))
    J11<-t(Z)%*%R%*%Z # Segundas derivadas para gamma
    S<-diag(as.vector(-di*(f(0,mu)*mu*(w+f(0,mu)-w*mu)*t^2*I+mu*I.)+
                        (1-di)*C*((C-mu)*pi-C*pi^2)))
    J22<-t(X)%*%S%*%X # Segundas derivadas para beta
    K<-diag(as.vector(di*f(0,mu)*mu*w*t^2*I))
    J12<-t(Z)%*%K%*%X # Segundas derivadas cruzadas
    as.matrix(rbind(cbind(J22,t(J12)),cbind(J12,J11)))} # Matriz Hessiana
  
  # Funcion log-likelihood 
  log.L<-function(beta){
    sum(di*(log(exp(Z%*%beta[(P+1):(P+s)])+f(0,exp(X%*%beta[1:P])))*I+log(f(Y,exp(X%*%beta[1:P])))*I.)
        +(1-di)*log(1-ppois(C-1,exp(X%*%beta[1:P])))-log(1+exp(Z%*%beta[(P+1):(P+s)])))}
  summary(maxLik(log.L,grad,Hess,start = beta))
}

################################################################################
#                              Variable Selection
################################################################################
#----
DATA_MATE3171<-DATA_MATE3171[,c(1,4:6,9:12,20,21)]
DATA_MATE3171<-na.omit(DATA_MATE3171)
DATA_MATE3171 <- DATA_MATE3171 %>%filter(IGS != 0,`APR MATE`!= 0,`GPA ESC SUP`!= 0)
attach(DATA_MATE3171)
################################################################################
#                                     AIC           
################################################################################
#----
DATA_MATE3171<-DATA_MATE3171[,c(1:4,6:8,10)] 
DATA_MATE3171$`TIPO ESC`<-relevel(as.factor(DATA_MATE3171$`TIPO ESC`),ref="PUBLICA") # (Optional) Reference level public 
DATA_MATE3171$`GENERO`<-as.factor(DATA_MATE3171$`GENERO`)
DATA_MATE3171$`ANO ADMISION`<-as.factor(DATA_MATE3171$`ANO ADMISION`)
attach(DATA_MATE3171)
################################################################################
#                               Correlation
cor(data.matrix(DATA_MATE3171))
################################################################################
################################################################################
#                             Initial values 
################################################################################
#----
fit<-zeroinfl(Fallas ~.|.-IGS, data=DATA_MATE3171 )
summary(fit)
betanew<-coef(fit)
################################################################################
#                             Hessian
################################################################################
#----
status<-1-Censura
X<-model.matrix(~.-Fallas,data=DATA_MATE3171) # Matrix para el modelo Poisson  
Z<-model.matrix(~.-Fallas-IGS,data=DATA_MATE3171) # Matrix para el modelo binario
fit<-Zipc.maxlik(X,Z,Fallas,status,betanew)
fit
################################################################################
#                             Prediction
################################################################################
#----
beta<-as.vector(coef(fit)[,1])
P<-ncol(X);s<-ncol(Z);n<-nrow(X)
mu<-exp(X%*%beta[1:P]) # Prediccion tipo count
pii<- plogis(Z%*%beta[(P+1):(P+s)]) # Prediccion tipo zero
resp<-(1-pii)*mu # Prediccion tipo respuesta 


predic<-data.frame(mu,pii,resp,status,Fallas)
predic<-predic %>% filter(status > 0) # Utilizando solamente los datos que no estan censurados
resp<-predic$resp;mu<-predic$mu;pii<-predic$pii;Fallas<-predic$Fallas
y <-0:max(Fallas)  # Prediccion tipo probabilidad
prob <- matrix(NA, nrow = length(pii), ncol = length(y))
for (i in 1:length(pii)) {
  prob[i,]<-(pii[i] * I(y == 0) + (1-pii[i]) * dpois(y, mu[i]))
}

################################################################################
#                           Confusion Matrix
################################################################################
prob
obs.zero <- Fallas<1;preds <-prob[,1]
P.corte<-0.7
preds.zero <- preds >= P.corte
confusionMatrix(as.factor(preds.zero),as.factor(obs.zero))
################################################################################
#                                P(y<2)
################################################################################
obs.zero <- Fallas<2;preds <-prob[,1]+prob[,2]
P.corte<-0.4
preds.zero <- preds >= P.corte
confusionMatrix(as.factor(preds.zero),as.factor(obs.zero))
################################################################################
#                          Confusion Matrix
################################################################################
U.corte<-0.8
predsmu.zero <- resp <= U.corte
table(predsmu.zero)
confusionMatrix(as.factor(predsmu.zero),as.factor(obs.zero))
