##%######################################################%##
#                                                          #
####     Script Análisis exploratorios Introducción     ####
####            a la estadística inferencial            ####
#                                                          #
##%######################################################%##


# R.1 ---------------------------------------------------------------------
 
# ESTE EJEMPLO ES DE DATOS DE EVERITT 04 P. 17. CONSISTE EN INFORMACIÓN SOBRE CONTAMINACIÓN AMBIENTAL EN EU EN ZONAS METROPOLITANAS.

# llamo los datos (ojo que Everitt los tiene en formato dat, si están como txt, hay que llamarlos con read.table)

airpoll<-source("chap2airpoll.dat")$value
ap <- data.frame(airpoll)
write.csv(ap, "apdf.csv")
adf <- read.csv("airpoll.csv", header = T, sep = ";") #tabla con dato faltante


attach(airpoll)
airpoll
names(airpoll)


# Exploración Univariada --------------------------------------------------

# Comenzamos por ver el vector de medias  y varianzas

#mean(airpoll)
#sd(airpoll)^2

summary(airpoll)
summary(airpoll$SO2)

# vease la diferencia entre la media y la mediana para reconocer desviaciones, calculese el intervalo intercuartiles (3er-1er).

windows()
boxplot(SO2, range=0, ylab="SO2") # en este caso, los "bigotes" del boxplot ubican el máximo (1) y el mínimo (278).

boxplot(SO2, ylab="SO2") #en este caso, la función se ejecuta con range = 1.5 por defecto.

 iqSO2<-69-11
 iqSO2

# Una buena regla de dedo para identificar datos atípicos es: que los puntos que caen mas allá del 3er+1.5(intercuartil) o mas bajo que 1er-1.5(intercuartil) son valores atípicos.
 
atipicossup<-69+(iqSO2*1.5)
  atipicossup

atipicosinf<-abs(11-(iqSO2*1.5))
  atipicosinf

  
hist(SO2,lwd=2)
abline(v = 156, col = "blue")


airpoll

# ¿Cuales son los valores atípicos para SO2?
# Ahora veamos lo que considera R como atípicos por default 

par(mfrow=c(1,3))
boxplot(SO2, range=0, ylab="SO2")
boxplot(SO2, ylab="SO2")
boxplot(SO2, range=1.5, ylab="SO2")
  
#Veamos las distribuciones de todas

par(mfrow=c(3,3))
hist(SO2,lwd=2); abline(v = c(53.77, 30), col = c("blue", "red"))
hist(Rainfall,lwd=2)
hist(Education,lwd=2)
hist(Popden,lwd=2)
hist(Nonwhite,lwd=2)
hist(NOX,lwd=2)
hist(Mortality,lwd=2)

# ¿reconocen desviaciones negativas o positivas? Son normales?

par(mfrow=c(3,3))                                                 
qqnorm(SO2, main="Q-Q plot SO2"); qqline(SO2, col = 2, lty = 2)
qqnorm(Rainfall, main="Q-Q plot Rainfall"); qqline(Rainfall, col = 2, lty = 2)
qqnorm(Education, main="Q-Q plot Education"); qqline(Education, col = 2, lty = 2)
qqnorm(Popden, main="Q-Q plot Popden"); qqline(Popden, col = 2, lty = 2)
qqnorm(Nonwhite, main="Q-Q plot Nonwhite"); qqline(Nonwhite, col = 2, lty = 2)
qqnorm(NOX, main="Q-Q plot NOX"); qqline(NOX, col = 2, lty = 2)
qqnorm(Mortality, main="Q-Q plot Mortality"); qqline(Mortality, col = 2, lty = 2)


# Relaciones bivariadas ---------------------------------------------------
# Veamos que relación hay entre las distintas variables. Aquí utilizo una función smoooth (regresión con pesos locales) que permite sugerir con los propios datos que tipo de relación pudieran tener. 

pairs(airpoll, panel=panel.smooth)

#veamos con mas detalle la relación SO2-mortalidad

nombres<-abbreviate(row.names(airpoll))
par(mfrow=c(1,1))
plot(SO2,Mortality,lwd=2,type="n")
text(SO2,Mortality,labels=nombres,lwd=2)

detach(airpoll)

# C.1


# R.2 valores faltantes ---------------------------------------------------

airpoll
airpoldf <- read.table("datofalta.txt")
airpoldf
attach(airpoldf)

# Lo mas fácil la media 

summary(Mortality)
sum(is.na(Mortality))
    
# Cual es el valor imputado? Cuales son los problemas asociados a esta imputación?

# regresión mortalidad y SO2  

par(mfrow=c(1,1)) 
plot(SO2,Mortality,lwd=2)
abline(v = 59, h = 940.2)
abline(v = 59, h = 921.9, col = "red")

regmort<-lm(Mortality~SO2)
summary(regmort)
m <-  (915.4720997 + (0.4266209*59))    


      

abline(lm(Mortality~SO2))

predict(regmort, list(SO2=58)) 
      
par(mfrow=c(1,2))
plot(SO2,Mortality,lwd=2) 
abline(regmort)
abline(v = 58, h = c(940.2161, 954.4211), col = "red")
      
plot(logSO2,logM,lwd=2) 
abline(lm(logM~logSO2))
abline(v = 58, h = 940.2161, col = "red")
      
      
logM<-log(Mortality)
logSO2<-log(SO2+7)
loglog<-lm(logM~logSO2)
summary(loglog)
plot(logSO2,logM,lwd=2) 
abline(lm(logM~logSO2))
       
#el valor de SO2 que corresponde al valor faltante de mortalidad es 58. Como hemos generado un modelo de logaritmos a ambos lados de la ecuación sacamos el log del (SO2+7)
 
log(58+7)

#Usamos la función predict para predecir el valor correspondiente de Mortalidad

predict(loglog, list(logSO2=4.174387)) 

# pero recordando que usamos logaritmos en el modelo, retrotransformamos con el antilog con base e (e elevado al numero que nos interesa retro transformar)

exp(6.861105)

# El valor predicho por regresión lineal es?   Cuales son los problemas  asociados a esta imputación?

#fin





