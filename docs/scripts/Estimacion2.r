##%######################################################%##
#                                                          #
####                Script Estimación 2                 ####
####     Introducción a la estadística inferencial      ####
#                                                          #
##%######################################################%##


# R.1 ---------------------------------------------------------------------
# Estimación de parámetros de tendencia central ---------------------------

# Sabemos que aunque los datos no se agrupen alrededor de un valor típico si hacemos muestras consecutivas los estadísticos de estas van a tener una tendencia central.

# Veamos un cuerpo de datos. Una variable y.

yvals <- read.table("yvalues.txt",header=T)
attach(yvals)
yvals

#Una manera muy simple de medir la tendencia central es ver cual es el valor más frecuente. Este se denomina MODA

yord<-sort(y)
yord

windows()
hist(y)

# ¿cual es la clase modal aquí?

# Pero ahora queremos saber la media aritmética (el promedio) que es la suma de todos los valores dividido por n. ¿Que hago?

total<-sum(y)
sum(y)

# pero ahora necesito saber cuantos valores son

n<-length(y)
n
media<- total/n
media

# pero ahora quiero tener una función que me sirva para siempre

media.aritmetica <- function(x) {
	sum(x)/length(x) }

#ya está, ahora probémosla

data<-c(3,4,6,7)
media.aritmetica(data)

media.aritmetica(y)

#¡Que bien! ¡R es fantástico! puedo calcular la media siempre que me plazca.

# En realidad la mayor parte de las funciones estadísticas están ya construidas en R. Por supuesto que la media es una de ellas. Solo quería mostrarles que no hay nada obscuro detrás de los objetos ya creados en R

mean(y)

# La media como medida de tendencia central tiene el serio problema de que es muy sensible a valores atípicos. vean lo siguiente.

dataat<-c(data,100)
dataat
mean(dataat)

#comparado con

mean(data)

# Una alternativa es la mediana, que es el valor de en medio, una vez que todos los valores han sido ordenados. Veamos dataat

dataat

# ¿Cual es la mediana?

median(dataat)

# es mucho mejor estimación de el centro que 24.

#y de data?

data

median(data)

# y para y?

median(y)
mean(y)

# Se parecen mucho porque no hay valores atípicos y porque la distribución es simétrica.

#Ahora pensemos en fenómenos que cambian multiplicativamente. ¿Conocen alguno?

# Uno de los más comunes en ecología es el crecimiento poblacional y por lo tanto la dispersión de organismos de una población. En dichos casos la media aritmética y/o la mediana suelen ser pésimos estimadores de la tendencia central. Veamos un ejemplo.

#El número de insectos en una serie de plantas vecinas es

insectos<-c(1,10,1000,10,1)

# ¿Cual es la mejor estimación de tendencia central?

windows()
hist(insectos)

mean(insectos)
median(insectos)

# C.1


# R.2 ---------------------------------------------------------------------

#Lo que se usa es la media geométrica que si se acuerdan hay dos maneras de calcularla. Para el ejemplo de los insectos ¿cual es la mas simple?

100000^0.2

# ¿y la otra?

exp(mean(log(insectos)))

detach(yvals)
rm(insectos)
ls()

#C.2

# R.3 ---------------------------------------------------------------------

# Medidas de dispersión

# Veamos un cuerpo de datos cualquiera y preguntémonos cómo podemos medir su dispersión.

y<-c(13,7,5,12,9,15,6,11,9,7,12)

# veamos cómo se ve

windows()
plot(y,ylim=c(0,20))

#Lo más fácil es decir de donde a donde va (el intervalo).

range(y)

# Pero esto tiene sus problemas.

# 1. No tiene relación con el parámetro población de intervalo
# 2.incrementa con la n.
# 3.Además es muy susceptible a valores atípicos
# 4. no considera a todos lo valores.

# Otra medida de dispersión muy importante es la varianza. Que está fundamentada en las desviaciones (o residuales) de cada valor con la media

y<-c(13,7,5,12,9,15,6,11,9,7,12)
abline(mean(y),0)
for (i in 1:11) lines(c(i,i),c(y[i],mean(y)))

# se usa la suma de cuadrados de la diferencia de cada valor con la media general como base. ¿Cómo lo calculamos?

y-mean(y)
(y-mean(y))^2
sum((y-mean(y))^2)

# Fantásitico, pero que sucede a SC cada vez que yo adiciono una nueva observación. ¿Que tenemos que hacer?

# Si divido entre n, se llama la desviación media de los cuadrados.

#C.3

# R.4 ---------------------------------------------------------------------

variance <- function (x)   sum((x-mean(x))^2)/(length(x)-1)
variance(y)

# Pero claro, ya está definido.

var(y)

#La relación entre la varianza de la muestra y el tamaño de muestra (n)

windows()
plot(c(0,32),c(0,15),type="n",xlab="Tamaño de muestra",ylab="Varianza")

# Lo que vamos a hacer es seleccionar aleatoriamente números de una población que tiene una distribución normal (media 10 y var 4). Esto lo vamos a hacer repetidas veces pero nuestra muestra va a ir incrementando su n desde 3 hasta 31. Vamos a sacar 30 muestras de cada tamaño de muestra. Es decir 30 muestras de 3 números, 30 muestras de 5 números etc. A cada muestra le vamos a calcular su varianza y las vamos a graficar.

for (tm in seq(3,31,2)) {
for( i in 1:30){
x<-rnorm(tm,mean=10,sd=2)
points(tm,var(x)) }}

# Ahora pueden ver que la varianza poblacional puede estar muy mal estimada con tamaños de muestra pequeños. Y a medida que aumentamos el tamaño de muestra la probabilidad de que la estimación está muy lejos del parámetro disminuye. Esta es una razón más para elegir muy cuidadosamente el tamaño de muestra. En unas clases vamos a hablar de esto en el contexto de pruebas de hipótesis.

#C.4

# R.5 ---------------------------------------------------------------------

#Hagamos el mismo proceso que hicimos antes para elegir muestras con tamaños de muestra que incrementan pero ahora veamos que sucede con nuestra medida de desconfianza de la estimación a medida que aumenta n

windows()
plot(c(0,32),c(0,2),type="n",xlab="Tamaño de muestra",ylab="Error estandar")
 for (tm in seq(3,31,2)) {
for( i in 1:30){
x<-rnorm(tm,mean=10,sd=2)
points(tm,sqrt(var(x)/tm)) }}

# Veamos un ejemplo de la concentración de ozono en unos invernaderos

ozono<- read.table("gardens.txt",header=T)
attach(ozono)
ozono

windows()
par(mfrow=c(1,3))
plot (gardenA)
plot (gardenB)
plot (gardenC)

MediaA<- mean(gardenA)
EEA<- sqrt(var(gardenA)/10)
MediaB<- mean(gardenB)
EEB<- sqrt(var(gardenB)/10)
MediaC<- mean(gardenC)
EEC<- sqrt(var(gardenC)/10)

MediaA
EEA
MediaB
EEB
MediaC
EEC

# ¿Para cual invernadero puedo yo confiar más de la estimación de la media?

# C.5

# R.6 ---------------------------------------------------------------------

# Vamos a calcular los intervalos de confianza para la media de los invernaderos usando la distribución de t. qt nos d? el valor de t para el cual hay cierta proporción a la izquierda.

qt(.975,9)

#calculemos los intervalos de confianza al 95% para el invernadero B

qt(.975,9)*sqrt(1.3333/10)

# el reporte diario, la concentración promedio de ozono en el invernadero B fue de 5.0?0.826(I.C.95%, n=10)

windows()
plot(gardenB)
abline(mean(gardenB),0)
abline((mean(gardenB)+0.826),0, lty=2)
abline((mean(gardenB)-0.826),0, lty=2)

# los dejo que ustedes calculen aquellos de los invernaderos A y C.

#C.6

# Fin ---------------------------------------------------------------------

