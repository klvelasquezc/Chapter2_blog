##%######################################################%##
#                                                          #
####                    Correlación                     ####
#                                                          #
##%######################################################%##


# R.1 ---------------------------------------------------------------------

# Vamos a ver una serie de datos para ver si existe una relación lineal entre ellos.

data<-read.table("twosample.txt",header=T)
attach(data)
data
plot(x,y)


# Se acuerdan que necesitamos primero para calcular el coeficiente de correlación de pearson?
# Las varianzas individuales

var(x)
var(y)

# ¿y que mas?
# La covarianza y estamos hechos

var(x,y)

# ahora calculamos r

var(x,y)/sqrt(var(x)*var(y))

# Ahora hagamoslo en automático

cor(x,y)

# y ahora hagamos la prueba de hipótesis

# Calculamos EE de r

EEr<-((1-(cor(x,y)^2))/(length(x)-2))^0.5
EEr

# Calculo t de la muestra

te<-cor(x,y)/EEr
te

# calculo t de tablas

qt(0.975,47)

#calculo la p

2*(1-pt(18.67914,47))

#ahora hagamoslo de manera automática

pearson<-cor.test(x,y)
pearson

# ¿que nos falta?. Pues no sabemos si cumplimos con los supuestos. Veamos el de normalidad

par(mfrow=c(1,2))
qqnorm(x, main="Q-Q plot x"); qqline(x, col = 2, lty = 2)
qqnorm(y, main="Q-Q plot y"); qqline(y, col = 2, lty = 2)

# que opciones tengo.

# 1. Hacer una prueba de sesgo y kurtosis para ver si estas desviaciones son significativas

# 2. Si son significativas, puedo intentar transformaciones o puedo utilizar muchas de las otras pruebas de correlación que son  robustas a la violación de este supuesto.
# vean Q y k p.76 y Crawley p.97-102.

# Fin


