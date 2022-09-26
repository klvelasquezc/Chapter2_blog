##%######################################################%##
#                                                          #
####      Prueba de Hipótesis de una población bis      ####
#                                                          #
##%######################################################%##

# C.2
# Vamos a escribir una función para calcular sesgo.

skew<-function(x){
m3<-sum((x-mean(x))^3)/length(x)
s3<-sqrt(var(x))^3
m3/s3  }


# veamos los datos skewdata

data<-read.table("skewdata.txt",header=T)
attach(data)
names(data)
hist(values)

# Calculemos el sesgo

skew(values)

# Ahora hagamos una prueba de t para determinar si hay sesgo respecto a lo esperado para una población de inferencia con distribución normal

skew(values)/sqrt(6/length(values))
1-pt(2.949,28)

# ¿Cual es la conclusión?.
# Efectivamente está más sesgado de los esperado normalmente, ¿que se puede hacer, si insistimos en cumplir con la normalidad?

skew(sqrt(values))/sqrt(6/length(values))
skew(log(values))/sqrt(6/length(values))


kurtosis<-function(x) {
m4<-sum((x-mean(x))^4)/length(x)
s4<-var(x)^2
m4/s4 - 3  }
kurtosis(values)
kurtosis(values)/sqrt(24/length(values))
