##%######################################################%##
#                                                          #
####              Pruebas de hipótesis Una              ####
####      población Introducción a la estadística       ####
#                                                          #
##%######################################################%##


# R.1 ---------------------------------------------------------------------

# Veamos un ejemplo de prueba de hipótesis para una población (Crawley p.64). Estos son datos de Michelson (1978), de medidas tomadas para estimar la velocidad de la luz. Que ahora sabemos es cercano (299,990 km p seg).
 
veluz<-read.table("light.txt",header=T)
attach(veluz)
names(veluz)
hist(speed)
summary(speed)

# A todos los valores se les ha restado 299 000 para facilitar su visualización ¿que pueden ver aquí? ¿tenemos valores atípicos?

boxplot(speed)

qqnorm(speed)
qqline(speed,lty=2)

#La muestra no se distribuye de acuerdo a lo esperado normalmente (cosa que una prueba de t para una población asume), pero además ¿Creen que la población de referencia se distribuya de manera normal?

# no puede porque tiene el problema de que la vel de la luz no puede tomar valores negativos

# Nuestra hipótesis es que los datos de Michelson difieren de el valor prevaleciente en esa época como la vel de la luz 299,990 km/seg. Como a todos los valores les quitaron 299 000 entonces ¿cual va a ser la referencia?

# Ahora, si cumpliéramos con el supuesto de que la pob. de referencia se distribuye normalmente, cómo resolvemos este problema?...ustedes ya lo saben hacer.

# Primero. ¿Este problema tiene una o dos colas que le pisen?
# Entonces ¿que harían?

# Bueno, pero como no cumplimos con el primero de los supuestos, necesitamos otra alternativa. Que se les ocurre?


# Por supuesto que también podemos usar una técnica de remuestreo con remplazo! ******Este es el punto 1. de la tarea de hoy. Hagan la prueba con bootstrap******

#Existe otra alternativa es una prueba llamada de rangos signados de Wilcoxon.
# C.1


# R.2 ---------------------------------------------------------------------

library(stats)

wilcox.test(speed,mu=990)

# La probabilidad de obtener la media de nuestra muestra (estadístico) en una población de las medias de muchas muestras con media 990 (parámetro) es del 0.2%. Como aceptamos una probabilidad de equivocarnos al rechazar una H0 cuando esta es cierta del 5% entonces, la rechazamos!


