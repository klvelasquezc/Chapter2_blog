##%######################################################%##
#                                                          #
####        Pruebas de hipótesis dos poblaciones bis    ####
#                                                          #
##%######################################################%##


# R.3 ---------------------------------------------------------------------

# Entonces ¿cual es la probabilidad de que mi muestra de frecuencias (o una mas extrema) pertenezca a una población donde el aislamiento de la planta y el plumaje del colibrí defensor son independientes?

1-pchisq(35.34,1)

# ¿Cual es el valor crítico para rechazar la H0 con un alfa de 0.05?

qchisq(0.95,1)

# Ahora vamos a hacerlo de manera automática en R

colibries<-matrix(c(38,14,11,51),nrow=2)
colibries
chisq.test(colibries)

# Noten que estos valores son un poco distintos porque se aplicó una corrección de Yates. Esta fue diseñada para frecuencias pequeñas, pero ahora existen pruebas mejores para frecuencias pequeñas (20% o más frecuencias menores a 5) como La Prueba exacta de Fisher Crawley p. 90 QyK p.388.

#Entonces le quito la corrección

chisq.test(colibries,correct=F)

# me da exactamente lo que calculamos a mano.

#C.3

#*****Aquí va la segunda parte de su tarea*****
# El problema se trata de un investigador que esta interesado en si las hormigas construyen preferentemente sus nidos en árboles de una de dos especies. Entonces muestreó 100 árboles de cada una de las  especies. Encontró que 60 árboles de la especie A  y 20 de la especie B tenían nidos.  ¿Esta muestra apoya la hipótesis de preferencia?


