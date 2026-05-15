#===============================================================================
# PROYECTO - Análisis de Regresión
# Área Quemada de Incendios Forestales
# Integrantes: 
#===============================================================================

# 1. PREPARACIÓN DE DATOS
#-------------------------------------------------------------------------------
# Aseguramos el formato de factores para las variables categóricas
datos <- read.csv("data/forestfires.csv")
datos$month <- as.factor(datos$month)
head(datos)

# 2. ANÁLISIS EXPLORATORIO DE DATOS (EDA)
#-------------------------------------------------------------------------------

summary(datos)

plot(datos$area, type = "l", 
     main = "Área quemada en el tiempo",
     ylab = "Área", xlab = "Índice Temporal (Meses)",
     col = "blue", lwd = 1.5)





