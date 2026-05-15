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



