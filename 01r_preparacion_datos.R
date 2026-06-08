# ============================================================
# BLOQUE 2: PREPARACIÓN PARA EL ANÁLISIS DE MODELOS
# - División entrenamiento/prueba 80/20
# - Agrupación de meses en estaciones
# - Decisión sobre "day" (F-test formal)
# ============================================================
# NOTA CRÍTICA sobre step() y overfitting (FEEDBACK ATENDIDO):
#   Todo el proceso de selección de variables se ejecuta
#   EXCLUSIVAMENTE sobre el conjunto de entrenamiento.
#   El conjunto de prueba se reserva únicamente para evaluar
#   el desempeño final. Esto evita que el criterio AIC/BIC
#   del paso a paso se "contamine" con información del test.
# ============================================================

# Requiere haber cargado "datos" desde 00_setup_y_eda.R
# source("00_setup_y_eda.R")

set.seed(2307)

# -- 2.1 División entrenamiento / prueba (80/20) ------------
n         <- nrow(datos)           # 517
n_train   <- floor(0.80 * n)       # 413
idx_train <- sample(seq_len(n), size = n_train, replace = FALSE)

train <- datos[ idx_train, ]
test  <- datos[-idx_train, ]

#esta info debemos incluirla como tabla o escrita para que tome menos espacio

cat("Conjunto entrenamiento:", nrow(train), "obs.\n")
cat("Conjunto prueba:       ", nrow(test),  "obs.\n")
cat("Proporción area==0 en train:", round(mean(train$area == 0), 3), "\n")
cat("Proporción area==0 en test: ", round(mean(test$area  == 0), 3), "\n")


# -- 2.2 Agrupación de meses en estaciones ------------------
# Las razones son las siguientes:
#   a) ESTADÍSTICA: enero (n=2), mayo (n=2), noviembre (n=1)
#      tienen estimaciones completamente inestables como dummies
#      individuales. Var(β̂_mes|X) = σ²/(n_mes × ...) se dispara.
#   b) Desde un punto de vista práctico: en Portugal, el régimen de incendios 
#      responde al ciclo estacional (verano seco/caluroso = máximo riesgo,
#      invierno húmedo = mínimo riesgo).
#   c) PARSIMONIA: 4 dummies de estación vs 11 de mes → ahorra
#      7 grados de libertad (del denominador de F) sin perder
#      estructura relevante.
#   La equivalencia estadística se verificará con F-test.

asignar_estacion <- function(mes) {
  # Hemisferio norte: Portugal
  verano   <- c("jun","jul","aug")
  otono    <- c("sep","oct","nov")
  invierno <- c("dec","jan","feb")
  primavera <- c("mar","apr","may")
  ifelse(mes %in% verano,    "verano",
  ifelse(mes %in% otono,     "otono",
  ifelse(mes %in% invierno,  "invierno",
                             "primavera")))
}

# Aplicar a los dos conjuntos
for (df_name in c("datos", "train", "test")) {
  d <- get(df_name)
  d$estacion <- factor(asignar_estacion(as.character(d$month)),
                        levels = c("invierno","primavera","verano","otono"))
  assign(df_name, d)
}

cat("\nDistribución por estación (datos completos):\n")
print(table(datos$estacion))
cat("Distribución por estación (train):\n")
print(table(train$estacion))


# -- 2.3 Test F: ¿month completo vs estacion? ---------------
# Si el modelo con estaciones no es significativamente peor que
# el de 12 meses, preferimos la versión parsimonosa.
# (Modelos anidados: estación ⊂ month)

m_tmp_est <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + estacion,
                data = train)
m_tmp_mes <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + month,
                data = train)

cat("\n--- F-test: modelo con estaciones vs. con meses individuales ---\n")
print(anova(m_tmp_est, m_tmp_mes))
# Como el p-value es > 0.05: no se rechaza H0: que el modelo con estaciones
# es suficiente para predecir → preferimos estación por parsimonia.
# Si fuera p-valor < 0.05: los meses individuales aportan → considerar
# mantener month


# -- 2.4 Test F: ¿incluir "day"? Mismo test que el anterior aunque sabemos que
# no tiene razón física por la que el día influya en el area quemada
m_tmp_sin_day <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + estacion,
                    data = train)
m_tmp_con_day <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + estacion + day,
                    data = train)

cat("\n--- F-test: ¿day aporta al modelo? ---\n")
print(anova(m_tmp_sin_day, m_tmp_con_day))
# Hipótesis: H0: β_lunes = β_martes = ... = β_domingo = 0
# Como p-valor > 0.05: no rechazamos H0 o sea day no aporta significativamente.
# Argumento adicional: no hay mecanismo físico que justifique
# que el día de la semana afecte el comportamiento del fuego.



# -- 2.5 ¿Transformar ISI? ----------------------------------
# ISI tiene sesgo positivo. Verificar si log(ISI+1) mejora
# la linealidad con log_area. Esto se evaluará formalmente
# con crPlots (preguntaré em el foro si lo podemos usar, en el código 04 están
# las razones de porque sirve)

cat("\nSesgos antes de entrar al modelo:\n")
cat("Sesgo ISI:  ", round(moments::skewness(train$ISI), 3), "\n")
cat("Sesgo DMC:  ", round(moments::skewness(train$DMC), 3), "\n")
cat("Sesgo DC:   ", round(moments::skewness(train$DC),  3), "\n")
cat("Sesgo FFMC: ", round(moments::skewness(train$FFMC),3), "\n")
# Si alguno tiene sesgo > 2 (también preguntaré en el foro), considerar transformación log.
# IMPORTANTE: transformar predictores NO viola los supuestos del
# modelo lineal clásico (el modelo sigue siendo lineal en parámetros).
# Solo cambia la interpretación de β_j (esto se habló en clases) para que no se nos
# olvide ponerlo en el informe.
