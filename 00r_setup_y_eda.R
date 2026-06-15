# ============================================================
# PROYECTO REGRESIÓN - GRUPO 9 - FOREST FIRES
# BLOQUE 0: Setup, carga de datos y EDA formal
# EYP2307 - Análisis de Regresión 2026
# ============================================================
# PROPÓSITO: Producir estadísticos formales que respalden
# decisiones tomadas visualmente en el avance.
# FEEDBACK ATENDIDO: se agrega skewness, kurtosis, y se llama
# al suavizador loess por su nombre correcto.
# ============================================================

# -- Paquetes ------------------------------------------------
# Instalar si no están disponibles:
#install.packages(c("car","lmtest","nortest","moments","MASS"))
#install.packages("stargazer")
library(car)      # vif(), crPlots(), avPlots()
library(lmtest)   # bptest(), dwtest()
library(nortest)  # lillie.test() - Kolmogorov-Smirnov con corrección Lilliefors
library(moments)  # skewness(), kurtosis()
library(MASS)     # boxcox(), stepAIC()
library(stargazer)# para obtener los códigos en latex en la consola
set.seed(2307)    # Reproducibilidad

# -- Carga de datos ------------------------------------------
datos <- read.csv("data/forestfires.csv")

# -- Codificación de factores (orden cronológico) ------------
datos$month <- factor(datos$month,
  levels = c("jan","feb","mar","apr","may","jun",
             "jul","aug","sep","oct","nov","dec"))
datos$day   <- factor(datos$day,
  levels = c("mon","tue","wed","thu","fri","sat","sun"))

# -- Transformación logarítmica de la respuesta --------------
# Justificación formal (ver sección 3.1 del informe):
# 1) La distribución de area es fuertemente asimétrica
# 2) 47.8% de observaciones tienen area = 0
# 3) La transformación log(y+1) estabiliza varianza y acerca
#    la distribución condicional a la normal requerida para
#    inferencia exacta (supuesto A5 del modelo lineal clásico)
datos$log_area <- log(datos$area + 1)

# -- Verificación de NAs ------------------------------------
stopifnot(sum(is.na(datos)) == 0)
# incluir tabla en el informe
tabla1 <- cat("Filas:", nrow(datos), "| Columnas:", ncol(datos), "| NAs:", sum(is.na(datos)), "\n")
# codigo en latex de la tabla

# Generar el código LaTeX
stargazer(tabla1, type = "latex", 
          title = "Modelos de Regresión para Áreas Quemadas",
          label = "tab:regresiones")


# ============================================================
# BLOQUE 1: EDA FORMAL CON ESTADÍSTICOS
# ============================================================

# 1.1 Distribución de la respuesta: estadísticos formales -----
cat("\n--- Estadísticos formales de la respuesta ---\n")

stats_area    <- c(media   = mean(datos$area),
                   mediana = median(datos$area),
                   DE      = sd(datos$area),
                   skew    = skewness(datos$area),
                   kurt    = kurtosis(datos$area),
                   prop0   = mean(datos$area == 0))

stats_logarea <- c(media   = mean(datos$log_area),
                   mediana = median(datos$log_area),
                   DE      = sd(datos$log_area),
                   skew    = skewness(datos$log_area),
                   kurt    = kurtosis(datos$log_area),
                   prop0   = mean(datos$log_area == 0))

tabla_resp <- rbind("area" = stats_area, "log(area+1)" = stats_logarea)
print(round(tabla_resp, 3))

# Un skewness(area) > 3-4 y kurtosis >> 3 justifica la transformación
# Skewness cercano a 0 en log(area+1) indica mejora clara (nose si sale en clases)

## Proporción de ceros 

prop_cero <- mean(datos$area == 0)
cat(sprintf("\nProporción con area = 0: %.1f%%\n", prop_cero * 100))

# 1.2 Justificación formal de excluir "rain" -----------------
cat("\n--- Variable rain: justificación de exclusión ---\n")
cat("Obs. con rain > 0:", sum(datos$rain > 0), "de", nrow(datos),
    "=", round(100*mean(datos$rain > 0), 1), "%\n")
cat("Varianza de rain:", round(var(datos$rain), 6), "\n")
# Con Snn ≈ 0, la varianza del coeficiente Var(β̂_rain|X) = σ²/[Srr(1-R²r)]
# se dispara (segun el feedback esto está bien asique hay que dejarlo). 
# Esto se reflejará en un VIF muy alto si se incluye.
# Además: Breusch-Pagan con rain == 0 casi siempre → su "efecto" no es estimable
# con precisión razonable. Se excluye ANTES de ajustar.


# 1.3 Justificación de excluir "day" -------------------------
# (Se verificará formalmente con F-test en el Bloque 3)
cat("\n--- Distribución de incendios por día ---\n")
print(table(datos$day))
# Argumento: no hay mecanismo físico por el que el día de la semana
# cambie el comportamiento del fuego. si es que hay alguna 
# variación puede deberse al  registro. Se verifica con F-test (pendiente).)

# Test de levene, para comprobar que las varianzas de day son homogeneas. 
levene_day <- leveneTest(log_area ~ day, data = datos)
print(levene_day)


# 1.4 Descriptivos de predictores numéricos ------------------
cat("\n--- Estadísticos de predictores numéricos ---\n")
pred_num <- c("FFMC","DMC","DC","ISI","temp","RH","wind")
tab_pred <- t(sapply(pred_num, function(v) {
  x <- datos[[v]]
  c(min    = min(x),
    max    = max(x),
    media  = round(mean(x), 2),
    DE     = round(sd(x), 2),
    skew   = round(skewness(x), 3))
}))
print(tab_pred)
# ISI tiene skewness > 2 → candidato a transformación logarítmica (al final no 
# es necesario porque (sesgo = 2.52): la transformación log(ISI+1) prácticamente
# no mejora la correlación (−0.010 → +0.009). El argumento para transformarlo igual
# es de corrección del modelo, no de R². ISI tiene valores extremos (hasta 56, 
# pero la mayoría están bajo 15) que generan puntos de alto leverage. 
# Reducir ese efecto mediante log mejora la estabilidad del estimador β̂ aunque
# no suba el R².)
# Verificar en partial plots si la no-linealidad importa


# 1.5 Resumen de log_area por mes (Tabla 2 del avance mejorada)--
cat("\n--- Media de log(area+1) por mes ---\n")
res_mes <- tapply(datos$log_area, datos$month, function(x)
  c(n = length(x), media = round(mean(x), 3), DE = round(sd(x), 3)))
print(do.call(rbind, res_mes))
# Los meses con n=1 o n=2 (nov: 1, jan: 2, may: 2) son problemáticos (lo pusimos 
# en el avance) → argumento para agrupar en estaciones (ver Bloque 2)

# Test de levene, para comprobar que las varianzas de day son homogeneas. 
levene_mes <- leveneTest(log_area ~ month, data = datos)
print(levene_mes)



# 1.6 Matriz de correlación (predictores numéricos + respuesta)--
cat("\n--- Correlaciones con log(area+1) y entre predictores ---\n")
vars_cor <- c("FFMC","DMC","DC","ISI","temp","RH","wind","log_area")
cor_mat  <- cor(datos[, vars_cor], use = "complete.obs")
print(round(cor_mat, 3))
# Pares con |r| > 0.5: posible colinealidad a vigilar con VIF (ya tenemos esta
# tabla en el avance)

# Crear la carpeta "figuras" automáticamente si no existe
if (!dir.exists("figuras")) {
  dir.create("figuras")
}


pdf("figuras/fig1_distribucion_respuesta.pdf", width = 9, height = 4)
dev.off() 


# 1.7 Gráfico: distribución de la respuesta ------------------
pdf("figuras/fig1_distribucion_respuesta.pdf", width = 9, height = 4)
par(mfrow = c(1, 2))
hist(datos$area,
     main = "Figura 1a: Distribución de area",
     xlab = "Área quemada (ha)",
     ylab = "Frecuencia",
     col  = "#4472C4", border = "white", breaks = 40)
text(500, 200,
     paste0("Sesgo = ", round(skewness(datos$area), 2),
            "\nCurt. = ", round(kurtosis(datos$area), 2)),
     cex = 0.85)

hist(datos$log_area,
     main = "Figura 1b: Distribución de log(area+1)",
     xlab = "log(área + 1)",
     ylab = "Frecuencia",
     col  = "#ED7D31", border = "white", breaks = 30)
text(3, 80,
     paste0("Sesgo = ", round(skewness(datos$log_area), 2),
            "\nCurt. = ", round(kurtosis(datos$log_area), 2)),
     cex = 0.85)
par(mfrow = c(1, 1))
dev.off()
cat("Figura 1 guardada en figuras/fig1_distribucion_respuesta.pdf\n")


# 1.8 Gráfico: dispersión y suavizador loess ------------------
# Esta ya la tenemos no se cual será mejor
pdf("figuras/fig2_dispersion_predictores.pdf", width = 12, height = 8)
par(mfrow = c(2, 4))
for (v in pred_num) {
  plot(datos[[v]], datos$log_area,
       xlab = v, ylab = "log(area+1)",
       pch  = 19, cex = 0.4, col = "#00000055",
       main = paste("log(area+1) ~", v))
  abline(lm(datos$log_area ~ datos[[v]]), col = "#4472C4", lwd = 1.5)
  # Suavizador no paramétrico loess (NO es interpolación)
  lo <- loess(datos$log_area ~ datos[[v]], span = 0.75)
  ox <- order(datos[[v]])
  lines(datos[[v]][ox], fitted(lo)[ox], col = "#ED7D31", lwd = 1.5, lty = 2)
}
par(mfrow = c(1, 1))
dev.off()
cat("Figura 2 guardada en figuras/fig2_dispersion_predictores.pdf\n")
