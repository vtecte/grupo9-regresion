# ============================================================
#  EDA COMPLETO – Forest Fires (EYP2307 Grupo 9)
#  Avance: 27 de mayo de 2026
# ============================================================
# Justificación general del EDA:
#   Antes de ajustar cualquier modelo lineal Y = Xβ + ε,
#   el EDA cumple dos roles: (1) verificar que E(Y|X) sea
#   razonablemente modelable de forma lineal (supuesto A1),
#   y (2) detectar problemas que afecten supuestos A2–A5
#   (sesgo, varianza no constante, normalidad).
# ============================================================

# ── 0. Paquetes ──────────────────────────────────────────────
#install.packages(c("ggplot2", "GGally", "corrplot", "gridExtra", "car"))
library(ggplot2)
library(GGally)
library(corrplot)
library(gridExtra)
library(car)        # para vif() más adelante; útil ya en EDA

# ── 1. Carga y preparación ───────────────────────────────────
# (construido sobre EDA_GM.R existente)
datos_incendios <- read.csv("data/forestfires.csv")
datos_limpios   <- datos_incendios

# Factores con niveles ordenados (mes y día)
# Ordenamos las variables categóricas month/day que, dentro del
# modelo lineal, se codifican con indicadoras (dummies)

datos_limpios$month <- factor(datos_limpios$month,
                               levels = c("jan","feb","mar","apr","may","jun",
                                          "jul","aug","sep","oct","nov","dec"))
datos_limpios$day   <- factor(datos_limpios$day,
                               levels = c("mon","tue","wed","thu","fri","sat","sun"))

# Transformación logarítmica de la respuesta
# Justificación: 'area' tiene distribución fuertemente asimétrica con
# cola derecha y muchos ceros. log(area+1) estabiliza la varianza
# (heterocedasticidad, semana 9) y acerca la distribución condicional
# a la normalidad (supuesto A5, semana 3). Es la misma lógica del
# ejemplo de unidad quirúrgica (y' = ln(y)) visto en semana 8.
datos_limpios$log_area <- log(datos_limpios$area + 1)

# Verificación de NAs
cat("Total NAs:", sum(is.na(datos_limpios)), "\n")
# El dataset no tiene faltantes, pero se reporta de todas formas
# en el avance para cumplir con el punto 3 del trabajo esperado.

# ── 2. Resumen general ───────────────────────────────────────
cat("\n── Dimensiones ──\n")
cat("Observaciones:", nrow(datos_limpios), "| Variables:", ncol(datos_limpios), "\n")

cat("\n── Resumen numérico ──\n")
summary(datos_limpios[, c("FFMC","DMC","DC","ISI","temp","RH","wind","rain",
                           "area","log_area")])

# ── 3. Distribución de la variable respuesta ─────────────────
# Justificación: antes de cualquier modelo debemos entender la forma de
# Y (y de log Y). Una asimetría fuerte viola el supuesto A5 y puede
# generar heterocedasticidad (Cov(ε|X) ≠ σ²Iₙ, semana 9).

p1 <- ggplot(datos_limpios, aes(x = area)) +
  geom_histogram(bins = 40, fill = "#4E84C4", color = "white", alpha = 0.85) +
  labs(title = "Distribución de area (escala original)",
       x = "Área quemada (ha)", y = "Frecuencia") +
  theme_bw()

p2 <- ggplot(datos_limpios, aes(x = log_area)) +
  geom_histogram(bins = 35, fill = "#D55E00", color = "white", alpha = 0.85) +
  labs(title = "Distribución de log(area + 1)",
       x = "log(área + 1)", y = "Frecuencia") +
  theme_bw()

p3 <- ggplot(datos_limpios, aes(sample = log_area)) +
  stat_qq(color = "#4E84C4") + stat_qq_line(color = "red") +
  labs(title = "QQ-plot de log(area + 1)",
       x = "Cuantiles teóricos N(0,1)", y = "Cuantiles muestrales") +
  theme_bw()

# Proporción de observaciones con area = 0
prop_cero <- mean(datos_limpios$area == 0)
cat(sprintf("\nProporción con area = 0: %.1f%%\n", prop_cero * 100))
# Este valor es importante para el plan de modelamiento: una fracción alta
# de ceros puede justificar la transformación log(area+1) o, si es muy
# grande, discutir modelos alternativos como variable de mezcla.

grid.arrange(p1, p2, p3, ncol = 3)

# ── 4. Variables numéricas: distribuciones ───────────────────
vars_num <- c("FFMC","DMC","DC","ISI","temp","RH","wind","rain")

par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
for (v in vars_num) {
  hist(datos_limpios[[v]], main = v, xlab = v,
       col = "#4E84C4", border = "white", breaks = 25)
}
par(mfrow = c(1, 1))
# Qué buscamos: asimetría fuerte o valores extremos que luego
# se manifestarán como alto leverage (hᵢᵢ grande, semana 10).

# Estadístico de lluvia (rain): casi todos en 0 → poca variación
cat("\nFrecuencia de rain > 0:", sum(datos_limpios$rain > 0), "\n")
# Si rain tiene varianza casi cero, aportará poco al modelo
# (columna de X casi constante → problema de identificabilidad, semana 1)

# ── 5. Variables categóricas ─────────────────────────────────
p_mes <- ggplot(datos_limpios, aes(x = month, y = log_area, fill = month)) +
  geom_boxplot(show.legend = FALSE, outlier.alpha = 0.4) +
  labs(title = "log(area+1) por mes",
       x = "Mes", y = "log(área + 1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_dia <- ggplot(datos_limpios, aes(x = day, y = log_area, fill = day)) +
  geom_boxplot(show.legend = FALSE, outlier.alpha = 0.4) +
  labs(title = "log(area+1) por día",
       x = "Día", y = "log(área + 1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p_mes, p_dia, ncol = 2)

# investigar
# Justificación: si las medias condicionales difieren entre meses/días,
# incluir esas dummies mejora la especificación de E(Y|X) = Xβ.
# Esto es exactamente el ANOVA de un factor escrito como regresión
# (semana 4, sección 17).

# ── 6. Scatter plots: predictores vs log(area+1) ─────────────
# Justificación: el supuesto A1 exige linealidad en E(Y|X).
# Los scatter plots permiten detectar no linealidades obvias antes
# de ajustar (semana 7 en diagnóstico: "si hay curvatura sistemática,
# la media está mal especificada").

plots_scatter <- lapply(vars_num, function(v) {
  ggplot(datos_limpios, aes_string(x = v, y = "log_area")) +
    geom_point(alpha = 0.35, color = "#4E84C4", size = 1.2) +
    geom_smooth(method = "loess", color = "red", se = FALSE, linewidth = 0.8) +
    geom_smooth(method = "lm",    color = "black", se = FALSE,
                linetype = "dashed", linewidth = 0.8) +
    labs(title = paste("log(area+1) ~", v),
         x = v, y = "log(área+1)") +
    theme_bw(base_size = 9)
})
do.call(grid.arrange, c(plots_scatter, ncol = 4))
# Línea roja (loess) vs línea negra punteada (lineal):
# Si difieren mucho → considerar transformación del predictor.

# ── 7. Matriz de correlación entre predictores ───────────────
# Justificación: alta correlación entre predictores es síntoma
# de colinealidad (semana 10). Los VIF > 10 son problemáticos.
# Aquí hacemos el diagnóstico visual preliminar; en el informe
# final se calculan VIF formalmente.

vars_cor <- c("FFMC","DMC","DC","ISI","temp","RH","wind","rain","log_area")
R <- cor(datos_limpios[, vars_cor], use = "complete.obs")

corrplot(R, method = "color", type = "upper", tl.col = "black",
         tl.cex = 0.85, addCoef.col = "black", number.cex = 0.65,
         col = colorRampPalette(c("#D55E00","white","#4E84C4"))(200),
         title = "Correlación entre variables numéricas", mar = c(0,0,1.5,0))

# Pares más correlacionados (valor absoluto > 0.5):
cor_pairs <- which(abs(R) > 0.5 & abs(R) < 1, arr.ind = TRUE)
cor_pairs <- cor_pairs[cor_pairs[,1] < cor_pairs[,2], , drop = FALSE]
cat("\nPares con |r| > 0.5:\n")
for (i in seq_len(nrow(cor_pairs))) {
  v1 <- rownames(R)[cor_pairs[i,1]]
  v2 <- colnames(R)[cor_pairs[i,2]]
  cat(sprintf("  %s – %s: r = %.3f\n", v1, v2, R[cor_pairs[i,1], cor_pairs[i,2]]))
}

# ── 8. GGpairs para visualización completa ───────────────────
ggpairs(datos_limpios[, c("FFMC","DMC","DC","ISI","temp","RH","wind","log_area")],
        lower  = list(continuous = wrap("points", alpha = 0.3, size = 0.6)),
        diag   = list(continuous = wrap("densityDiag", color = "#4E84C4")),
        upper  = list(continuous = wrap("cor", size = 3)),
        title  = "Matriz de dispersión (predictores meteorológicos + log_area)")

# ── 9. Detección preliminar de outliers en la respuesta ──────
# Justificación: observaciones extremas en Y pueden ser influyentes
# si además tienen leverage alto (semana 10–11: Cook, DFFITS).
# Aquí los identificamos antes de modelar.

q_lo <- quantile(datos_limpios$log_area, 0.25)
q_hi <- quantile(datos_limpios$log_area, 0.75)
iqr  <- q_hi - q_lo
outliers_Y <- datos_limpios[datos_limpios$log_area > q_hi + 3 * iqr, ]
cat(sprintf("\nOutliers extremos en log_area (> Q3 + 3·IQR): %d observaciones\n",
            nrow(outliers_Y)))
if (nrow(outliers_Y) > 0) print(outliers_Y[, c("X","Y","month","day","log_area","area")])

# ── 10. Distribución espacial (coordenadas X, Y) ─────────────
# Justificación: X e Y son coordenadas de cuadrícula del parque.
# Si hay patrón espacial en los residuos → autocorrelación espacial
# (análoga a autocorrelación temporal, semana 11).
ggplot(datos_limpios, aes(x = X, y = Y, color = log_area, size = log_area)) +
  geom_jitter(alpha = 0.6, width = 0.15, height = 0.15) +
  scale_color_gradient(low = "lightyellow", high = "#D55E00",
                       name = "log(área+1)") +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  labs(title = "Distribución espacial de incendios",
       x = "Coordenada X (cuadrícula)", y = "Coordenada Y (cuadrícula)") +
  theme_bw()

# ── 11. Frecuencias de incendios por mes y día ───────────────
p_frec_mes <- ggplot(datos_limpios, aes(x = month)) +
  geom_bar(fill = "#4E84C4", alpha = 0.85) +
  labs(title = "Incendios por mes", x = "Mes", y = "Frecuencia") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_frec_dia <- ggplot(datos_limpios, aes(x = day)) +
  geom_bar(fill = "#D55E00", alpha = 0.85) +
  labs(title = "Incendios por día", x = "Día", y = "Frecuencia") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

grid.arrange(p_frec_mes, p_frec_dia, ncol = 2)
# Observación esperada: agosto y septiembre concentran la mayoría.
# Esto anticipa que la dummy de mes puede ser relevante.

# ── 12. Tabla de estadísticos descriptivos ───────────────────
desc_tabla <- data.frame(
  Media    = round(sapply(vars_num, function(v) mean(datos_limpios[[v]])), 2),
  Mediana  = round(sapply(vars_num, function(v) median(datos_limpios[[v]])), 2),
  DE       = round(sapply(vars_num, function(v) sd(datos_limpios[[v]])), 2),
  Min      = round(sapply(vars_num, function(v) min(datos_limpios[[v]])), 2),
  Max      = round(sapply(vars_num, function(v) max(datos_limpios[[v]])), 2),
  Asimetria = round(sapply(vars_num, function(v) {
    x <- datos_limpios[[v]]
    n <- length(x); m <- mean(x); s <- sd(x)
    (n / ((n-1)*(n-2))) * sum(((x - m)/s)^3)
  }), 2)
)
print(desc_tabla)

# ── MODELAMIENTO  ────────────────────────────────────────────
# ── 1. Pregunta de investigación ─────────────────────────────
# ¿Qué variables meteorológicas y de los índices FWI permiten
# predecir el área quemada (en escala log) en incendios del
# Parque Natural de Montesinho?
# Objetivo: PREDICTIVO (semana 2 del curso: "construir reglas
# para nuevos perfiles x0 e interesa desempeño fuera de muestra").

# ── 2. Modelo 0 – Sólo intercepto (referencia) ───────────────
# Justificación: es el modelo reducido de la prueba F global
# (semana 5). Todo modelo candidato se comparará contra éste
# mediante H0: β1 = ... = βp-1 = 0.
m0 <- lm(log_area ~ 1, data = datos_limpios)
cat("── Modelo 0 (intercepto) ──\n")
cat("SCE:", deviance(m0), "| gl:", df.residual(m0), "\n\n")

# ── 3. Modelo 1 – Índices FWI solamente ──────────────────────
# Justificación: FFMC, DMC, DC e ISI son los cuatro índices del
# sistema canadiense de peligro de incendios. Agruparlos en un
# modelo separado permite contrastar si el sistema FWI, en sí
# solo, tiene poder predictivo (prueba F parcial, semana 5).
m1 <- lm(log_area ~ FFMC + DMC + DC + ISI, data = datos_limpios)
cat("── Modelo 1 (FWI) ──\n")
print(summary(m1))

# Matriz de diseño X₁
cat("Dimensiones de X₁:", dim(model.matrix(m1)), "\n\n")

# ── 4. Modelo 2 – Variables meteorológicas ───────────────────
# Justificación: temp, RH, wind son los predictores "de campo"
# directamente medibles. Compararlos con los índices FWI permite
# evaluar si la información meteorológica básica es suficiente
# (modelos anidados dentro del modelo completo, semana 5 sección 7).
m2 <- lm(log_area ~ temp + RH + wind, data = datos_limpios)
cat("── Modelo 2 (meteorológicas) ──\n")
print(summary(m2))

# ── 5. Modelo 3 – Modelo completo con numéricas ──────────────
# Justificación: punto de partida para selección de variables
# (semana 8). Se ajusta con todos los predictores numéricos para
# calcular VIF (semana 10) y decidir cuáles mantener.
# rain se incluye aquí aunque casi siempre = 0; el VIF o la
# prueba t dirán si debe salir.
m3 <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + rain,
         data = datos_limpios)
cat("── Modelo 3 (todas las numéricas) ──\n")
print(summary(m3))

# VIF preliminar (diagnóstico de colinealidad, semana 10)
cat("\nVIF del Modelo 3:\n")
print(vif(m3))
# Regla orientadora: VIF_j > 10 → problema serio (semana 10, sección 8).

# ── 6. Modelo 4 – Con dummies de mes ─────────────────────────
# Justificación: si los boxplots del EDA muestran diferencias entre
# meses, las dummies de mes mejoran la especificación de E(Y|X)
# (semana 1, sección 8: ANOVA de un factor como regresión).
# Se usa mes como factor; R crea automáticamente p-1 = 11 indicadoras.
m4 <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + month,
         data = datos_limpios)
cat("\n── Modelo 4 (numéricas + mes) ──\n")
print(summary(m4))

# Prueba F parcial: ¿agregar month mejora significativamente?
cat("\nF parcial m3 vs m4 (¿aporta month?):\n")
print(anova(m3, m4))
# H0: todos los coeficientes de month = 0.
# Si p < 0.05 → el mes aporta información más allá de los índices FWI.

# ── 7. Modelo 5 – Selección automática (punto de partida) ────
# Justificación: con K=8 predictores numéricos hay 2^8 = 256 modelos
# posibles. Usamos selección paso a paso (semana 8, sección 17)
# para acotar el espacio. Los candidatos resultantes se analizan
# con criterios penalizados (R²_aj, Cp, AIC).
m_full <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + rain,
             data = datos_limpios)
m_step <- step(m_full, direction = "both",
               scope = list(lower = ~ 1, upper = m_full),
               trace = 1)
cat("\n── Modelo seleccionado por step() ──\n")
print(summary(m_step))

# ── 8. Comparación de modelos con criterios penalizados ──────
modelos <- list(m0 = m0, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m_step = m_step)
comp <- data.frame(
  Modelo  = names(modelos),
  p       = sapply(modelos, function(m) length(coef(m))),
  R2      = round(sapply(modelos, function(m) {
    s <- summary(m); ifelse(is.null(s$r.squared), 0, s$r.squared) }), 4),
  R2_adj  = round(sapply(modelos, function(m) {
    s <- summary(m); ifelse(is.null(s$adj.r.squared), 0, s$adj.r.squared) }), 4),
  AIC     = round(sapply(modelos, AIC), 2),
  SCE     = round(sapply(modelos, deviance), 4)
)
comp$Cp <- round(comp$SCE / (deviance(m3) / df.residual(m3)) +
                   2 * comp$p - nrow(datos_limpios), 2)
cat("\n── Tabla comparativa de modelos ──\n")
print(comp[order(comp$AIC), ])
# Esta tabla replica la estructura de la Tabla de la semana 8 (unidad quirúrgica).
# Se busca R²_aj alto + Cp ≈ p + AIC bajo → candidatos promisorios.

# ── 9. Plan de diagnóstico (informe final) ───────────────────
# Una vez elegido el modelo final, se ejecutarán los siguientes
# diagnósticos (semanas 7, 9, 10, 11):

# 9a. Residuos vs valores ajustados (semana 7: linealidad y homocedasticidad)
# plot(modelo_final, which = 1)

# 9b. QQ-plot de residuos estandarizados (semana 7: normalidad)
# plot(modelo_final, which = 2)

# 9c. Scale-Location (semana 9: heterocedasticidad)
# plot(modelo_final, which = 3)

# 9d. Residuos vs leverage / distancia de Cook (semana 11)
# plot(modelo_final, which = 5)

# 9e. VIF formal (semana 10)
# car::vif(modelo_final)

# 9f. Test de Breusch-Pagan (semana 9, sección 15)
# lmtest::bptest(modelo_final)

# 9g. DFBETAS y DFFITS (semana 11, secciones 5 y 8)
# influence.measures(modelo_final)

# ── FIN DEL PLAN DE MODELAMIENTO ─────────────────────────────