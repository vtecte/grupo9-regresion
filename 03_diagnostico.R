# ============================================================
# BLOQUE 4: DIAGNÓSTICO COMPLETO DEL MODELO SELECCIONADO
# ============================================================
# Cubrir los 4 supuestos diagnósticos de la rúbrica:
#   (a) Linealidad / especificación de la media
#   (b) Normalidad de residuos
#   (c) Homocedasticidad
#   (d) Colinealidad
# Y los diagnósticos adicionales:
#   (e) Autocorrelación (aunque datos no son series de tiempo)
#   (f) Observaciones influyentes (Cook, DFFITS, leverage)
# ============================================================
# Requiere: M_sel ajustado sobre train (Bloque 3)

n_tr <- nrow(train)
p_s  <- length(coef(M_sel))


# ============================================================
# (a) LINEALIDAD: partial residual plots (crPlots)
# ============================================================
# crPlots muestra, para cada predictor continuo xⱼ,
# el gráfico de eⱼ* = ej + β̂j·xij vs. xij.
# Una relación lineal indica que la media está bien especificada.
# Una curvatura sugiere que falta un término (ej. xⱼ²).
# REFERENCIA TEÓRICA: Clase Semana 7 (diagnóstico de media)
# A diferencia de los ggplots que usamos en el EDA donde se muestra la regresion
# simple de la respuesta con respecto a una sola variable, es decir, hay una
# relacion marginal (se ignora el efecto de otras variables)
# Los crplots tienen una relación condiciional o sea evaluamos la relación dado
# que las otras variables ya están en el modelo.

# Ejemplos para aclarar: Control de Confusión,"No controla nada. Si el viento
# influye en el área y está correlacionado con la temperatura, el gráfico de 
# ggplot estará ""contaminado"" por el efecto del viento.",
# en cambio crplots "Controla por todos los demás predictores. Muestra el efecto
# ""puro"" de la temperatura tras restar el efecto lineal del viento y la humedad."

cat("--- Gráfico de residuos componente + residuo (crPlots) ---\n")
pdf("figuras/fig3_crplots.pdf", width = 12, height = 8)
crPlots(M_sel,
        main = "Figura 3: Residuos componente + residuo (crPlots)",
        smooth = list(smoother = loessLine))
dev.off()

# Si crPlots de "temp" muestra curvatura clara:
# → Ajustar modelo con temp²: lm(... + temp + I(temp^2) ...)
# → Verificar con F-test si β₂_temp es significativo
# NOTA: I(temp^2) sigue siendo lineal en parámetros (A1 se cumple)

# PRUEBA FORMAL de no-linealidad en temp (Box-Tidwell):
# Si se sospecha xⱼ^λ con λ ≠ 1:
# boxTidwell(log_area ~ temp, other.x = ~ FFMC + DMC + ..., data = train)


# ============================================================
# (b) NORMALIDAD DE RESIDUOS
# ============================================================
# Se necesita para justificar inferencia exacta (pruebas t, F,
# IC e IP). Supuesto A5: ε|X ~ Nₙ(0, σ²Iₙ).

# QQ-plot de residuos estandarizados
pdf("figuras/fig4_qqplot_residuos.pdf", width = 6, height = 6)
r_std <- rstandard(M_sel)    # residuos estandarizados ri = ei / (S√(1-hii))
qqnorm(r_std,
       main  = "Figura 4: QQ-plot de residuos estandarizados",
       xlab  = "Cuantiles N(0,1) teóricos",
       ylab  = "Cuantiles muestrales",
       pch   = 19, cex = 0.6, col = "#00000077")
qqline(r_std, col = "red", lwd = 1.5)
dev.off()

# Test formal: Kolmogorov-Smirnov con corrección Lilliefors
# H₀: los residuos estandarizados siguen una N(0,1)
# Se usa Lilliefors porque media y varianza son estimadas
cat("\n--- Test de Lilliefors (normalidad de residuos estandarizados) ---\n")
test_norm <- lillie.test(r_std)
print(test_norm)
# p-valor < 0.05: evidencia contra normalidad
# Si se rechaza normalidad, los IC y IP exactos son aproximados.
# Acción: reportar con cautela, mencionar que los IC/IP son
# aproximados para muestras grandes (TCL justifica aproximación).


# ============================================================
# (c) HOMOCEDASTICIDAD
# ============================================================
# Supuesto A3: Cov(ε|X) = σ²Iₙ

# Gráfico: residuos vs ajustados (diagnóstico visual central)
pdf("figuras/fig5_residuos_ajustados.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Panel izquierdo: residuos ordinarios vs ajustados
plot(fitted(M_sel), resid(M_sel),
     xlab = "Valores ajustados Ŷ",
     ylab = "Residuos eᵢ",
     main = "Figura 5a: Residuos vs Ajustados",
     pch = 19, cex = 0.5, col = "#00000077")
abline(h = 0, col = "red", lty = 2)
lines(lowess(fitted(M_sel), resid(M_sel)), col = "#ED7D31", lwd = 2)

# Panel derecho: √|residuos estandarizados| vs ajustados
# (Scale-Location plot: mejor para detectar heterocedasticidad)
sqrt_std <- sqrt(abs(rstandard(M_sel)))
plot(fitted(M_sel), sqrt_std,
     xlab = "Valores ajustados Ŷ",
     ylab = "√|Residuos estand.|",
     main = "Figura 5b: Scale-Location",
     pch = 19, cex = 0.5, col = "#00000077")
lines(lowess(fitted(M_sel), sqrt_std), col = "#ED7D31", lwd = 2)
par(mfrow = c(1, 1))
dev.off()

# Test formal: Breusch-Pagan
# H₀: σᵢ² = σ² para todo i (homocedasticidad)
# Se regresa uᵢ = eᵢ²/σ̂² - 1 sobre las covariables del modelo
# LM = n·R²_aux ~ χ²_q bajo H₀
cat("\n--- Test de Breusch-Pagan (homocedasticidad) ---\n")
test_bp <- bptest(M_sel)
print(test_bp)
# Si p-valor < 0.05: heterocedasticidad detectada
# → Ver Bloque 4b para posibles remedios


# ============================================================
# (c) ESTADÍSTICO bᵢ (diagnóstico visual complementario)
# ============================================================
# Bajo homocedasticidad: E(bᵢ|X) = σ²
# bᵢ = eᵢ² / (1 - hᵢᵢ)  → corrige varianza por leverage
h_ii <- hatvalues(M_sel)
b_i  <- resid(M_sel)^2 / (1 - h_ii)

pdf("figuras/fig6_bi_ajustados.pdf", width = 6, height = 5)
plot(fitted(M_sel), b_i,
     xlab = "Valores ajustados Ŷ",
     ylab = expression(b[i] == e[i]^2 / (1-h[ii])),
     main = "Figura 6: Estadístico bᵢ (dispersión corregida por leverage)",
     pch  = 19, cex = 0.5, col = "#00000077")
abline(h = mean(b_i), col = "red", lty = 2)
lines(lowess(fitted(M_sel), b_i), col = "#4472C4", lwd = 2)
dev.off()
# Una tendencia creciente en bᵢ confirma heterocedasticidad


# ============================================================
# (d) COLINEALIDAD: VIF y número de condición
# ============================================================
cat("\n--- VIF del modelo seleccionado ---\n")
vif_sel <- vif(M_sel)
print(round(vif_sel, 3))
# Para variables continuas: VIF > 10 es problemático
# Para variables categóricas: GVIF^(1/(2·Df)) > √10 ≈ 3.16
cat("VIF > 10 en variables continuas:", 
    sum(vif_sel[, 1] > 10, na.rm = TRUE), "\n")

# Número de condición (diagnóstico global)
X_mat <- model.matrix(M_sel)
X_sin_inter <- X_mat[, -1]   # sin intercepto para centrar
svd_X <- svd(scale(X_sin_inter, center = TRUE, scale = FALSE))
kappa <- max(svd_X$d) / min(svd_X$d)
cat("Número de condición κ₂(X):", round(kappa, 2), "\n")
# κ₂ < 10: bien condicionado
# 10 < κ₂ < 30: colinealidad moderada
# κ₂ > 30: colinealidad severa


# ============================================================
# (e) AUTOCORRELACIÓN
# ============================================================
# Aunque los datos NO son una serie temporal clásica, los datos
# del parque están ordenados espacialmente y podrían tener
# dependencia en el orden de registro.
# Argumento: "dado que no existe una secuencia temporal explícita
# de observaciones, no se espera estructura fuerte. Sin embargo,
# el test se realiza como verificación."

cat("\n--- Test de Durbin-Watson (autocorrelación de primer rezago) ---\n")
test_dw <- dwtest(M_sel)
print(test_dw)
# D ≈ 2 → no hay autocorrelación
# D < 2 → posible autocorrelación positiva
# D > 2 → posible autocorrelación negativa
# Si se detecta: podría indicar observaciones del mismo sector
# del parque agrupadas en el dataset (estructura espacial)

# Función de autocorrelación (ACF) de residuos
pdf("figuras/fig7_acf_residuos.pdf", width = 7, height = 5)
acf(resid(M_sel),
    main  = "Figura 7: ACF de residuos del modelo seleccionado",
    xlab  = "Rezago",
    ylab  = "Autocorrelación",
    lag.max = 20)
dev.off()
# Barras dentro de las bandas azules punteadas → no hay autocorrelación detectada


# ============================================================
# (f) OBSERVACIONES INFLUYENTES
# ============================================================
# Tres preguntas distintas (Clase Semana 10):
#   ¿Es raro Yᵢ dado Xᵢ?  → |tᵢ| (outlier en respuesta)
#   ¿Es raro Xᵢ?          → hᵢᵢ (leverage)
#   ¿Cambia el ajuste?     → Dᵢ (distancia de Cook)
# La influencia requiere AMBAS: residuo relevante Y leverage relevante.

r_stud <- rstudent(M_sel)   # residuos estudentizados externos
cook_d <- cooks.distance(M_sel)
dff    <- dffits(M_sel)

# Umbrales descriptivos
umbral_h  <- 2 * p_s / n_tr    # leverage: 2p/n
umbral_t  <- 3                  # outlier: |tᵢ| > 3
umbral_dk <- 4 / n_tr           # Cook: 4/n (regla descriptiva)
umbral_df <- 2 * sqrt(p_s / n_tr)  # DFFITS: 2√(p/n)

cat("\n--- Diagnóstico de observaciones influyentes ---\n")
cat("Leverage alto (hᵢᵢ > 2p/n =", round(umbral_h, 4), "):",
    sum(h_ii > umbral_h), "\n")
cat("Outliers en respuesta (|tᵢ| > 3):",
    sum(abs(r_stud) > umbral_t), "\n")
cat("Cook > 4/n =", round(umbral_dk, 4), ":",
    sum(cook_d > umbral_dk), "\n")
cat("|DFFITS| > 2√(p/n) =", round(umbral_df, 4), ":",
    sum(abs(dff) > umbral_df), "\n")

# Top 5 más influyentes
cat("\nTop 5 por distancia de Cook:\n")
top5 <- order(cook_d, decreasing = TRUE)[1:5]
print(data.frame(
  obs_train = top5,
  cook   = round(cook_d[top5], 5),
  h_ii   = round(h_ii[top5], 4),
  t_i    = round(r_stud[top5], 3),
  dffits = round(dff[top5], 3),
  log_area_real = round(train$log_area[top5], 3)
))

# Gráfico resumen de influencia
pdf("figuras/fig8_influencia.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

# Leverage
plot(h_ii, type = "h",
     ylab = expression(h[ii]),
     main = "Figura 8a: Leverage por observación",
     col  = ifelse(h_ii > umbral_h, "#ED7D31", "#4472C4"))
abline(h = umbral_h, col = "red", lty = 2)
legend("topright", c(paste0("h > ", round(umbral_h,3)), "h normal"),
       col = c("#ED7D31","#4472C4"), pch = 15, bty = "n", cex = 0.8)

# Residuos estudentizados externos
plot(r_stud, type = "h",
     ylab = expression(t[i]),
     main = "Figura 8b: Residuos estudentizados externos",
     col  = ifelse(abs(r_stud) > umbral_t, "#ED7D31", "#4472C4"))
abline(h = c(-umbral_t, umbral_t), col = "red", lty = 2)

# Distancia de Cook
plot(cook_d, type = "h",
     ylab = expression(D[i]),
     main = "Figura 8c: Distancia de Cook",
     col  = ifelse(cook_d > umbral_dk, "#ED7D31", "#4472C4"))
abline(h = umbral_dk, col = "red", lty = 2)
# Identificar los puntos más influyentes
text(top5, cook_d[top5], labels = top5, pos = 3, cex = 0.7, col = "darkred")

# Mapa de diagnóstico: leverage vs |residuo|
plot(h_ii, abs(r_stud),
     xlab = expression(h[ii]),
     ylab = expression("|"*t[i]*"|"),
     main = "Figura 8d: Leverage vs. |Residuo estudentizado|",
     pch  = 19, cex = 0.5, col = "#00000055")
abline(v = umbral_h, h = umbral_t, col = "red", lty = 2)
# Cuadrante superior derecho = casos potencialmente influyentes
text(top5, abs(r_stud[top5]), labels = top5, pos = 3, cex = 0.7, col = "darkred")
par(mfrow = c(1, 1))
dev.off()

cat("Figuras de influencia guardadas en figuras/fig8_influencia.pdf\n")


# ============================================================
# TABLA RESUMEN DE SUPUESTOS (para el informe, Tabla X)
# ============================================================
cat("\n=== RESUMEN DE DIAGNÓSTICOS DEL MODELO SELECCIONADO ===\n")
cat("Supuesto  | Test              | Estadístico | p-valor | Conclusión\n")
cat("----------+-------------------+-------------+---------+--------------------\n")
cat(sprintf("A2 (exog.) | crPlots / residuos| visual      | -       | Ver Fig. 3 y 5\n"))
cat(sprintf("A3 (homoc.)| Breusch-Pagan     | LM=%.4f   | %.4f  | %s\n",
            test_bp$statistic, test_bp$p.value,
            ifelse(test_bp$p.value > 0.05, "No se rechaza H₀", "Se rechaza H₀")))
cat(sprintf("A5 (normal)| Lilliefors        | D=%.4f    | %.4f  | %s\n",
            test_norm$statistic, test_norm$p.value,
            ifelse(test_norm$p.value > 0.05, "No se rechaza H₀", "Se rechaza H₀")))
cat(sprintf("(autocor.) | Durbin-Watson     | DW=%.4f   | %.4f  | %s\n",
            test_dw$statistic, test_dw$p.value,
            ifelse(test_dw$p.value > 0.05, "No se rechaza H₀", "D≠2, revisar ACF")))
