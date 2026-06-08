# ============================================================
# BLOQUE 6: EVALUACIÓN DE AJUSTE, INFERENCIA Y PREDICCIÓN
# ============================================================
# Cubrir puntos 6 y 7 de la rúbrica del informe:
#   - Coeficientes, IC, pruebas t y F
#   - R², R²aj, ANOVA de regresión
#   - Intervalos de confianza para la media E(Y|x₀)
#   - Intervalos de predicción para una nueva observación Y₀
#   - Evaluación de desempeño fuera de muestra (RMSE, MAE, cobertura IP)
# ============================================================

# Requiere: M_sel, M0, M_met, M_fwi, M_full, train, test


# ============================================================
# 6.1 TABLA DE INFERENCIA DEL MODELO SELECCIONADO
# ============================================================
cat("=== Resumen del modelo seleccionado ===\n")
sum_sel <- summary(M_sel)
print(sum_sel)

# Extraer tabla de coeficientes con IC 95%
cat("\n--- Coeficientes con IC 95% ---\n")
coef_ic <- cbind(
  coef(M_sel),
  confint(M_sel, level = 0.95)
)
colnames(coef_ic) <- c("Estimación", "IC 2.5%", "IC 97.5%")
print(round(coef_ic, 4))
# Interpretación de cada βⱼ:
# Un aumento de 1 unidad en xⱼ (manteniendo el resto fijo) está
# asociado a un cambio promedio de βⱼ en log(area+1).
# En escala original: multiplicar área por exp(βⱼ).


# ============================================================
# 6.2 ANOVA DE REGRESIÓN Y F GLOBAL
# ============================================================
cat("\n--- ANOVA de regresión ---\n")
print(anova(M_sel))
# Interpretación:
# SCT = SCR + SCE
# F global prueba H₀: β₁ = β₂ = ... = β_{p-1} = 0

cat("\n--- Métricas de ajuste ---\n")
cat("R²:          ", round(sum_sel$r.squared,     4), "\n")
cat("R² ajustado: ", round(sum_sel$adj.r.squared,  4), "\n")
cat("F global:    ", round(sum_sel$fstatistic[1],  3),
    "  gl:", sum_sel$fstatistic[2], ",", sum_sel$fstatistic[3],
    "  p-valor:", format.pval(pf(sum_sel$fstatistic[1],
                                  sum_sel$fstatistic[2],
                                  sum_sel$fstatistic[3],
                                  lower.tail = FALSE), digits = 4), "\n")
cat("S (σ̂):       ", round(sum_sel$sigma, 4), "(escala log(area+1))\n")


# ============================================================
# 6.3 INTERVALOS DE CONFIANZA PARA LA MEDIA CONDICIONAL
# ============================================================
# IC para E(Y|x₀) usando la fórmula del curso:
# μ̂₀ ± t_{1-α/2}(n-p) · S · √(x₀ᵀ(XᵀX)⁻¹x₀)
# En R: predict(..., interval = "confidence")

# Perfil de ejemplo: observación "típica" del parque
x0_tipico <- data.frame(
  FFMC    = median(train$FFMC),
  DMC     = median(train$DMC),
  DC      = median(train$DC),
  ISI     = median(train$ISI),
  temp    = median(train$temp),
  RH      = median(train$RH),
  wind    = median(train$wind),
  estacion = "verano"    # ajustar según fórmula de M_sel
)

cat("\n--- IC 95% para E(log(area+1)) en perfil típico ---\n")
ic_media <- predict(M_sel, newdata = x0_tipico,
                    interval = "confidence", level = 0.95)
print(ic_media)
cat("En escala original (área, ha): exp(IC) - 1 =",
    round(exp(ic_media) - 1, 3), "\n")


# ============================================================
# 6.4 INTERVALOS DE PREDICCIÓN PARA UNA NUEVA OBSERVACIÓN
# ============================================================
# IP para Y₀ = x₀ᵀβ + ε₀ usando:
# Ŷ₀ ± t_{1-α/2}(n-p) · S · √(1 + x₀ᵀ(XᵀX)⁻¹x₀)
# El "+1" bajo la raíz agrega la variabilidad del nuevo error ε₀
# (ver Clase Semana 6: IP siempre más ancho que IC para la media)

cat("\n--- IP 95% para nueva observación en perfil típico ---\n")
ip_nueva <- predict(M_sel, newdata = x0_tipico,
                    interval = "prediction", level = 0.95)
print(ip_nueva)
cat("En escala original (área, ha): exp(IP) - 1 =",
    round(exp(ip_nueva) - 1, 3), "\n")

cat("\nDiferencia entre IP e IC (ancho):\n")
cat("Ancho IC 95%: ", round(ic_media[3] - ic_media[2], 4), "\n")
cat("Ancho IP 95%: ", round(ip_nueva[3] - ip_nueva[2], 4), "\n")
# El IP es más ancho porque incorpora la variabilidad de la nueva obs.


# ============================================================
# 6.5 EVALUACIÓN PREDICTIVA EN TEST SET
# ============================================================
cat("\n=== Evaluación de desempeño fuera de muestra (test set) ===\n")

# Predicciones puntuales sobre test
pred_log  <- predict(M_sel, newdata = test)
real_log  <- test$log_area

# Métricas en escala log(area+1)
rmse_log  <- sqrt(mean((pred_log - real_log)^2))
mae_log   <- mean(abs(pred_log - real_log))
mape_log  <- mean(abs((pred_log - real_log) / (real_log + 1e-8))) * 100

cat("--- Escala log(area+1) ---\n")
cat("RMSE: ", round(rmse_log, 4), "\n")
cat("MAE:  ", round(mae_log,  4), "\n")
cat("MAPE: ", round(mape_log, 2), "%\n")

# Métricas en escala original (ha)
pred_orig <- exp(pred_log) - 1
real_orig <- test$area
rmse_orig <- sqrt(mean((pred_orig - real_orig)^2))
mae_orig  <- mean(abs(pred_orig - real_orig))

cat("\n--- Escala original (área, ha) ---\n")
cat("RMSE: ", round(rmse_orig, 2), "ha\n")
cat("MAE:  ", round(mae_orig,  2), "ha\n")


# ============================================================
# 6.6 COBERTURA DE INTERVALOS DE PREDICCIÓN
# ============================================================
# La teoría dice: con el modelo correcto y normalidad, el IP 95%
# debería contener el 95% de nuevas observaciones.
# Una cobertura empírica muy distinta de 0.95 sugiere problemas.

ip_test <- predict(M_sel, newdata = test,
                   interval = "prediction", level = 0.95)
dentro  <- (real_log >= ip_test[, "lwr"]) & (real_log <= ip_test[, "upr"])
cobertura <- mean(dentro)
cat("\n--- Cobertura empírica del IP 95% en test ---\n")
cat("Cobertura: ", round(cobertura, 4),
    "(esperado: 0.95)\n")
# Si cobertura << 0.95: los IP son demasiado angostos
#   → posible no-normalidad, heterocedasticidad o falta de ajuste
# Si cobertura >> 0.95: los IP son demasiado anchos (modelo conservador)


# ============================================================
# 6.7 GRÁFICOS DE PREDICCIÓN
# ============================================================
pdf("figuras/fig9_predicciones_test.pdf", width = 12, height = 5)
par(mfrow = c(1, 2))

# Panel 1: predicho vs real (escala log)
lim  <- range(c(real_log, pred_log, ip_test))
plot(real_log, pred_log,
     xlab = "log(area+1) observado (test)",
     ylab = "log(area+1) predicho",
     main = "Figura 9a: Predicciones vs. Valores reales",
     pch  = 19, cex = 0.6, col = "#00000077",
     xlim = lim, ylim = lim)
abline(0, 1, col = "red", lty = 2, lwd = 1.5)
legend("topleft",
       legend = c(paste0("RMSE = ", round(rmse_log, 3)),
                  paste0("MAE  = ", round(mae_log,  3))),
       bty = "n", cex = 0.85)

# Panel 2: con intervalos de predicción al 95%
plot(real_log, pred_log,
     xlab = "log(area+1) observado (test)",
     ylab = "Predicción e IP 95%",
     main = "Figura 9b: Intervalos de predicción 95%",
     pch  = 19, cex = 0.5, col = "#00000077",
     ylim = lim)
arrows(real_log, ip_test[,"lwr"], real_log, ip_test[,"upr"],
       code = 3, angle = 90, length = 0.02,
       col  = ifelse(dentro, "#4472C455", "#ED7D3155"))
points(real_log, pred_log, pch = 19, cex = 0.5, col = "#00000099")
abline(0, 1, col = "red", lty = 2, lwd = 1.5)
legend("topleft",
       legend = c(paste0("Cobertura = ", round(cobertura*100, 1), "%"),
                  "Dentro del IP", "Fuera del IP"),
       col = c("black", "#4472C4","#ED7D31"), pch = c(NA, 15, 15),
       bty = "n", cex = 0.85)
par(mfrow = c(1, 1))
dev.off()
cat("Figura 9 guardada en figuras/fig9_predicciones_test.pdf\n")


# ============================================================
# 6.8 TABLA FINAL COMPARATIVA (para el informe)
# ============================================================
modelos_eval <- list(M0=M0, M_met=M_met, M_fwi=M_fwi, M_sel=M_sel, M_full=M_full)
nombres      <- c("M0","M_met","M_fwi","M_sel","M_full")

rmse_train <- sapply(modelos_eval, function(m)
  sqrt(mean(resid(m)^2)))
rmse_test_v <- sapply(modelos_eval, function(m)
  sqrt(mean((test$log_area - predict(m, newdata=test))^2)))

tabla_final <- data.frame(
  Modelo   = nombres,
  Vars     = sapply(modelos_eval, function(m) length(coef(m))),
  R2       = sapply(modelos_eval, function(m) round(summary(m)$r.squared, 4)),
  R2aj     = sapply(modelos_eval, function(m) round(summary(m)$adj.r.squared, 4)),
  AIC      = sapply(modelos_eval, function(m) round(AIC(m), 1)),
  RMSE_tr  = round(rmse_train, 4),
  RMSE_test= round(rmse_test_v, 4),
  Brecha   = round(rmse_test_v - rmse_train, 4)  # indicador de sobreajuste
)

cat("\n=== TABLA FINAL: Comparación completa de modelos ===\n")
print(tabla_final, row.names = FALSE)
cat("\nNota: 'Brecha' = RMSE_test - RMSE_train. Valores altos indican sobreajuste.\n")
cat("El modelo seleccionado busca minimizar RMSE_test con la menor brecha posible.\n")
