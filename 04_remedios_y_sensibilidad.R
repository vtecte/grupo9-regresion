# ============================================================
# BLOQUE 5: REMEDIOS Y ANÁLISIS DE SENSIBILIDAD
# ============================================================
# La rúbrica exige: "propone remedios razonables" y
# "análisis de sensibilidad: qué cambia ante decisiones
# razonables de limpieza, transformación o retiro de casos".
# ============================================================

# Requiere: M_sel, train, test, top5 (del Bloque 4)


# ============================================================
# 5.1 REMEDIO ANTE HETEROCEDASTICIDAD: WLS o transformación
# ============================================================
# Si el test de Breusch-Pagan rechazó homocedasticidad:
#
# OPCIÓN A: WLS (mínimos cuadrados ponderados)
# Idea: Var(εᵢ|X) = σ²/wᵢ → observaciones más variables
# reciben menor peso. WLS = MCO en modelo transformado por W^(1/2).
# Problema práctico: los pesos exactos son desconocidos.
# Estimación de pesos: wᵢ ∝ 1/Var(εᵢ), que se aproxima con
# el cuadrado de los ajustados o de algún predictor sospechoso.

# Paso 1: estimar la función de varianza
log_e2   <- log(resid(M_sel)^2 + 1e-6)  # log de residuos cuadrados
fit_var  <- lm(log_e2 ~ fitted(M_sel))   # regresión auxiliar de log(e²) ~ Ŷ
cat("Regresión auxiliar de log(e²) ~ Ŷ para estimar pesos WLS:\n")
summary(fit_var)
# Si la pendiente es significativa → varianza crece con Ŷ

# Paso 2: pesos estimados
w_hat <- 1 / exp(fitted(fit_var))        # wᵢ = 1/Var_estimada

# Paso 3: ajustar WLS
M_wls <- lm(log_area ~ . - area - month - day - rain - estacion,   
            # sustituir ~ formula(M_sel) para usar exactamente las mismas variables:
            # M_wls <- lm(formula(M_sel), data = train, weights = w_hat)
            data = train, weights = w_hat)
# CORRECTO: usar la misma fórmula que M_sel
M_wls <- lm(formula(M_sel), data = train, weights = w_hat)

cat("\n--- Comparación M_sel vs M_wls ---\n")
cat("M_sel: R²aj =", round(summary(M_sel)$adj.r.squared, 4),
    "| AIC =", round(AIC(M_sel), 2), "\n")
cat("M_wls: R²aj =", round(summary(M_wls)$adj.r.squared, 4),
    "| AIC =", round(AIC(M_wls), 2), "\n")

# Verificar si Breusch-Pagan mejora con WLS
cat("\nBreusch-Pagan sobre M_wls:\n")
print(bptest(M_wls))


# OPCIÓN B: Transformación de la respuesta (Box-Cox)
# Nota: Ya aplicamos log(area+1). Si la heterocedasticidad
# persiste, podría indicar que la escala transformada aún no
# es completamente adecuada, o que hay heterocedasticidad
# por algún predictor específico (no resolvible con transformación
# de Y).
# En cualquier caso, reportar el diagnóstico honestamente.


# ============================================================
# 5.2 ANÁLISIS DE SENSIBILIDAD POR OBSERVACIONES INFLUYENTES
# ============================================================
# Pregunta: ¿las conclusiones principales cambian si retiramos
# los casos más influyentes?
# NOTA: esto NO implica eliminarlos del modelo final.
# Solo se reporta la sensibilidad para ser honestos.

cat("\n--- Análisis de sensibilidad: M_sel con y sin top-5 influyentes ---\n")

# Identificar índices originales (en train) de los casos top5
# top5 ya definido en el Bloque 4

# Reajustar sin los 5 casos más influyentes
M_sin5 <- update(M_sel, data = train[-top5, ])

# Comparar coeficientes
cat("\nComparación de coeficientes:\n")
comp_coef <- cbind(
  "M_sel"   = round(coef(M_sel), 4),
  "M_sin5"  = round(coef(M_sin5)[names(coef(M_sel))], 4),
  "Cambio%" = round(100*(coef(M_sin5)[names(coef(M_sel))] - coef(M_sel)) /
                    (abs(coef(M_sel)) + 1e-8), 2)
)
print(comp_coef)

# DFBETAS: qué coeficiente cambia más al retirar cada obs.
cat("\n--- DFBETAS (top observaciones influyentes) ---\n")
dfb    <- dfbetas(M_sel)
umbral_dfb <- 2 / sqrt(n_tr)
cat("Umbral DFBETAS: 2/√n =", round(umbral_dfb, 4), "\n")

# Mostrar solo observaciones que superan el umbral en algún coeficiente
obs_dfb_alto <- which(apply(abs(dfb), 1, max) > umbral_dfb)
cat("Observaciones con |DFBETAS| > 2/√n en algún coeficiente:",
    length(obs_dfb_alto), "\n")
if (length(obs_dfb_alto) > 0) {
  print(round(dfb[obs_dfb_alto, , drop = FALSE], 3))
}
# Interpretación: si un coeficiente clave (ej. β_temp) tiene
# DFBETAS grande para un solo caso, la conclusión sobre ese
# predictor es sensible a ese dato.


# ============================================================
# 5.3 VERIFICAR TRANSFORMACIÓN DE temp (curvatura)
# ============================================================
# Si crPlots mostró curvatura en temp, comparar modelo con y sin temp²

M_sel_quad <- update(M_sel, . ~ . + I(temp^2))
cat("\n--- F-test: M_sel vs M_sel + temp² ---\n")
print(anova(M_sel, M_sel_quad))
# Si p < 0.05: el término cuadrático aporta significativamente
# → Reportar que la media puede tener curvatura en temperatura
# → Incluir I(temp^2) en el modelo final si es el caso

cat("Comparación AIC:\n")
cat("M_sel:          AIC =", round(AIC(M_sel), 2), "\n")
cat("M_sel + temp²:  AIC =", round(AIC(M_sel_quad), 2), "\n")


# ============================================================
# 5.4 TABLA DE SENSIBILIDAD PARA EL INFORME
# ============================================================
# Esta tabla responde al punto 8 de la rúbrica:
# "Análisis de sensibilidad"
cat("\n=== TABLA RESUMEN DE SENSIBILIDAD ===\n")
cat("Decisión                  | R²aj   | RMSE_test | Conclusión\n")
cat("--------------------------+--------+-----------+-----------\n")
cat(sprintf("Modelo seleccionado       | %.4f | %.4f    | Modelo base\n",
            summary(M_sel)$adj.r.squared,
            sqrt(mean((test$log_area - predict(M_sel, newdata=test))^2))))
cat(sprintf("Sin top-5 influyentes     | %.4f | %.4f    | Ver nota\n",
            summary(M_sin5)$adj.r.squared,
            sqrt(mean((test$log_area - predict(M_sin5, newdata=test))^2))))
cat(sprintf("Con WLS                   | %.4f | %.4f    | Ver nota\n",
            summary(M_wls)$adj.r.squared,
            sqrt(mean((test$log_area - predict(M_wls, newdata=test))^2))))
cat(sprintf("Con temp²                 | %.4f | %.4f    | Ver nota\n",
            summary(M_sel_quad)$adj.r.squared,
            sqrt(mean((test$log_area - predict(M_sel_quad, newdata=test))^2))))
cat("\nNota: si las conclusiones principales (qué variables predicen,\n")
cat("si el modelo es globalmente significativo) no cambian bajo estas\n")
cat("variaciones, el análisis es robusto.\n")
