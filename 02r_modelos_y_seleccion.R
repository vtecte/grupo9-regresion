# ============================================================
# BLOQUE 3: MODELOS CANDIDATOS Y SELECCIÓN BACKWARDS
# ============================================================
# ARGUMENTO DEL PROFESOR (dos cosas que hay que defender):
#
#  (A) "Por qué puedo hacer este modelo":
#      → Verificar supuestos A1–A4 (y A5 si se usa inferencia exacta):
#        A1: linealidad en parámetros (sí, por construcción)
#        A2: exogeneidad E(ε|X)=0 (gráficos de residuos)
#        A3: homocedasticidad (Breusch-Pagan)
#        A4: rk(X) = p (verificar, VIF)
#        A5: normalidad condicional (Lilliefors)
#        Además: n/p >> 1 (relación obs/parámetros razonable)
#
#  (B) "Por qué otros modelos son peores":
#      → F-tests parciales (modelos anidados)
#      → AIC, R²aj, RMSE en test set (modelos no anidados)
#      → Parsimonia: R²aj puede bajar si se agregan predictores irrelevantes
#
# ESTRATEGIA: 3+1 modelos con propósito claro
#   M0:   modelo nulo (solo intercepto) → baseline
#   M_met: predictores meteorológicos directos → hipótesis 1
#   M_fwi: índices FWI → hipótesis 2
#   M_full: todos los candidatos → punto de partida backwards
#   M_sel: resultado del backwards sobre M_full en TRAIN → modelo final
# ============================================================

# Requiere: train y test cargados desde 01_preparacion_datos.R

# NOTA: El modelo lineal que se asume en todos los casos es:
#   Y = Xβ + ε,  ε|X ~ N(0, σ²Iₙ)
# donde Y = log(area+1), Xβ = E(Y|X) (lineal en β)
#Recordemos incluir la interpretacion del coef de mallows Cp y la fórmula, tambié
# del AIC que incluye un factor que penaliza la complejidad del modelo (agregar 
# variables)

# ============================================================
# 3.1 Modelo nulo: M0
# ============================================================
# Justificación: sirve únicamente como baseline. Cualquier modelo
# con predictores debe superar a M0 para tener valor añadido.
# Esto lo hace La prueba F global implícitamente.

M0 <- lm(log_area ~ 1, data = train)
cat("M0 (nulo): p =", length(coef(M0)), "| AIC:", round(AIC(M0), 2), "\n")


# ============================================================
# 3.2 Modelo meteorológico: M_met
# ============================================================
# Justificación: temp, RH y wind son variables observadas
# directamente Representan
# las condiciones ambientales en el momento del incendio.
# Excluimos rain: 97.6% de observaciones tienen rain=0
# Excluimos day: no aporta significativamente
# Especificación del modelo:
#   log(area+1) = β₀ + β₁·temp + β₂·RH + β₃·wind + ε

M_met <- lm(log_area ~ temp + RH + wind, data = train)
cat("\nM_met (meteorológico): p =", length(coef(M_met)),
    "| AIC:", round(AIC(M_met), 2),
    "| R²aj:", round(summary(M_met)$adj.r.squared, 4), "\n")
summary(M_met)


# ============================================================
# 3.3 Modelo FWI: M_fwi
# ============================================================
# Justificación: FFMC, DMC, DC e ISI son los componentes del
# sistema FWI canadiense, diseñados específicamente para capturar
# peligro de incendio. Resumen temperatura, humedad y viento en
# índices interpretables de riesgo. Si estos índices están bien
# calibrados, deberían capturar todo lo relevante de las variables
# meteorológicas brutas (en el informe pongamos lo que sale en la página de que 
# incluyen variables meteorológicas pero no sabemos como (son cajas negras).

# Especificación del modelo:
#   log(area+1) = β₀ + β₁·FFMC + β₂·DMC + β₃·DC + β₄·ISI + ε

M_fwi <- lm(log_area ~ FFMC + DMC + DC + ISI, data = train)
cat("\nM_fwi (índices FWI): p =", length(coef(M_fwi)),
    "| AIC:", round(AIC(M_fwi), 2),
    "| R²aj:", round(summary(M_fwi)$adj.r.squared, 4), "\n")
summary(M_fwi)


# ============================================================
# 3.4 Modelo completo: M_full
# ============================================================
# Justificación: punto de partida para la selección backwards.
# Incluye todos los candidatos numéricos  + estación
# Se excluyen rain (cuasi-constante) y day (sin aporte).
# Especificación:
#   log(area+1) = β₀ + β₁·FFMC + β₂·DMC + β₃·DC + β₄·ISI
#               + β₅·temp + β₆·RH + β₇·wind
#               + Σ βₖ·I(estacion=k) + ε

M_full <- lm(log_area ~ FFMC + DMC + DC + ISI + temp + RH + wind + estacion,
             data = train)
cat("\nM_full (completo): p =", length(coef(M_full)),
    "| AIC:", round(AIC(M_full), 2),
    "| R²aj:", round(summary(M_full)$adj.r.squared, 4), "\n")
cat("Dimensiones de X:", nrow(model.matrix(M_full)), "x", ncol(model.matrix(M_full)), "\n")
cat("Ratio n/p:", round(nrow(train)/length(coef(M_full)), 1),
    "(mínimo recomendado ≈ 10)\n") #confirmar que el minimo recomendado es ese
summary(M_full)

# Verificar colinealidad del modelo completo ANTES de selección
cat("\nVIF del modelo completo:\n")
print(round(vif(M_full), 2))

# GVIF^(1/(2·Df)) < √10 ≈ 3.16 para variables categóricas (quizás esto haya que 
# SACARLO porque no me acuerdo de haberlo visto)

# esta table tiene que ir SI o SI (por confirmar el GVIF^(1/(2*Df)))


# ============================================================
# 3.5 Selección backwards: M_sel
# ============================================================
# JUSTIFICACIÓN de usar backwards (y no forward ni stepwise):
#   - La pregunta es "QUÉ variables predicen el área": partimos
#     del modelo más rico razonable y eliminamos lo innecesario.
#   - Forward puede omitir variables que solo son importantes
#     conjuntamente (problema de poder en escenarios correlacionados).
#   - Backwards parte de la solución más completa y elimina lo que
#     no aporta marginalmente, dado el resto.

#
# CRITERIO: AIC (equivalente a Mallows Cp cuando σ² se estima).
#   Cp de Mallows ≈ AIC en escala de complejidad vs. ajuste.
#   Se elige el modelo donde Cp ≈ p 
#
# PROTECCIÓN CONTRA OVERFITTING (FEEDBACK ATENDIDO):
#   step() se aplica ÚNICAMENTE a los datos de entrenamiento.
#   La fórmula resultante se fija y se evalúa sobre TEST sin
#   re-ejecutar el proceso de selección. El test set permanece
#   completamente "ciego" a la selección.

cat("\n--- Selección backwards con AIC (solo sobre TRAIN) ---\n")
M_sel <- step(M_full,
              direction = "backward",
              trace     = 1)        # trace=1 muestra pasos; usar 0 para silencio

cat("\n--- Modelo seleccionado ---\n")
print(formula(M_sel))
cat("p seleccionado:", length(coef(M_sel)), "\n")
cat("AIC:", round(AIC(M_sel), 2), "\n")
cat("R²aj:", round(summary(M_sel)$adj.r.squared, 4), "\n")
summary(M_sel)

#parece que elimina el invierno preguntar a matías y revisar el código


# ============================================================
# 3.6 F-tests entre modelos anidados
# ============================================================
# Para justificar "por qué otros modelos son peores"
# (argumentos formales con distribución F bajo H₀)

cat("\n=== F-test: M0 vs M_sel (¿algún predictor aporta?) ===\n")
print(anova(M0, M_sel))
# H₀: todos los βj del modelo seleccionado = 0
# Si se rechaza → el modelo seleccionado aporta significativamente ( lo cual es 
# nuestro caso valor-p muy chico)

# Si M_met está anidado en M_sel:
# (depende de qué variables sobrevivieron la selección)
# Ajustar según resultado de step()

cat("\n=== F-test: M_met vs M_sel (¿FWI + estacion aportan sobre solo clima?) ===\n")
# SOLO si M_met ⊆ M_sel (modelos anidados)
# Verificar con: all(names(coef(M_met)) %in% names(coef(M_sel)))
if (all(attr(terms(M_met), "term.labels") %in% attr(terms(M_sel), "term.labels"))) {
  print(anova(M_met, M_sel))
} else {
  cat("M_met no está anidado en M_sel. Comparar con AIC/RMSE-test.\n")
}

# no se debe hacer F test porque no está anidado

cat("\n=== F-test: M_fwi vs M_sel (¿variables meteorológicas + estacion aportan sobre solo FWI?) ===\n")
if (all(attr(terms(M_fwi), "term.labels") %in% attr(terms(M_sel), "term.labels"))) {
  print(anova(M_fwi, M_sel))
} else {
  cat("M_fwi no está anidado en M_sel. Comparar con AIC/RMSE-test.\n")
}

cat("\n=== F-test: M_sel vs M_full (¿vale la complejidad extra del completo?) ===\n")
print(anova(M_sel, M_full))
# NO se rechaza H₀ → M_sel es suficiente (preferimos parsimonia).
# Este argumento justifica M_sel sobre M_full.


# ============================================================
# 3.7 Tabla comparativa de modelos (para el informe)
# ============================================================
calcular_rmse <- function(modelo, datos_nuevos) {
  pred <- predict(modelo, newdata = datos_nuevos)
  sqrt(mean((datos_nuevos$log_area - pred)^2))
}

modelos     <- list(M0 = M0, M_met = M_met, M_fwi = M_fwi,
                    M_sel = M_sel, M_full = M_full)
nom_modelos <- c("M0 (nulo)", "M_met (clima)", "M_fwi (FWI)",
                 "M_sel (seleccionado)", "M_full (completo)")

tabla_comp <- data.frame(
  Modelo    = nom_modelos,
  p         = sapply(modelos, function(m) length(coef(m))),
  SCE_train = sapply(modelos, function(m) round(sum(resid(m)^2), 3)),
  R2aj      = sapply(modelos, function(m) round(summary(m)$adj.r.squared, 4)),
  AIC       = sapply(modelos, function(m) round(AIC(m), 2)),
  RMSE_train= sapply(modelos, function(m) round(calcular_rmse(m, train), 4)),
  RMSE_test = sapply(modelos, function(m) round(calcular_rmse(m, test), 4))
)

cat("\n=== Tabla comparativa de modelos ===\n")
print(tabla_comp, row.names = FALSE)
# LEER: Un modelo sobreajustado tendrá RMSE_train << RMSE_test.
# El modelo seleccionado debería mostrar la mejor relación
# RMSE_test vs. parsimonia (p).
