# grupo9-regresion

## Este repositorio contiene:

1. /data: Carpeta que almacena base de datos fireforest.cvs.
2. /figuras: Carpeta que almacena los gráficos obtenidos en los códigos.
3. 00r_setup_y_eda.R: Preparación de los datos, transformaciones iniciales y EDA.
4. 01r_preparacion_datos.R: División 80/20 (train - test) de los datos, agrupación en Estaciones y ultimas desiciones de descarte de coeficientes.
5. 02r_modelos_y_seleccion.R: Planteamiento de los modelos, Test F para modelos anidados.
6. 03_diagnostico.R: Diagnóstico del modelo elegido; Comprobación de supuestos.
7. 04_remedios_y_sensibilidad.R: Análisis de sensibilidad, ante problematicas encontradas en código anterior.
8. 05_evaluacion_prediccion.R: Prueba de modelo elegido en subconjunto de datos test.
9. grupo9_regresion.Rproj: Enviroment para trabajar en R
10. informe_final.zip: Informe con plantilla. Para visualizar subir como proyecto en Overleaf. 

Pueden revisar la investigación original [A Data Mining Approach to Predict Forest Fires using Meteorological Data](
https://fileserver-az.core.ac.uk/download/pdf/55609027.pdf)

Para una reproducción optima del contenido del repositorio, abrir `grupo9_regresion.Rproj` y luego los códigos de 00 a 05 en orden. Si no se pueden visualizar las gráficas en R, correr las líneas de código del gráfico despúes de la función pdf() y antes de dev.off(), ejemplo,

```R
1 | pdf("figuras/fig1_distribucion_respuesta.pdf", width = 9, height = 4)
2 | par(mfrow = c(1, 2))
3 | hist(datos$area,
4 |      main = "Figura 1a: Distribución de area",
5 |      xlab = "Área quemada (ha)",
6 |      ylab = "Frecuencia",
7 |      col  = "#4472C4", border = "white", breaks = 40)
8 | text(500, 200,
9 |      paste0("Sesgo = ", round(skewness(datos$area), 2),
10|             "\nCurt. = ", round(kurtosis(datos$area), 2)),
11|      cex = 0.85)
12|
13| hist(datos$log_area,
14|      main = "Figura 1b: Distribución de log(area+1)",
15|      xlab = "log(área + 1)",
16|      ylab = "Frecuencia",
17|      col  = "#ED7D31", border = "white", breaks = 30)
18| text(3, 80,
19|      paste0("Sesgo = ", round(skewness(datos$log_area), 2),
20|             "\nCurt. = ", round(kurtosis(datos$log_area), 2)),
21|      cex = 0.85)
22| par(mfrow = c(1, 1))
23| dev.off()
24| cat("Figura 1 guardada en figuras/fig1_distribucion_respuesta.pdf\n")
```
Se debería correr el codigo entre la línea 2 y la línea 22. En caso de querer ver solo un gráfico por ejemplo la distribución de área, correr el código entre la linea 3 y la línea 8. 
