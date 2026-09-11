# Dashboard SARIMA — Remesas de trabajadores en Colombia

Dashboard interactivo construido en R (Shiny) para analizar el comportamiento de las remesas de trabajadores hacia Colombia mediante modelos de series de tiempo, y su relación con la tasa de cambio (TRM).

🔗 **Demo en vivo:** https://mariana24i.shinyapps.io/DashboardRemesas/

## Datos

Serie mensual de remesas hacia Colombia, 2000–2020 (252 observaciones), y serie de TRM para el mismo periodo. Fuente: Banco de la República.

## Metodología

**Modelado de series de tiempo (metodología Box-Jenkins):**
- Transformación logarítmica y diferenciación de la serie para lograr estacionariedad.
- Pruebas de estacionariedad ADF (Augmented Dickey-Fuller) y KPSS sobre la serie original y diferenciada.
- Análisis de correlogramas ACF/PACF para identificar el orden del modelo.
- Estimación y comparación de 4 modelos SARIMA candidatos mediante AIC, BIC, AICc y el test de Ljung-Box.
- Modelo seleccionado: **SARIMA(1,1,2)(2,0,0)[12]**.
- Diagnóstico de supuestos del modelo ganador: normalidad de residuos (Jarque-Bera), ruido blanco (Ljung-Box), QQ-plot.
- Pronóstico a 12 meses (2021) con intervalos de confianza del 80% y 95%, comparado contra los datos reales y evaluado con MAE, RMSE y MAPE.

**Análisis remesas vs. TRM:**
- Correlación de Pearson calculada en niveles y en diferencias, para descartar una relación espuria.
- Test de cointegración de Engle-Granger (regresión de largo plazo + prueba ADF sobre los residuos).
- Test de causalidad de Granger en ambas direcciones.

## Resultados principales

- La serie original no es estacionaria (confirmado por ADF/KPSS); tras diferenciar en logaritmos, sí lo es.
- El pronóstico 2021 capturó la tendencia general, aunque subestimó algunos meses.
- La correlación remesas–TRM es fuerte en niveles (r ≈ 0,82) pero casi nula en diferencias (p ≈ 0,34), lo que sugiere una posible correlación espuria; se corroboró mediante el test de cointegración de Engle-Granger.
- Estacionalidad marcada en diciembre, coherente con un mayor envío de dinero en fin de año.

## Herramientas

R · Shiny (shinydashboard) · paquetes: `forecast`, `tseries`, `urca`, `lmtest`, `tidyverse`, `readxl`, `zoo` · visualización: `plotly`, `ggplot2` · tablas interactivas: `DT`.
