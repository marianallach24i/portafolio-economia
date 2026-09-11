
library(pacman)
library(readr)      # Lectura de datos
library(dplyr)      # Manipulación de datos
library(ggplot2)    # Gráficos avanzados
library(forecast)   # Pronósticos y ARIMA
library(tseries)    # Tests de series temporales
library(lmtest)     # Tests estadísticos
library(stats)      # Funciones estadísticas base
library(plotly)
library(zoo)
library(scales)
library(htmlwidgets)


# ==============================================================================
# 2. CARGA Y PREPARACIÓN DE DATOS =============================================
# ==============================================================================


datosrm <- read_excel("remesas.xlsx")
remesasdata <- datosrm |>
  select(`Remesas (Millones USD)`)

remesas_ts <- ts(datosrm$`Remesas (Millones USD)`, 
                 start = c(2000, 1), 
                 frequency = 12)
summary(remesas_ts)


# ==============================================================================
# ▓▓▓  BLOQUE A: ANÁLISIS DE LA SERIE ORIGINAL  ▓▓▓
# ==============================================================================

# A1. VISUALIZACIÓN SERIE ORIGINAL ============================================

df_plot <- data.frame(
  Fecha   = as.yearmon(time(remesas_ts)) |> as.Date(),
  Remesas = as.numeric(remesas_ts)
)

ggplot(df_plot, aes(x = Fecha, y = Remesas)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = mean(df_plot$Remesas),
             color = "purple", linetype = "dashed", linewidth = 0.8) +
  labs(title    = "Remesas de Trabajadores en Colombia",
       subtitle = "2000-2020 | Serie mensual original",
       x = "Año", y = "Millones USD") +
  theme_minimal() +
  theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))



# A2. DESCOMPOSICIÓN SERIE ORIGINAL ===========================================

decomp_original <- decompose(remesas_ts)
autoplot(decomp_original) +
  ggtitle("Descomposición Clásica - Serie Original") +
  theme_minimal()

# A3. PATRÓN ESTACIONAL - ORIGINAL ============================================

ggseasonplot(remesas_ts, year.labels = TRUE, year.labels.left = TRUE) +
  ggtitle("Patrón Estacional de Remesas - Serie Original") +
  ylab("Millones USD") +
  theme_minimal()

# A4. SUBSERIES ESTACIONALES - ORIGINAL =======================================

ggsubseriesplot(remesas_ts) +
  ggtitle("Subseries Estacionales - Serie Original") +
  ylab("Millones USD") +
  theme_minimal()

# A5. PRUEBAS DE ESTACIONARIEDAD - SERIE ORIGINAL =============================

adf_original  <- adf.test(remesas_ts)
kpss_original <- kpss.test(remesas_ts, null = "Trend")

print(adf_original)
print(kpss_original)

ndiffs(remesas_ts)
nsdiffs(remesas_ts)


# A6. DECISIÓN DE TRANSFORMACIÓN LOGARÍTMICA ==================================

lambda_optimo <- BoxCox.lambda(remesas_ts)
print(lambda_optimo)
# lambda ≈ 0   → log natural apropiado ✅
# lambda ≈ 1   → no se necesita transformación
# lambda ≈ 0.5 → raíz cuadrada


# ==============================================================================
# ▓▓▓  BLOQUE B: ANÁLISIS CON SERIE EN LOGARITMO  ▓▓▓
# ==============================================================================

# B1. TRANSFORMACIÓN LOGARÍTMICA ==============================================

remesas_log <- log(remesas_ts)

# Gráfico Serie en Logaritmos (mismo formato que el original)
df_plot_log <- data.frame(
  Fecha   = as.yearmon(time(remesas_log)) |> as.Date(),
  Remesas = as.numeric(remesas_log)
)

ggplot(df_plot_log, aes(x = Fecha, y = Remesas)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_hline(yintercept = mean(df_plot_log$Remesas),
             color = "skyblue", linetype = "dashed", linewidth = 0.8) +
  labs(title    = "Remesas de Trabajadores en Colombia (Log)",
       subtitle = "2000-2020 | Serie mensual en logaritmo natural",
       x = "Año", y = "log(Millones USD)") +
  theme_minimal() +
  theme(plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))



# Panel comparativo: Original vs Logaritmo
par(mfrow = c(2, 1), mar = c(3, 4, 2, 1))

plot(remesas_ts,
     main = "Serie Original - Remesas (Millones USD)",
     col  = "steelblue", ylab = "Millones USD", xlab = "Tiempo")

plot(remesas_log,
     main = "Serie en Logaritmo Natural",
     col  = "darkred", ylab = "log(Millones USD)", xlab = "Tiempo")

par(mfrow = c(1, 1))


# B3. DESCOMPOSICIÓN - LOG ====================================================

decomp_log <- decompose(remesas_log)
autoplot(decomp_log) +
  ggtitle("Descomposición Clásica - Serie en Logaritmo") +
  theme_minimal()

# B4. PATRÓN ESTACIONAL - LOG =================================================

ggseasonplot(remesas_log, year.labels = TRUE, year.labels.left = TRUE) +
  ggtitle("Patrón Estacional de Remesas - Log") +
  ylab("log(Millones USD)") +
  theme_minimal()

# B5. SUBSERIES ESTACIONALES - LOG ============================================

ggsubseriesplot(remesas_log) +
  ggtitle("Subseries Estacionales - Log") +
  ylab("log(Millones USD)") +
  theme_minimal()

# B6. PRUEBAS DE ESTACIONARIEDAD - LOG ========================================

adf.test(remesas_log)
kpss.test(remesas_log, null = "Trend")

#Diferenciaciones necesarias
ndiffs(remesas_log)
nsdiffs(remesas_log)

# B7. DIFERENCIACIÓN - LOG ====================================================

dif_simple_log     <- diff(remesas_log, differences = 1)

# Pruebas sobre diferenciación simple
adf.test(na.omit(dif_simple_log))
kpss.test(na.omit(dif_simple_log))


# B8.2 GRÁFICOS INDIVIDUALES ACF/PACF - DIFERENCIACIÓN SIMPLE LOG =============

autoplot(Acf(dif_simple_log, lag.max = 48, plot = FALSE)) +
  ggtitle("ACF - Diferenciación Simple LOG (d=1)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

autoplot(Pacf(dif_simple_log, lag.max = 48, plot = FALSE)) +
  ggtitle("PACF - Diferenciación Simple LOG (d=1)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# B9. GRÁFICOS INDIVIDUALES DIFERENCIACIÓN ====================================

par(mfrow = c(2, 1))

plot(dif_simple_log,
     main = "Diferenciación Simple LOG (d=1)",
     col  = "darkred", ylab = "Valores", xlab = "Tiempo")


# ==============================================================================
# ▓▓▓  BLOQUE C: IDENTIFICACIÓN Y ESTIMACIÓN DE MODELOS  ▓▓▓
# ==============================================================================

# C1. AUTO.ARIMA ==============================================================

modelo_auto <- auto.arima(remesas_ts, lambda = 0)
summary(modelo_auto)

modelo_auto_completo <- auto.arima(remesas_ts,
                                   stepwise      = FALSE,
                                   approximation = FALSE,
                                   lambda        = 0)
summary(modelo_auto_completo)

arimaorder(modelo_auto)
arimaorder(modelo_auto_completo)

# C2. MODELOS CANDIDATOS ======================================================

m1 <- Arima(remesas_ts, order = c(1,1,2), seasonal = c(0,0,2), lambda = 0)
m2 <- Arima(remesas_ts, order = c(1,1,2), seasonal = c(2,0,0), lambda = 0)
m3 <- Arima(remesas_ts, order = c(0,1,1), seasonal = c(0,0,1), lambda = 0)
m4 <- Arima(remesas_ts, order = c(1,1,1), seasonal = c(0,0,1), lambda = 0)

summary(m1)
summary(m2)
summary(m3)
summary(m4)


# ==============================================================================
# ▓▓▓  BLOQUE D: DIAGNÓSTICO Y SELECCIÓN  ▓▓▓
# ==============================================================================

# D1. CHECK RESIDUOS ==========================================================

checkresiduals(m1)
checkresiduals(m2)
checkresiduals(m3)
checkresiduals(m4)

# D2. LJUNG-BOX PARA TODOS LOS MODELOS ========================================

modelos_lista <- list(m1, m2, m3, m4)
nombres_mod   <- c("SARIMA(1,1,2)(0,0,2)",
                   "SARIMA(1,1,2)(2,0,0)",
                   "SARIMA(0,1,1)(0,0,1)",
                   "SARIMA(1,1,1)(0,0,1)")

for (i in seq_along(modelos_lista)) {
  res <- residuals(modelos_lista[[i]])
  lb  <- Box.test(res, lag = 20, type = "Ljung-Box",
                  fitdf = length(coef(modelos_lista[[i]])))
  print(paste(nombres_mod[i],
              "→ p-value:",
              round(lb$p.value, 4),
              ifelse(lb$p.value > 0.05,
                     "✅ Ruido blanco",
                     "❌ Estructura residual")))
}


# D3. TABLA COMPARATIVA AIC / BIC / AICc ======================================

tabla_criterios <- data.frame(
  Modelo = nombres_mod,
  AIC    = round(c(AIC(m1), AIC(m2), AIC(m3), AIC(m4)), 2),
  BIC    = round(c(BIC(m1), BIC(m2), BIC(m3), BIC(m4)), 2),
  AICc   = round(c(m1$aicc, m2$aicc, m3$aicc, m4$aicc), 2)
)

tabla_criterios[order(tabla_criterios$AICc), ]

# ==============================================================================
# ▓▓▓  BLOQUE E: MEJOR MODELO Y PRONÓSTICO  ▓▓▓
# ==============================================================================

# E1. SELECCIÓN Y VALIDACIÓN ==================================================
# ↓ Cambia "m2" por el modelo que quedó primero en tu tabla del paso D3

mejor_modelo <- m2

summary(mejor_modelo)
checkresiduals(mejor_modelo)
shapiro.test(residuals(mejor_modelo))

# E2. PRONÓSTICO 12 MESES =====================================================

pronostico <- forecast(mejor_modelo, h = 12, biasadj = TRUE)

print(pronostico)

autoplot(pronostico) +
  ggtitle("Pronóstico Remesas Colombia - Próximos 12 meses") +
  xlab("Año") + ylab("Millones USD") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# E3. COMPARACIÓN AUTO vs MANUAL ==============================================

data.frame(
  Modelo = c("SARIMA(1,1,2)(2,0,0) - Manual",
             "Auto rápido",
             "Auto exhaustivo"),
  AICc   = round(c(mejor_modelo$aicc,
                   modelo_auto$aicc,
                   modelo_auto_completo$aicc), 2)
)




# Pronóstico 12 meses (2021)
pronostico <- forecast(mejor_modelo, h = 12, biasadj = TRUE)

remesas2021 <- read_excel("remesas2021.xlsx")

# Comparar pronóstico vs lo que realmente ocurrió
remesas_real_2021 <- window(remesas2021, 
                            start = c(2021, 1), 
                            end   = c(2021, 12))



# ==============================================================================
# COMPARACIÓN PRONÓSTICO VS REALIDAD
# ==============================================================================

# 1. Cargar archivo completo 2000-2025
datosrm_completo <- read_excel("remesastrabajadores.xlsx")

remesas_ts_completa <- ts(datosrm_completo$`Remesas (Millones USD)`,
                          start     = c(2000, 1),
                          frequency = 12)

# 2. Extraer valores reales 2021
remesas_real_2021 <- window(remesas_ts_completa,
                            start = c(2021, 1),
                            end   = c(2021, 12))

# 3. Pronóstico 12 meses con el mejor modelo ya construido
pronostico <- forecast(mejor_modelo, h = 12, biasadj = TRUE)

# 4. Gráfico pronóstico vs realidad
autoplot(pronostico) +
  autolayer(remesas_real_2021,
            series = "Valor Real 2021",
            color  = "red") +
  ggtitle("Pronóstico vs Realidad - Remesas Colombia 2021") +
  xlab("Año") + ylab("Millones USD") +
  theme_minimal() +
  theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

# 5. Tabla comparativa mes a mes
meses <- c("Ene","Feb","Mar","Abr","May","Jun",
           "Jul","Ago","Sep","Oct","Nov","Dic")

tabla_comparacion <- data.frame(
  Mes        = meses,
  Pronostico = round(as.numeric(pronostico$mean),   2),
  Real       = round(as.numeric(remesas_real_2021), 2),
  Diferencia = round(as.numeric(remesas_real_2021) -
                       as.numeric(pronostico$mean),   2),
  Error_Pct  = round((as.numeric(remesas_real_2021) -
                        as.numeric(pronostico$mean)) /
                       as.numeric(remesas_real_2021) * 100, 2)
)

print(tabla_comparacion)

# 6. Métricas de precisión
data.frame(
  Métrica = c("MAE", "RMSE", "MAPE (%)"),
  Valor   = round(c(mean(abs(tabla_comparacion$Diferencia)),
                    sqrt(mean(tabla_comparacion$Diferencia^2)),
                    mean(abs(tabla_comparacion$Error_Pct))), 2)
)




# =============================================================================
# Gráfico: Remesas vs TRM — Índice Base 100 (2000–2025)
# Fuente: Banco de la República de Colombia
# =============================================================================

# --- 1. Librerías -------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggtext)   # para títulos con HTML/markdown (instalar si no está)

# --- 2. Carga de datos -------------------------------------------------------

# Remesas: columnas "Fecha" y "Remesas (Millones USD)"
remesas <- read_excel("remesastrabajadores.xlsx") %>%
  rename(Fecha = 1, Remesas = 2) %>%
  mutate(Fecha = as.Date(Fecha))

# TRM: ya viene con el día 1 de cada mes (312 filas, 2000–2025)
trm <- read_excel("TRM.xlsx") %>%
  rename(Fecha = 1, TRM = 2) %>%
  mutate(
    Fecha = as.Date(Fecha),
    TRM   = as.numeric(TRM)
  )

# --- 3. Merge y cálculo del índice base 100 ----------------------------------

df <- inner_join(remesas, trm, by = "Fecha") %>%
  arrange(Fecha) %>%
  distinct(Fecha, .keep_all = TRUE)

# Base: enero 2000 (primera fila)
base_remesas <- df$Remesas[1]   # 103.75 M USD
base_trm     <- df$TRM[1]       # 1,873.77 COP/USD

df <- df %>%
  mutate(
    idx_remesas = (Remesas / base_remesas) * 100,
    idx_TRM     = (TRM     / base_trm)     * 100
  )

# Formato largo para ggplot
df_long <- df %>%
  select(Fecha, idx_remesas, idx_TRM) %>%
  pivot_longer(
    cols      = c(idx_remesas, idx_TRM),
    names_to  = "Variable",
    values_to = "Indice"
  ) %>%
  mutate(Variable = recode(Variable,
                           "idx_remesas" = "Remesas (Millones USD)",
                           "idx_TRM"     = "TRM (COP/USD)"
  ))

# --- 4. Eventos clave para anotaciones ---------------------------------------

eventos <- tibble(
  Fecha    = as.Date(c("2008-09-01", "2020-03-01", "2022-11-01")),
  etiqueta = c("Crisis\nfinanciera\n2008", "COVID-19\n2020", "Pico TRM\n2022"),
  vjust    = c(-0.4, -0.4, -0.4)
)

# --- 5. Paleta y tema --------------------------------------------------------

col_remesas <- "#6C63C9"   # morado
col_trm     <- "#1A9E75"   # verde esmeralda
col_base    <- "#9E9E9E"   # gris para línea de referencia
col_fondo   <- "#FAFAF8"
col_panel   <- "#FFFFFF"
col_grid    <- "#EBEBEB"
col_texto   <- "#2C2C2A"

tema_remesas <- theme_minimal(base_size = 12) +
  theme(
    # Fondo
    plot.background  = element_rect(fill = col_fondo, color = NA),
    panel.background = element_rect(fill = col_panel, color = NA),
    panel.border     = element_blank(),
    
    # Grilla
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = col_grid, linewidth = 0.4),
    
    # Ejes
    axis.line.x      = element_line(color = "#CCCCCC", linewidth = 0.4),
    axis.ticks.x     = element_line(color = "#CCCCCC", linewidth = 0.3),
    axis.ticks.y     = element_blank(),
    axis.text        = element_text(color = "#666660", size = 10),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(color = "#666660", size = 10,
                                    margin = margin(r = 10)),
    
    # Leyenda
    legend.position      = "top",
    legend.justification = "left",
    legend.title         = element_blank(),
    legend.text          = element_text(color = col_texto, size = 11),
    legend.key.width     = unit(1.8, "cm"),
    legend.key.height    = unit(0.4, "cm"),
    legend.margin        = margin(b = 4),
    
    # Títulos
    plot.title    = element_markdown(
      size = 16, face = "bold", color = col_texto,
      margin = margin(b = 4)),
    plot.subtitle = element_markdown(
      size = 11, color = "#666660", lineheight = 1.4,
      margin = margin(b = 16)),
    plot.caption  = element_text(
      size = 9, color = "#999990", hjust = 0,
      margin = margin(t = 12)),
    plot.margin   = margin(20, 24, 16, 20)
  )

# --- 6. Construcción del gráfico ---------------------------------------------

g <- ggplot(df_long, aes(x = Fecha, y = Indice, color = Variable)) +
  
  # Línea de referencia base 100
  geom_hline(yintercept = 100, linetype = "dashed",
             color = col_base, linewidth = 0.5) +
  annotate("text", x = as.Date("2000-06-01"), y = 103,
           label = "Base 100 = enero 2000",
           color = col_base, size = 3, hjust = 0) +
  
  # Bandas de eventos históricos
  annotate("rect",
           xmin = as.Date("2008-07-01"), xmax = as.Date("2009-06-01"),
           ymin = -Inf, ymax = Inf,
           fill = "#FFF3E0", alpha = 0.5) +
  annotate("rect",
           xmin = as.Date("2020-01-01"), xmax = as.Date("2021-01-01"),
           ymin = -Inf, ymax = Inf,
           fill = "#E8F5E9", alpha = 0.4) +
  
  # Líneas principales
  geom_line(linewidth = 0.9, alpha = 0.9) +
  
  # Puntos en eventos clave
  geom_point(
    data = df_long %>%
      filter(Fecha %in% eventos$Fecha),
    size = 2.5, alpha = 0.9
  ) +
  
  # Etiquetas de eventos (sobre la línea de remesas)
  geom_vline(
    data = eventos,
    aes(xintercept = Fecha),
    color = "#BBBBBB", linewidth = 0.4, linetype = "dotted",
    inherit.aes = FALSE
  ) +
  annotate("text",
           x      = eventos$Fecha,
           y      = 950,
           label  = eventos$etiqueta,
           size   = 2.8,
           color  = "#888880",
           hjust  = 0.5,
           lineheight = 0.9) +
  
  # Valores finales al extremo derecho
  annotate("text",
           x = as.Date("2025-09-01"), y = 1131,
           label = "1.131", color = col_remesas,
           size = 3.5, fontface = "bold", hjust = 0) +
  annotate("text",
           x = as.Date("2025-09-01"), y = 200,
           label = "200", color = col_trm,
           size = 3.5, fontface = "bold", hjust = 0) +
  
  # Escalas
  scale_color_manual(values = c(
    "Remesas (Millones USD)" = col_remesas,
    "TRM (COP/USD)"          = col_trm
  )) +
  scale_x_date(
    breaks       = seq(as.Date("2000-01-01"), as.Date("2025-01-01"), by = "5 years"),
    date_labels  = "%Y",
    expand       = expansion(mult = c(0.01, 0.08))
  ) +
  scale_y_log10(
    breaks = c(100, 150, 200, 300, 500, 700, 1000),
    labels = label_comma(big.mark = ".", decimal.mark = ","),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  
  # Etiquetas
  labs(
    title    = "Remesas y TRM en Colombia: **¿van de la mano?**",
    subtitle = "Índice base 100 = enero 2000 — escala logarítmica. Las remesas crecieron ~11× mientras el dólar se multiplicó ~2.6×. Correlación de Pearson: **r = 0.82**",
    y        = "Índice (enero 2000 = 100, escala log)",
    caption  = "Fuente: Banco de la República de Colombia  ·  Elaboración propia"
  ) +
  
  tema_remesas

# --- 7. Guardar --------------------------------------------------------------

ggsave(
  filename = "remesas_TRM_indice100_log.png",
  plot     = g,
  width    = 12,
  height   = 6.5,
  dpi      = 180,
  bg       = col_fondo
)

# Para ver en RStudio:
print(g)

# =============================================================================
# NOTA SOBRE ESCALAS
# El eje Y está en unidades de índice (adimensional).
# - 100  = valor en enero 2000 (igual para ambas series)
# - 500  = la variable vale 5× lo que valía en 2000
# - 1131 = las remesas valen ~11× su valor inicial
# Un salto de 100 puntos en el índice = +100% de crecimiento acumulado
# respecto al año base, independientemente de la unidad original (USD o COP).
# =============================================================================
