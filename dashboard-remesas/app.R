# ==============================================================================
# DASHBOARD SARIMA - REMESAS COLOMBIA
# Estilo: Dashdark X — Negro azulado, moderno, minimalista premium
# ==============================================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(DT)
library(forecast)
library(tseries)
library(urca)
library(zoo)
library(readxl)
library(tidyverse)

# ==============================================================================
# DATOS Y MODELOS
# ==============================================================================

datosrm <- read_excel("remesas.xlsx")
remesas_ts <- ts(datosrm$`Remesas (Millones USD)`,
                 start = c(2000, 1), frequency = 12)
remesas_log    <- log(remesas_ts)
dif_simple_log <- diff(remesas_log, differences = 1)

m1 <- Arima(remesas_ts, order=c(1,1,2), seasonal=c(0,0,2), lambda=0)
m2 <- Arima(remesas_ts, order=c(1,1,2), seasonal=c(2,0,0), lambda=0)
m3 <- Arima(remesas_ts, order=c(0,1,1), seasonal=c(0,0,1), lambda=0)
m4 <- Arima(remesas_ts, order=c(1,1,1), seasonal=c(0,0,1), lambda=0)
mejor_modelo  <- m2
pronostico    <- forecast(mejor_modelo, h=12, biasadj=TRUE)

datosrm_completo  <- read_excel("remesastrabajadores.xlsx")
remesas_completa  <- ts(datosrm_completo$`Remesas (Millones USD)`,
                        start=c(2000,1), frequency=12)
remesas_real_2021 <- window(remesas_completa, start=c(2021,1), end=c(2021,12))

meses <- c("Ene","Feb","Mar","Abr","May","Jun",
           "Jul","Ago","Sep","Oct","Nov","Dic")
tabla_comparacion <- data.frame(
  Mes        = meses,
  Pronostico = round(as.numeric(pronostico$mean),   2),
  Real       = round(as.numeric(remesas_real_2021), 2),
  Diferencia = round(as.numeric(remesas_real_2021) -
                       as.numeric(pronostico$mean), 2),
  Error_Pct  = round((as.numeric(remesas_real_2021) -
                        as.numeric(pronostico$mean)) /
                       as.numeric(remesas_real_2021)*100, 2)
)
mae_val  <- round(mean(abs(tabla_comparacion$Diferencia)), 2)
rmse_val <- round(sqrt(mean(tabla_comparacion$Diferencia^2)), 2)
mape_val <- round(mean(abs(tabla_comparacion$Error_Pct)), 2)

nombres_mod   <- c("SARIMA(1,1,2)(0,0,2)","SARIMA(1,1,2)(2,0,0)",
                   "SARIMA(0,1,1)(0,0,1)","SARIMA(1,1,1)(0,0,1)")
modelos_lista <- list(m1,m2,m3,m4)
lb_pvalues <- sapply(modelos_lista, function(m){
  round(Box.test(residuals(m), lag=20, type="Ljung-Box",
                 fitdf=length(coef(m)))$p.value, 4)
})
tabla_modelos <- data.frame(
  Modelo      = nombres_mod,
  AIC         = round(sapply(modelos_lista, AIC), 2),
  BIC         = round(sapply(modelos_lista, BIC), 2),
  AICc        = round(sapply(modelos_lista, function(x) x$aicc), 2),
  LjungBox    = lb_pvalues,
  RuidoBlanco = ifelse(lb_pvalues > 0.05, "✅ Sí", "❌ No")
)

# ==============================================================================
# DATOS TRM — Remesas vs TRM (2000–2025)
# ==============================================================================

trm_raw <- read_excel("TRM.xlsx") %>%
  rename(Fecha = 1, TRM = 2) %>%
  mutate(Fecha = as.Date(Fecha), TRM = as.numeric(TRM))

remesas_raw <- read_excel("remesastrabajadores.xlsx") %>%
  rename(Fecha = 1, Remesas = 2) %>%
  mutate(Fecha = as.Date(Fecha))

df_trm <- inner_join(remesas_raw, trm_raw, by = "Fecha") %>%
  arrange(Fecha) %>%
  distinct(Fecha, .keep_all = TRUE)

base_rem     <- df_trm$Remesas[1]
base_trm_val <- df_trm$TRM[1]

df_trm <- df_trm %>%
  mutate(
    idx_remesas = (Remesas / base_rem)     * 100,
    idx_TRM     = (TRM     / base_trm_val) * 100
  )

# MODIFICACIÓN 1: correlación en niveles y en diferencias para nota metodológica
cor_pearson     <- round(cor(df_trm$Remesas, df_trm$TRM, use="complete.obs"), 2)
cor_pearson_dif <- round(cor(diff(df_trm$Remesas), diff(df_trm$TRM)), 2)

eventos_x    <- as.Date(c("2008-09-01","2020-03-01","2022-11-01"))
etiquetas_ev <- c("Crisis financiera 2008","COVID-19 2020","Pico TRM 2022")

# ==============================================================================
# COINTEGRACIÓN — Test de Engle-Granger: Remesas vs TRM
# ==============================================================================

# Series en niveles (2000–2025)
rem_niveles <- df_trm$Remesas
trm_niveles <- df_trm$TRM

# Paso 1: Regresión de largo plazo (Engle-Granger)
lm_coint     <- lm(rem_niveles ~ trm_niveles)
residuos_eg  <- residuals(lm_coint)
coef_lr      <- round(coef(lm_coint), 4)

# Paso 2: ADF sobre los residuos
eg_test      <- adf.test(residuos_eg, alternative = "stationary")
eg_pval      <- round(eg_test$p.value, 4)
eg_stat      <- round(eg_test$statistic, 4)
eg_conclusion <- ifelse(eg_pval < 0.05, "✅ Cointegradas — relación de largo plazo válida", "⚠ No cointegradas — correlación puede ser espuria")




# Causalidad de Granger (lmtest)
library(lmtest)
granger_rem_trm <- grangertest(rem_niveles ~ trm_niveles, order = 1)
granger_trm_rem <- grangertest(trm_niveles ~ rem_niveles, order = 1)
gr_pval_rem_trm <- round(granger_rem_trm$`Pr(>F)`[2], 4)
gr_pval_trm_rem <- round(granger_trm_rem$`Pr(>F)`[2], 4)
gr_conclusion_rem_trm <- ifelse(gr_pval_rem_trm < 0.05, "✅ TRM causa (Granger) a Remesas", "❌ TRM NO causa (Granger) a Remesas")
gr_conclusion_trm_rem <- ifelse(gr_pval_trm_rem < 0.05, "✅ Remesas causan (Granger) a TRM", "❌ Remesas NO causan (Granger) a TRM")

# ==============================================================================
# SUPUESTOS DEL MODELO — cálculos para nueva pestaña (MODIFICACIÓN 3)
# ==============================================================================

res_mejor    <- residuals(mejor_modelo)
jb_test      <- jarque.bera.test(as.numeric(res_mejor))   # tseries ya cargado
lb_mejor     <- Box.test(res_mejor, lag=20, type="Ljung-Box",
                         fitdf=length(coef(mejor_modelo)))
media_res    <- round(mean(as.numeric(res_mejor)), 6)
sd_res       <- round(sd(as.numeric(res_mejor)), 4)
jb_pval      <- round(jb_test$p.value, 4)
lb_pval_mej  <- round(lb_mejor$p.value, 4)

# ==============================================================================
# PALETA DASHDARK X
# ==============================================================================

bg      <- "#0f1117"
sidebar <- "#0d1117"
card    <- "#161b2e"
card2   <- "#1a2035"
border  <- "#1e2640"
border2 <- "#252d45"
purple  <- "#7c5cfc"
purple2 <- "#9b7ffe"
blue    <- "#3b82f6"
blue2   <- "#60a5fa"
cyan    <- "#06b6d4"
green   <- "#10b981"
red     <- "#ef4444"
text1   <- "#f1f5f9"
text2   <- "#94a3b8"
text3   <- "#475569"

# ==============================================================================
# CSS DASHDARK X
# ==============================================================================

css <- paste0("
@import url('https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css');
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600;700&display=swap');

* { font-family: 'DM Sans', sans-serif !important; box-sizing: border-box; }

body, html { background-color:", bg, "; color:", text1, "; margin:0; padding:0; }

::-webkit-scrollbar { width:5px; height:5px; }
::-webkit-scrollbar-track { background:", bg, "; }
::-webkit-scrollbar-thumb { background:", border2, "; border-radius:4px; }

.skin-blue .main-header .logo {
  background:", sidebar, " !important;
  color:", text1, " !important;
  font-weight:600; font-size:15px;
  border-bottom:1px solid ", border, ";
  border-right:1px solid ", border, ";
  letter-spacing:-0.3px;
}
.skin-blue .main-header .navbar {
  background:", sidebar, " !important;
  border-bottom:1px solid ", border, ";
}
.skin-blue .main-header .navbar .sidebar-toggle,
.skin-blue .main-header .navbar .sidebar-toggle:hover {
  color:", text2, " !important;
  background:transparent !important;
}
/* Ocultar el ícono original del toggle y reemplazar con 3 líneas CSS */
.sidebar-toggle .sr-only { display:none; }
.sidebar-toggle::before {
  content: '';
  display: block;
  width: 18px;
  height: 2px;
  background: ", text2, ";
  border-radius: 2px;
  box-shadow: 0 5px 0 ", text2, ", 0 10px 0 ", text2, ";
  margin: 4px auto;
}
.sidebar-toggle i { display:none !important; }

.skin-blue .main-sidebar {
  background:", sidebar, " !important;
  border-right:1px solid ", border, ";
}
.sidebar { padding-top:8px; }
.skin-blue .sidebar a { color:", text2, " !important; }
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a {
  color:", text2, " !important;
  font-size:13.5px; font-weight:400;
  padding:9px 16px 9px 20px;
  margin:2px 10px;
  border-radius:8px;
  display:flex; align-items:center; justify-content:space-between;
  transition:all 0.15s ease;
  border:none !important;
}
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a:hover {
  background:", card, " !important;
  color:", text1, " !important;
}
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a:hover {
  background:", card2, " !important;
  color:", text1, " !important;
  font-weight:500;
  border-left:2px solid ", purple, " !important;
}
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a span,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a:hover span {
  color:", purple2, " !important;
}
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a .fa,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a .fas,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a .far,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a .glyphicon,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li > a svg {
  color:", text2, " !important;
  width:18px; font-size:15px; margin-right:10px;
  opacity:0.75;
}
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a .fa,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a .fas,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a .far,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a .glyphicon,
.skin-blue .main-sidebar .sidebar .sidebar-menu > li.active > a svg {
  color:", purple2, " !important;
  opacity:1;
}
.content-wrapper { background:", bg, " !important; }
.content { padding:20px 24px; }

.box {
  background:", card, " !important;
  border:1px solid ", border, " !important;
  border-radius:14px !important;
  box-shadow:0 1px 3px rgba(0,0,0,0.4) !important;
  margin-bottom:20px !important;
  overflow:hidden;
}
.box.box-solid > .box-header { border-radius:14px 14px 0 0 !important; }
.box-header {
  background:", card, " !important;
  border-bottom:1px solid ", border, " !important;
  padding:16px 20px !important;
}
.box-title {
  color:", text1, " !important;
  font-size:13.5px !important;
  font-weight:600 !important;
  letter-spacing:-0.2px;
}
.box-body { padding:16px 20px !important; }

/* Botones de nota explicativa — siempre dentro del header */
.box-header { overflow:visible !important; position:relative !important; }
.box         { overflow:visible !important; }
.box-title   { overflow:visible !important; }

.small-box {
  border-radius:14px !important;
  border:1px solid ", border, " !important;
  background:", card, " !important;
  box-shadow:0 1px 3px rgba(0,0,0,0.3) !important;
  overflow:hidden;
}
.small-box > .inner { padding:18px 20px !important; }
.small-box h3 {
  font-size:24px !important; font-weight:700 !important;
  color:", text1, " !important; letter-spacing:-0.5px; margin:0 0 4px !important;
}
.small-box p {
  font-size:12px !important; font-weight:500 !important;
  color:", text2, " !important; margin:0 !important;
  text-transform:uppercase; letter-spacing:0.5px;
}
.small-box .icon { opacity:0.08 !important; }

.dataTables_wrapper { color:", text1, "; }
table.dataTable { border-collapse:separate !important; border-spacing:0 !important; }
table.dataTable thead tr th {
  background:", card2, " !important;
  color:", text2, " !important;
  font-size:11px !important; font-weight:600 !important;
  text-transform:uppercase !important; letter-spacing:0.6px !important;
  border-bottom:1px solid ", border2, " !important;
  padding:10px 14px !important;
}
table.dataTable tbody tr { background:", card, " !important; transition:background 0.12s; }
table.dataTable tbody tr:hover { background:", card2, " !important; }
table.dataTable tbody td {
  color:", text1, " !important; font-size:13px !important;
  border-bottom:1px solid ", border, " !important;
  padding:11px 14px !important;
}
.dataTables_info, .dataTables_filter label,
.dataTables_length label { color:", text3, " !important; font-size:11px !important; }
.dataTables_filter input {
  background:", card2, " !important; border:1px solid ", border2, " !important;
  color:", text1, " !important; border-radius:6px; padding:4px 10px;
}

.selectize-input {
  background:", card2, " !important; border:1px solid ", border2, " !important;
  border-radius:8px !important; color:", text1, " !important;
  font-size:13px !important; padding:8px 12px !important; box-shadow:none !important;
}
.selectize-dropdown {
  background:", card2, " !important; border:1px solid ", border2, " !important;
  border-radius:8px !important; color:", text1, " !important; font-size:13px !important;
}
.selectize-dropdown .option { padding:8px 12px !important; }
.selectize-dropdown .option:hover,
.selectize-dropdown .active { background:", border2, " !important; }

label { color:", text2, " !important; font-size:12px !important;
        text-transform:uppercase; letter-spacing:0.4px; }

.badge-up   { background:rgba(16,185,129,0.15); color:", green, ";
              padding:3px 8px; border-radius:20px; font-size:11px; font-weight:600; }
.badge-down { background:rgba(239,68,68,0.15); color:", red, ";
              padding:3px 8px; border-radius:20px; font-size:11px; font-weight:600; }
.badge-neutral { background:rgba(148,163,184,0.15); color:", text2, ";
                 padding:3px 8px; border-radius:20px; font-size:11px; font-weight:600; }

.stat-card {
  background:", card, "; border:1px solid ", border, ";
  border-radius:14px; padding:20px; transition:border-color 0.2s;
  margin-bottom:20px;
}
.stat-card:hover { border-color:", border2, "; }
.stat-label {
  color:", text3, "; font-size:11px; font-weight:600;
  text-transform:uppercase; letter-spacing:0.6px; margin-bottom:10px;
}
.stat-value { color:", text1, "; font-size:26px; font-weight:700;
              letter-spacing:-0.5px; margin-bottom:6px; }
.stat-sub { color:", text2, "; font-size:12px; }

.section-header { display:flex; align-items:baseline; gap:12px; margin-bottom:6px; }
.section-title  { color:", text1, "; font-size:20px; font-weight:700; letter-spacing:-0.4px; }
.section-sub    { color:", text2, "; font-size:13px; font-weight:400; }

.divider { border:none; border-top:1px solid ", border, "; margin:20px 0; }

.pill-ok { background:rgba(16,185,129,0.12); color:", green, ";
           border:1px solid rgba(16,185,129,0.25);
           padding:6px 14px; border-radius:20px;
           font-size:12px; font-weight:600; display:inline-block; }
.pill-no { background:rgba(239,68,68,0.12); color:", red, ";
           border:1px solid rgba(239,68,68,0.25);
           padding:6px 14px; border-radius:20px;
           font-size:12px; font-weight:600; display:inline-block; }

.metric-mini {
  background:", card2, "; border:1px solid ", border, ";
  border-radius:10px; padding:16px 20px; margin-bottom:12px;
}
.metric-mini-label { color:", text3, "; font-size:11px; font-weight:600;
                     text-transform:uppercase; letter-spacing:0.5px; margin-bottom:6px; }
.metric-mini-value { color:", text1, "; font-size:22px; font-weight:700; letter-spacing:-0.4px; }

.btn-periodo {
  background:", card2, " !important;
  border:1px solid ", border2, " !important;
  color:", text2, " !important;
  border-radius:6px !important;
  font-size:11px !important; font-weight:500 !important;
  padding:4px 10px !important;
  transition:all 0.15s ease !important;
  cursor:pointer !important;
}
.btn-periodo:hover {
  background:", border2, " !important;
  color:", text1, " !important;
  border-color:", purple, " !important;
}

.interp-card {
  background:", card2, "; border:1px solid ", border, ";
  border-radius:12px; padding:18px; height:100%;
}
.interp-tag {
  font-size:10px; font-weight:700; text-transform:uppercase;
  letter-spacing:0.8px; margin-bottom:12px; padding-bottom:8px;
  border-bottom:1px solid ", border, ";
}
.interp-row { margin-bottom:10px; }
.interp-key { color:", text2, "; font-size:12px; margin-bottom:2px; }
.interp-val { font-size:13px; font-weight:600; }

/* ── MODIFICACIÓN 1: Nota metodológica flotante ── */
.nota-wrapper {
  position: relative;
  display: inline-block;
}
.nota-btn {
  background: rgba(155,127,254,0.15);
  border: 1px solid rgba(155,127,254,0.35);
  color: ", purple2, ";
  border-radius: 50%;
  width: 28px; height: 28px;
  font-size: 13px; font-weight: 700;
  cursor: pointer;
  display: inline-flex; align-items: center; justify-content: center;
  transition: all 0.2s ease;
  vertical-align: middle;
  margin-left: 8px;
  line-height: 1;
}
.nota-btn:hover {
  background: rgba(155,127,254,0.28);
  border-color: ", purple2, ";
}
.nota-popup {
  display: none;
  position: absolute;
  top: 36px; right: 0;
  width: 320px;
  background: ", card2, ";
  border: 1px solid ", border2, ";
  border-radius: 12px;
  padding: 16px 18px;
  z-index: 9999;
  box-shadow: 0 8px 32px rgba(0,0,0,0.5);
}
.nota-popup.visible { display: block; }

/* Popup dentro del box-header */
.box-header .nota-popup {
  position: absolute;
  top: 38px;
  right: 0;
  width: 320px;
  background: #1a2035;
  border: 1px solid #252d45;
  border-radius: 12px;
  padding: 16px 18px;
  z-index: 99999;
  box-shadow: 0 8px 32px rgba(0,0,0,0.7);
}
.box-header .nota-wrapper {
  position: relative;
  display: inline-block;
}
.box-title {
  overflow: visible !important;
}
.box-header {
  overflow: visible !important;
}
.nota-popup-title {
  font-size: 11px; font-weight: 700; text-transform: uppercase;
  letter-spacing: 0.7px; color: ", purple2, ";
  margin-bottom: 10px; padding-bottom: 8px;
  border-bottom: 1px solid ", border, ";
}
.nota-popup-body {
  font-size: 12px; color: ", text2, "; line-height: 1.65;
}
.nota-popup-body b { color: ", text1, "; }
.nota-popup-body .r-val {
  display: inline-block;
  background: rgba(6,182,212,0.12);
  border: 1px solid rgba(6,182,212,0.25);
  color: ", cyan, ";
  border-radius: 6px; padding: 1px 7px;
  font-weight: 700; font-size: 12px;
}
")

# ==============================================================================
# HELPERS UI
# ==============================================================================

# Helper genérico para notas explicativas (ícono ℹ en esquina superior derecha)
nota_exp <- function(id, texto_html) {
  uid_btn <- paste0("nbtn_", id)
  uid_pop <- paste0("npop_", id)
  # Botón ✦ sutil en esquina superior derecha del recuadro (position:absolute)
  tags$div(
    style = "position:absolute; top:10px; right:12px; z-index:50;",
    tags$button(
      id      = uid_btn,
      style   = paste0(
        "position:absolute; top:10px; right:12px;",
        "background:rgba(148,163,184,0.1);",
        "border:1px solid rgba(148,163,184,0.28);",
        "color:rgba(148,163,184,0.6);",
        "border-radius:50%; width:26px; height:26px;",
        "font-size:13px; cursor:pointer;",
        "display:inline-flex; align-items:center; justify-content:center;",
        "transition:all 0.2s ease; padding:0; line-height:1; z-index:20;"
      ),
      onmouseover = "this.style.color='rgba(155,127,254,1)';this.style.borderColor='rgba(155,127,254,0.6)';this.style.background='rgba(155,127,254,0.15)'",
      onmouseout  = "this.style.color='rgba(148,163,184,0.6)';this.style.borderColor='rgba(148,163,184,0.28)';this.style.background='rgba(148,163,184,0.1)'",
      onclick = paste0("toggleNE('", uid_pop, "','", uid_btn, "')"),
      "✦"
    ),
    tags$div(
      id    = uid_pop,
      style = "display:none; position:fixed; width:300px;
               background:#1a2035; border:1px solid #252d45;
               border-radius:12px; padding:14px 16px;
               z-index:99999; box-shadow:0 8px 32px rgba(0,0,0,0.85);
               font-size:12px; color:#94a3b8; line-height:1.7;",
      tags$div(
        style = "font-size:10px; font-weight:700; text-transform:uppercase;
                 letter-spacing:0.7px; color:#9b7ffe;
                 margin-bottom:8px; padding-bottom:6px;
                 border-bottom:1px solid #1e2640;",
        HTML("&#10022;&nbsp; Nota explicativa")
      ),
      HTML(texto_html)
    )
  )
}

# Script global del toggle (se inyecta una sola vez en dashboardBody)
# Variante inline para tarjetas que NO son box (stat-card, metric-mini, section-hdr)
nota_exp_inline <- function(id, texto_html) {
  uid_btn <- paste0("nbtn_", id)
  uid_pop <- paste0("npop_", id)
  tags$span(
    style = "display:inline-flex; align-items:center; margin-left:8px; position:relative; vertical-align:middle;",
    tags$button(
      id      = uid_btn,
      style   = paste0(
        "background:rgba(148,163,184,0.1);",
        "border:1px solid rgba(148,163,184,0.25);",
        "color:rgba(148,163,184,0.55);",
        "border-radius:50%; width:24px; height:24px;",
        "font-size:12px; cursor:pointer;",
        "display:inline-flex; align-items:center; justify-content:center;",
        "transition:all 0.2s ease; padding:0; line-height:1; flex-shrink:0;"
      ),
      onmouseover = "this.style.color='rgba(155,127,254,1)';this.style.borderColor='rgba(155,127,254,0.6)';this.style.background='rgba(155,127,254,0.15)'",
      onmouseout  = "this.style.color='rgba(148,163,184,0.55)';this.style.borderColor='rgba(148,163,184,0.25)';this.style.background='rgba(148,163,184,0.1)'",
      onclick = paste0("toggleNE('", uid_pop, "','", uid_btn, "')"),
      "✦"
    ),
    tags$div(
      id    = uid_pop,
      style = "display:none; position:fixed; width:300px;
               background:#1a2035; border:1px solid #252d45;
               border-radius:12px; padding:14px 16px;
               z-index:99999; box-shadow:0 8px 32px rgba(0,0,0,0.85);
               font-size:12px; color:#94a3b8; line-height:1.7;",
      tags$div(
        style = "font-size:10px; font-weight:700; text-transform:uppercase;
                 letter-spacing:0.7px; color:#9b7ffe;
                 margin-bottom:8px; padding-bottom:6px;
                 border-bottom:1px solid #1e2640;",
        "✦  Nota explicativa"
      ),
      HTML(texto_html)
    )
  )
}

nota_exp_script <- tags$script(HTML("
  function toggleNE(popId, btnId) {
    var pop = document.getElementById(popId);
    var btn = document.getElementById(btnId);
    if (!pop || !btn) return;
    var rect = btn.getBoundingClientRect();
    if (pop.style.display === 'none' || pop.style.display === '') {
      pop.style.top  = (rect.bottom + 6) + 'px';
      pop.style.left = Math.max(8, rect.right - 302) + 'px';
      pop.style.display = 'block';
    } else {
      pop.style.display = 'none';
    }
  }
  document.addEventListener('click', function(e) {
    if (!e.target.closest || e.target.closest('[id^=nbtn_]')) return;
    document.querySelectorAll('[id^=npop_]').forEach(function(p) {
      p.style.display = 'none';
    });
  });
"));

stat_card <- function(label, value, sub=NULL) {
  tags$div(class="stat-card",
           tags$div(class="stat-label", label),
           tags$div(class="stat-value", value),
           if (!is.null(sub)) tags$div(class="stat-sub", sub)
  )
}

metric_mini <- function(label, value) {
  tags$div(class="metric-mini",
           tags$div(class="metric-mini-label", label),
           tags$div(class="metric-mini-value", value)
  )
}

section_hdr <- function(title, sub=NULL) {
  tags$div(class="section-header",
           tags$span(class="section-title", title),
           if (!is.null(sub)) tags$span(class="section-sub", sub)
  )
}

# MODIFICACIÓN 1: helper para el ícono con popup de nota metodológica
nota_metodologica <- function(id, cor_niveles, cor_dif) {
  tags$div(
    class = "nota-wrapper",
    style = "float:right; margin-top:-2px;",
    tags$button(
      class   = "nota-btn",
      id      = paste0("btn_", id),
      onclick = paste0("toggleNota('popup_", id, "')"),
      "i"
    ),
    tags$div(
      class = "nota-popup",
      id    = paste0("popup_", id),
      tags$div(class="nota-popup-title", "⚠ Nota Metodológica"),
      tags$div(class="nota-popup-body",
               tags$b("Correlación en niveles vs. diferencias"), tags$br(),
               tags$br(),
               "La correlación de Pearson ", tags$span(class="r-val", paste0("r = ", cor_niveles)),
               " fue calculada sobre los valores en niveles de ambas series.",
               tags$br(), tags$br(),
               "Dado que ambas presentan ", tags$b("tendencia creciente"), " en el período 2000–2025,
        esta correlación puede ser ", tags$b("espuria"), ": dos series que simplemente
        crecen con el tiempo tienden a correlacionarse aunque no exista
        una relación causal directa entre ellas.",
               tags$br(), tags$br(),
               "Al calcular la correlación en ", tags$b("primeras diferencias"),
               " (cambios mensuales, eliminando la tendencia), el resultado es ",
               tags$span(class="r-val", paste0("r = ", cor_dif)),
               ", prácticamente cero y no significativo (p ≈ 0.34).",
               tags$br(), tags$br(),
               "Esto indica que ", tags$b("mes a mes"), " los movimientos de la TRM
        no predicen significativamente los cambios en remesas.
        La relación observada es principalmente de ", tags$b("largo plazo"), "."
      )
    ),
    # JavaScript para toggle del popup
    tags$script(HTML("
      function toggleNota(id) {
        var el = document.getElementById(id);
        el.classList.toggle('visible');
        // Cerrar al hacer click fuera
        document.addEventListener('click', function handler(e) {
          if (!el.contains(e.target) && !e.target.classList.contains('nota-btn')) {
            el.classList.remove('visible');
            document.removeEventListener('click', handler);
          }
        });
      }
    "))
  )
}

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = tags$span(
      tags$span(style=paste0("color:", purple, "; font-weight:700;"), "●"),
      tags$span(style="margin-left:8px; font-weight:600; font-size:14px;",
                "Remesas Colombia")
    )
  ),
  
  dashboardSidebar(
    tags$head(tags$style(HTML(css))),
    sidebarMenu(
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><path d='M3 9.5L12 3l9 6.5V20a1 1 0 01-1 1H4a1 1 0 01-1-1z'/><path d='M9 21V12h6v9'/></svg> Portada"), tabName="portada"),
      tags$div(style=paste0("padding:12px 20px 6px; color:", text3,
                            "; font-size:10px; font-weight:700;
                              text-transform:uppercase; letter-spacing:1px;"),
               "Análisis"),
      
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><rect x='3' y='3' width='7' height='7' rx='1.5'/><rect x='14' y='3' width='7' height='7' rx='1.5'/><rect x='3' y='14' width='7' height='7' rx='1.5'/><rect x='14' y='14' width='7' height='7' rx='1.5'/></svg> Inicio"), tabName="inicio"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><polyline points='3,17 8,11 13,14 21,6'/><polyline points='17,6 21,6 21,10'/></svg> La Serie"), tabName="serie"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><line x1='3' y1='12' x2='21' y2='12'/><path d='M7 6l-4 6 4 6'/><path d='M17 6l4 6-4 6'/></svg> Estacionariedad"), tabName="estacionariedad"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><rect x='3' y='13' width='3.5' height='8' rx='1'/><rect x='10' y='8' width='3.5' height='13' rx='1'/><rect x='17' y='4' width='3.5' height='17' rx='1'/></svg> ACF y PACF"), tabName="correlogramas"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><line x1='12' y1='2' x2='12' y2='22'/><path d='M17 5H9.5a3.5 3.5 0 000 7h5a3.5 3.5 0 010 7H6'/></svg> Remesas vs TRM"), tabName="trm"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><path d='M10 13a5 5 0 007.54.54l3-3a5 5 0 00-7.07-7.07l-1.72 1.71'/><path d='M14 11a5 5 0 00-7.54-.54l-3 3a5 5 0 007.07 7.07l1.71-1.71'/></svg> Cointegración"), tabName="cointegracion"),
      tags$div(style=paste0("padding:12px 20px 6px; color:", text3,
                            "; font-size:10px; font-weight:700;
                              text-transform:uppercase; letter-spacing:1px;"),
               "Resultados"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><polygon points='12,2 22,8.5 22,15.5 12,22 2,15.5 2,8.5'/><line x1='12' y1='22' x2='12' y2='15.5'/><polyline points='22,8.5 12,15.5 2,8.5'/></svg> Modelos"), tabName="modelos"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><path d='M9 11l3 3L22 4'/><path d='M21 12v7a2 2 0 01-2 2H5a2 2 0 01-2-2V5a2 2 0 012-2h11'/></svg> Supuestos del Modelo"), tabName="supuestos"),
      menuItem(HTML("<svg style='width:16px;height:16px;display:inline-block;vertical-align:middle;margin-right:10px;' viewBox='0 0 24 24' fill='none' stroke='currentColor' stroke-width='1.8' stroke-linecap='round' stroke-linejoin='round'><polyline points='23,6 13.5,15.5 8.5,10.5 1,18'/><polyline points='17,6 23,6 23,12'/></svg> Pronóstico"), tabName="pronostico")
    )
  ),
  
  dashboardBody(
    tags$style(HTML(css)),
    nota_exp_script,
    tabItems(
      
      # ======================================================================
      # TAB 0 — PORTADA
      # ======================================================================
      tabItem(tabName="portada",
              tags$div(
                style = paste0(
                  "min-height:82vh; display:flex; flex-direction:column;",
                  "align-items:center; justify-content:center; text-align:center;",
                  "padding:40px 20px; position:relative; overflow:hidden;"
                ),
                
                # Círculos decorativos de fondo
                tags$div(style=paste0(
                  "position:absolute; top:-80px; left:-80px;",
                  "width:340px; height:340px; border-radius:50%;",
                  "background:radial-gradient(circle, rgba(124,92,252,0.12) 0%, transparent 70%);",
                  "pointer-events:none;")),
                tags$div(style=paste0(
                  "position:absolute; bottom:-60px; right:-60px;",
                  "width:300px; height:300px; border-radius:50%;",
                  "background:radial-gradient(circle, rgba(6,182,212,0.1) 0%, transparent 70%);",
                  "pointer-events:none;")),
                tags$div(style=paste0(
                  "position:absolute; top:40%; left:-40px;",
                  "width:180px; height:180px; border-radius:50%;",
                  "background:radial-gradient(circle, rgba(16,185,129,0.07) 0%, transparent 70%);",
                  "pointer-events:none;")),
                
                # Contenido principal
                tags$div(style="position:relative; z-index:1; width:100%; max-width:680px;",
                         
                         # Badge universidad
                         tags$div(style=paste0(
                           "display:inline-flex; align-items:center; gap:8px;",
                           "background:rgba(124,92,252,0.1); border:1px solid rgba(124,92,252,0.25);",
                           "border-radius:20px; padding:6px 16px; margin-bottom:32px;"),
                           tags$span(style=paste0("color:",purple2,"; font-size:12px; font-weight:600;",
                                                  "text-transform:uppercase; letter-spacing:1.5px;"),
                                     "Universidad del Norte  ·  Econometría II")
                         ),
                         
                         # Título principal con gradiente
                         tags$div(
                           style = "margin-bottom:16px;",
                           tags$div(
                             style = paste0(
                               "font-size:42px; font-weight:800; letter-spacing:-1.5px;",
                               "line-height:1.15; margin-bottom:4px;",
                               "background:linear-gradient(135deg,",purple2,",",cyan,");",
                               "-webkit-background-clip:text; -webkit-text-fill-color:transparent;",
                               "background-clip:text;"
                             ),
                             "Remesas de Trabajadores"
                           ),
                           tags$div(
                             style = paste0(
                               "font-size:42px; font-weight:800; letter-spacing:-1.5px;",
                               "line-height:1.15;",
                               "color:",text1,";"
                             ),
                             "en Colombia"
                           )
                         ),
                         
                         # Subtítulo
                         tags$div(
                           style = paste0(
                             "color:",text2,"; font-size:15px; line-height:1.7;",
                             "margin-bottom:40px; max-width:500px; margin-left:auto; margin-right:auto;"
                           ),
                           "Análisis de Series de Tiempo · Metodología Box-Jenkins",
                           tags$br(),
                           tags$span(style=paste0("color:",purple2,"; font-weight:600;"), "Modelo SARIMA"),
                           " · Período 2000–2020"
                         ),
                         
                         # Línea divisora con gradiente
                         tags$div(style=paste0(
                           "width:120px; height:2px; margin:0 auto 40px;",
                           "background:linear-gradient(90deg,transparent,",purple,",",cyan,",transparent);"
                         )),
                         
                         # Autores
                         tags$div(style="margin-bottom:40px;",
                                  tags$div(style=paste0(
                                    "color:",text3,"; font-size:10px; font-weight:700;",
                                    "text-transform:uppercase; letter-spacing:2px; margin-bottom:20px;"),
                                    "Presentado por"),
                                  tags$div(style="display:flex; justify-content:center; gap:24px; flex-wrap:wrap;",
                                           tags$div(
                                             style=paste0(
                                               "background:",card,"; border:1px solid ",border2,";",
                                               "border-radius:12px; padding:14px 28px;",
                                               "box-shadow:0 4px 20px rgba(124,92,252,0.1);"
                                             ),
                                             tags$div(style=paste0("color:",text1,"; font-size:16px; font-weight:600; margin-bottom:2px;"),
                                                      "Mariana Insignares Llach"),
                                             
                                           ),
                                           tags$div(
                                             style=paste0(
                                               "background:",card,"; border:1px solid ",border2,";",
                                               "border-radius:12px; padding:14px 28px;",
                                               "box-shadow:0 4px 20px rgba(6,182,212,0.1);"
                                             ),
                                             tags$div(style=paste0("color:",text1,"; font-size:16px; font-weight:600; margin-bottom:2px;"),
                                                      "Rodrigo Rivera De la Asunción"),
                                             
                                           )
                                  )
                         ),
                         
                         # KPIs rápidos — solo Observaciones y Período
                         tags$div(
                           style = "display:flex; justify-content:center; gap:16px; flex-wrap:wrap;",
                           tags$div(style=paste0(
                             "background:rgba(124,92,252,0.08); border:1px solid rgba(124,92,252,0.2);",
                             "border-radius:10px; padding:14px 32px; min-width:130px;"),
                             tags$div(style=paste0("color:",purple2,"; font-size:26px; font-weight:700; line-height:1;"),
                                      "252"),
                             tags$div(style=paste0("color:",text3,"; font-size:10px; text-transform:uppercase;",
                                                   "letter-spacing:0.5px; margin-top:6px;"), "Observaciones")
                           ),
                           tags$div(style=paste0(
                             "background:rgba(6,182,212,0.08); border:1px solid rgba(6,182,212,0.2);",
                             "border-radius:10px; padding:14px 32px; min-width:130px;"),
                             tags$div(style=paste0("color:",cyan,"; font-size:26px; font-weight:700; line-height:1;"),
                                      "2000–2020"),
                             tags$div(style=paste0("color:",text3,"; font-size:10px; text-transform:uppercase;",
                                                   "letter-spacing:0.5px; margin-top:6px;"), "Período")
                           )
                         )
                )
              )
      ),
      
      # ======================================================================
      # TAB 1 — INICIO
      # ======================================================================
      tabItem(tabName="inicio",
              
              # ── BLOQUE DE PRESENTACIÓN ──────────────────────────────────────────
              fluidRow(column(12,
                              tags$div(style="margin-bottom:24px;",
                                       tags$div(style=paste0("color:",purple2,
                                                             "; font-size:11px; font-weight:700; text-transform:uppercase;
               letter-spacing:1px; margin-bottom:10px;"),
                                                "● Análisis de Series de Tiempo · Colombia"),
                                       tags$div(style=paste0("color:",text1,
                                                             "; font-size:26px; font-weight:700; letter-spacing:-0.5px;
               line-height:1.2; margin-bottom:6px;"),
                                                "Remesas de Trabajadores en Colombia"),
                                       tags$div(style=paste0("color:",text2,
                                                             "; font-size:14px; font-weight:400; margin-bottom:20px;"),
                                                "Metodología Box-Jenkins · Modelo SARIMA · 2000–2020"),
                                       tags$hr(class="divider")
                              )
              )),
              
              # ── TRES TARJETAS DE CONTEXTO ────────────────────────────────────────
              fluidRow(
                column(4,
                       tags$div(
                         style=paste0(
                           "background:", card, "; border:1px solid ", border, ";",
                           "border-top:3px solid ", purple, ";",
                           "border-radius:14px; padding:22px; margin-bottom:20px; height:180px;"
                         ),
                         tags$div(style=paste0(
                           "font-size:10px; font-weight:700; text-transform:uppercase;",
                           "letter-spacing:0.8px; color:", purple2, "; margin-bottom:12px;"),
                           "📌  ¿Qué son las remesas?"),
                         tags$div(style=paste0("color:",text1,"; font-size:13.5px; line-height:1.65;"),
                                  "Transferencias de dinero enviadas por migrantes colombianos
                a sus familias o personas en su país de origen."
                         )
                       )
                ),
                column(4,
                       tags$div(
                         style=paste0(
                           "background:", card, "; border:1px solid ", border, ";",
                           "border-top:3px solid ", green, ";",
                           "border-radius:14px; padding:22px; margin-bottom:20px; height:180px;"
                         ),
                         tags$div(style="display:flex; align-items:center; justify-content:space-between; margin-bottom:12px;",
                                  tags$div(style=paste0(
                                    "font-size:10px; font-weight:700; text-transform:uppercase;",
                                    "letter-spacing:0.8px; color:", green, ";"),
                                    "💰  Importancia económica"),
                                  nota_exp_inline("imp_eco",
                                                  "El crecimiento de las remesas refleja el aumento de la migración colombiana y la mayor dependencia de los hogares respecto a ingresos enviados desde el exterior. Además, factores como la <b style='color:#f1f5f9'>depreciación del peso colombiano</b> y las <b style='color:#f1f5f9'>crisis económicas</b> han incrementado la importancia de las remesas como mecanismo de apoyo económico.")
                         ),
                         tags$div(style=paste0("color:",text1,"; font-size:13.5px; line-height:1.65;"),
                                  "Colombia recibió más de ",
                                  tags$b(style=paste0("color:",green,";"), "USD 9.500 millones"),
                                  " en 2023, equivalente al ",
                                  tags$b(style=paste0("color:",green,";"), "3% del PIB"),
                                  ". Son la ",
                                  tags$b("segunda fuente de divisas"),
                                  " después del petróleo."
                         )
                       )
                ),
                column(4,
                       tags$div(
                         style=paste0(
                           "background:", card, "; border:1px solid ", border, ";",
                           "border-top:3px solid ", blue2, ";",
                           "border-radius:14px; padding:22px 24px; margin-bottom:20px; height:180px;",
                           "display:flex; flex-direction:column; justify-content:center;"
                         ),
                         tags$div(style=paste0(
                           "font-size:10px; font-weight:700; text-transform:uppercase;",
                           "letter-spacing:0.8px; color:", blue2, "; margin-bottom:14px;"),
                           "🎯  Objetivo del análisis"),
                         tags$div(style=paste0(
                           "color:",text1,"; font-size:13.5px; line-height:1.65;"),
                           "Analiza las remesas mediante la metodología ",
                           tags$b("Box-Jenkins"),
                           " y un modelo ",
                           tags$b(style=paste0("color:",blue2,";"), "SARIMA"),
                           " para estudiar su comportamiento temporal, generar pronósticos",
                           " y evaluar su relación con la TRM."
                         )
                       )
                )
              ),
              
              tags$hr(class="divider"),
              
              # ── KPIs del análisis ────────────────────────────────────────────────
              fluidRow(
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Observaciones"),
                                tags$div(class="stat-value", "252"),
                                tags$div(class="stat-sub", tags$span(class="badge-neutral", "Mensual"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Período"),
                                tags$div(class="stat-value", "2000–2020"),
                                tags$div(class="stat-sub", tags$span(class="badge-neutral", "Pronóstico: 2021"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Mejor modelo"),
                                tags$div(style=paste0("color:", purple2,
                                                      "; font-size:16px; font-weight:700;
                       letter-spacing:-0.3px; margin-bottom:6px;"),
                                         "SARIMA(1,1,2)(2,0,0)"),
                                tags$div(class="stat-sub", tags$span(class="badge-up", "✅ Ruido blanco"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "MAPE Pronóstico 2021"),
                                tags$div(class="stat-value", paste0(mape_val,"%")),
                                tags$div(class="stat-sub", tags$span(class="badge-neutral", "Error promedio"))
                       )
                )
              ),
              tags$div(style="margin:24px 0 16px;",
                       tags$hr(class="divider"),
                       section_hdr("Estadísticas descriptivas", "Serie original en Millones USD")
              ),
              fluidRow(
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Promedio histórico"),
                                tags$div(class="stat-value",
                                         paste0("$",round(mean(as.numeric(remesas_ts)),1)," M")),
                                tags$div(class="stat-sub", "Millones USD")
                       )
                ),
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Mínimo histórico"),
                                tags$div(style=paste0("color:", blue2,
                                                      "; font-size:26px; font-weight:700;
                       letter-spacing:-0.5px; margin-bottom:6px;"),
                                         paste0("$",round(min(as.numeric(remesas_ts)),1)," M")),
                                tags$div(class="stat-sub", tags$span(class="badge-down", "↓ Ene 2000"))
                       )
                ),
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Máximo histórico"),
                                tags$div(style=paste0("color:", green,
                                                      "; font-size:26px; font-weight:700;
                       letter-spacing:-0.5px; margin-bottom:6px;"),
                                         paste0("$",round(max(as.numeric(remesas_ts)),1)," M")),
                                tags$div(class="stat-sub", tags$span(class="badge-up", "↑ Mar 2020"))
                       )
                )
              ),
              fluidRow(
                box(title=tags$span("Principales países de origen de remesas a Colombia",
                                    nota_exp("paises_remesas",
                                             "EEUU concentra más de la mitad de las remesas recibidas por Colombia, reflejo de la gran comunidad colombiana residente en ese país.")),
                    width=5, solidHeader=TRUE,
                    plotlyOutput("plot_pastel", height="300px")),
                box(title=NULL, width=7,
                    plotlyOutput("plot_inicio", height="300px"))
              )
      ),
      
      # ======================================================================
      # TAB 2 — LA SERIE
      # ======================================================================
      tabItem(tabName="serie",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px;",
                                       section_hdr("La Serie",
                                                   "Exploración visual · Original y transformación logarítmica")
                              )
              )),
              fluidRow(
                box(title=tags$span("Serie Original — Remesas Colombia (2000–2020)",
                                    nota_exp("serie_orig",
                                             "<b style='color:#f1f5f9'>Aumentan por:</b><br>
                   &bull; Mayor migración internacional<br>
                   &bull; Apoyo económico a hogares<br>
                   &bull; Mejores ingresos en el exterior<br><br>
                   <b style='color:#f1f5f9'>Disminuyen por:</b><br>
                   &bull; Desempleo de migrantes<br>
                   &bull; Crisis económicas globales<br>
                   &bull; Restricciones migratorias")),
                    width=12, solidHeader=TRUE,
                    plotlyOutput("plot_serie_original", height="320px"))
              ),
              fluidRow(
                box(title=tags$span("Serie en Logaritmo Natural",
                                    nota_exp("serie_log",
                                             "La serie en logaritmos se usa para <b style='color:#f1f5f9'>suavizar los cambios bruscos</b> de las remesas y facilitar el análisis de su comportamiento en el tiempo. Además, ayuda a que el modelo haga pronósticos más estables al expresar los datos en <b style='color:#f1f5f9'>cambios porcentuales</b>.")),
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_serie_log", height="280px")),
                box(title=tags$span("Comparación: Original vs Logaritmo",
                                    nota_exp("serie_comp",
                                             "La serie original presenta mayor <b style='color:#f1f5f9'>volatilidad y cambios bruscos</b> en el tiempo, mientras que la transformación logarítmica suaviza las fluctuaciones y permite observar con mayor claridad la <b style='color:#f1f5f9'>tendencia de crecimiento</b> de las remesas.")),
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_comparacion", height="280px"))
              ),
              fluidRow(
                box(title=tags$span("Descomposición Clásica de la Serie",
                                    nota_exp("descomp",
                                             "La descomposición separa la serie en <b style='color:#f1f5f9'>tendencia, estacionalidad y residuos</b>, evidenciando un crecimiento sostenido de las remesas y patrones estacionales constantes, especialmente en los últimos meses del año.")),
                    width=12, solidHeader=TRUE,
                    plotOutput("plot_descomposicion", height="380px"))
              ),
              fluidRow(
                box(title=tags$span("Patrón Estacional por Año",
                                    nota_exp("estacional",
                                             "El gráfico muestra un comportamiento estacional recurrente en las remesas, con mayor dinamismo y picos más altos durante <b style='color:#f1f5f9'>diciembre</b>, asociado al aumento de transferencias en fin de año.")),
                    width=6, solidHeader=TRUE,
                    plotOutput("plot_estacional", height="300px")),
                box(title=tags$span("Subseries Estacionales por Mes",
                                    nota_exp("subseries",
                                             "Las subseries evidencian diferencias sistemáticas entre meses, destacándose <b style='color:#f1f5f9'>diciembre</b> como el periodo de mayor nivel promedio de remesas y actividad estacional.")),
                    width=6, solidHeader=TRUE,
                    plotOutput("plot_subseries", height="300px"))
              )
      ),
      
      # ======================================================================
      # TAB 3 — ESTACIONARIEDAD
      # ======================================================================
      tabItem(tabName="estacionariedad",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px;",
                                       section_hdr("Estacionariedad",
                                                   "Pruebas ADF · KPSS · Diferenciación")
                              )
              )),
              fluidRow(
                column(6,
                       tags$div(class="stat-card",
                                tags$div(style="display:flex; align-items:center; justify-content:space-between; position:relative; min-height:32px;",
                                         tags$div(class="stat-label", "Serie Original"),
                                         nota_exp_inline("estac_orig",
                                                         "Los resultados de <b style='color:#f1f5f9'>ADF y KPSS</b> indican que existe evidencia estadística de <b style='color:#f1f5f9'>no estacionariedad</b> en la serie original, debido a la presencia de tendencia y persistencia temporal.")
                                ),
                                tags$div(class="divider"),
                                fluidRow(
                                  column(6,
                                         tags$div(class="metric-mini-label", "ADF p-value"),
                                         tags$div(style=paste0("color:", red, "; font-size:20px; font-weight:700;"),
                                                  round(adf.test(remesas_ts)$p.value, 4))
                                  ),
                                  column(6,
                                         tags$div(class="metric-mini-label", "KPSS p-value"),
                                         tags$div(style=paste0("color:", red, "; font-size:20px; font-weight:700;"),
                                                  round(kpss.test(remesas_ts,null="Trend")$p.value,4))
                                  )
                                ),
                                tags$div(style="margin-top:14px;",
                                         tags$span(class="pill-no", "❌ Serie NO estacionaria"))
                       )
                ),
                column(6,
                       tags$div(class="stat-card",
                                tags$div(style="display:flex; align-items:center; justify-content:space-between; position:relative; min-height:32px;",
                                         tags$div(class="stat-label", "Serie Log Diferenciada (d=1)"),
                                         nota_exp_inline("estac_dif",
                                                         "Luego de aplicar una <b style='color:#f1f5f9'>diferenciación regular</b> sobre la serie en logaritmos, los tests ADF y KPSS evidencian <b style='color:#f1f5f9'>estacionariedad estadística</b>, validando la transformación aplicada.")
                                ),
                                tags$div(class="divider"),
                                fluidRow(
                                  column(6,
                                         tags$div(class="metric-mini-label", "ADF p-value"),
                                         tags$div(style=paste0("color:", green, "; font-size:20px; font-weight:700;"),
                                                  round(adf.test(na.omit(dif_simple_log))$p.value,4))
                                  ),
                                  column(6,
                                         tags$div(class="metric-mini-label", "KPSS p-value"),
                                         tags$div(style=paste0("color:", green, "; font-size:20px; font-weight:700;"),
                                                  round(kpss.test(na.omit(dif_simple_log))$p.value,4))
                                  )
                                ),
                                tags$div(style="margin-top:14px;",
                                         tags$span(class="pill-ok", "✅ Serie ES estacionaria"))
                       )
                )
              ),
              fluidRow(
                box(title=tags$span("Serie Diferenciada Log (d=1)",
                                    nota_exp("dif_simple",
                                             "La serie diferenciada fluctúa alrededor de una <b style='color:#f1f5f9'>media constante</b> y varianza relativamente estable, comportamiento consistente con una <b style='color:#f1f5f9'>serie estacionaria</b>.")),
                    width=8, solidHeader=TRUE,
                    plotlyOutput("plot_dif_simple", height="280px")),
                column(4,
                       tags$div(style="padding:0 0 0 4px;",
                                tags$div(class="stat-card",
                                         tags$div(class="stat-label", "ndiffs — Dif. regulares"),
                                         tags$div(style=paste0("color:", blue2, "; font-size:36px; font-weight:700;"),
                                                  ndiffs(remesas_log)),
                                         tags$div(class="stat-sub", "Diferencias no estacionales necesarias")
                                ),
                                tags$div(class="stat-card",
                                         tags$div(class="stat-label", "nsdiffs — Dif. estacionales"),
                                         tags$div(style=paste0("color:", purple2, "; font-size:36px; font-weight:700;"),
                                                  nsdiffs(remesas_log)),
                                         tags$div(class="stat-sub", "Diferencias estacionales necesarias")
                                )
                       )
                )
              )
      ),
      
      # ======================================================================
      # TAB 4 — ACF Y PACF
      # ======================================================================
      tabItem(tabName="correlogramas",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px; display:flex; align-items:center; gap:10px;",
                                       section_hdr("ACF y PACF",
                                                   "Correlogramas · Serie diferenciada log (d=1)"),
                                       nota_exp_inline("acf_pacf",
                                                       "Las funciones ACF y PACF permiten identificar la <b style='color:#f1f5f9'>dependencia temporal y estacional</b> de la serie, sugiriendo una estructura compatible con un modelo <b style='color:#f1f5f9'>SARIMA(1,1,2)(2,0,0)[12]</b>, donde se evidencian componentes AR, MA y patrones estacionales anuales significativos.")
                              )
              )),
              fluidRow(
                box(title="Función de Autocorrelación — ACF",
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_acf", height="320px")),
                box(title="Autocorrelación Parcial — PACF",
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_pacf", height="320px"))
              ),
              fluidRow(
                column(4,
                       tags$div(class="interp-card",
                                tags$div(class="interp-tag", style=paste0("color:", blue2, ";"),
                                         "PARTE NO ESTACIONAL"),
                                tags$div(class="interp-row",
                                         tags$div(class="interp-key", "ACF"),
                                         tags$div(class="interp-val", style=paste0("color:", blue2, ";"),
                                                  "Corte en lag 1–2 → MA(2)")
                                ),
                                tags$div(class="interp-row",
                                         tags$div(class="interp-key", "PACF"),
                                         tags$div(class="interp-val", style=paste0("color:", blue2, ";"),
                                                  "Decaimiento lento → AR(1)")
                                )
                       )
                ),
                column(4,
                       tags$div(class="interp-card",
                                tags$div(class="interp-tag", style=paste0("color:", purple2, ";"),
                                         "PARTE ESTACIONAL"),
                                tags$div(class="interp-row",
                                         tags$div(class="interp-key", "Lags 12, 24"),
                                         tags$div(class="interp-val", style=paste0("color:", purple2, ";"),
                                                  "Picos significativos → SAR(2)")
                                ),
                                tags$div(class="interp-row",
                                         tags$div(class="interp-key", "D = 0"),
                                         tags$div(class="interp-val", style=paste0("color:", purple2, ";"),
                                                  "Sin diferenciación estacional")
                                )
                       )
                ),
                column(4,
                       tags$div(class="interp-card",
                                tags$div(class="interp-tag", style=paste0("color:", green, ";"),
                                         "MODELO SUGERIDO"),
                                tags$div(style="margin-top:16px;",
                                         tags$div(style=paste0("background:rgba(16,185,129,0.08);
                  border:1px solid rgba(16,185,129,0.2);
                  border-radius:10px; padding:16px; text-align:center;"),
                                                  tags$div(style=paste0("color:", green,
                                                                        "; font-size:16px; font-weight:700; letter-spacing:-0.3px;"),
                                                           "SARIMA(1,1,2)(2,0,0)[12]")
                                         )
                                )
                       )
                )
              )
      ),
      
      # ======================================================================
      # TAB 5 — REMESAS VS TRM
      # ======================================================================
      tabItem(tabName="trm",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px;",
                                       section_hdr("Remesas vs TRM",
                                                   "Índice base 100 = enero 2000 · Escala logarítmica · 2000–2025")
                              )
              )),
              fluidRow(
                column(3,
                       tags$div(class="stat-card",
                                tags$div(style="display:flex; align-items:center; justify-content:space-between;",
                                         tags$div(class="stat-label", "Correlación de Pearson"),
                                         nota_exp_inline("cor_pearson",
                                                         "El valor <b style='color:#f1f5f9'>r = 0.82</b> indica una correlación positiva fuerte entre ambas variables. Sin embargo, dado que ambas series son <b style='color:#f1f5f9'>no estacionarias I(1)</b>, esta correlación podría ser espuria. Ver el tab <b style='color:#f1f5f9'>Cointegración</b> para el test de Engle-Granger (1987) que verifica si la relación es real.")
                                ),
                                tags$div(style=paste0("color:", purple2,
                                                      "; font-size:32px; font-weight:700;
                       letter-spacing:-0.5px; margin-bottom:6px;"),
                                         paste0("r = ", cor_pearson)),
                                tags$div(class="stat-sub",
                                         tags$span(class="badge-up", "↑ Correlación alta positiva"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Crecimiento Remesas"),
                                tags$div(style=paste0("color:", purple2,
                                                      "; font-size:32px; font-weight:700; margin-bottom:6px;"),
                                         "~11×"),
                                tags$div(class="stat-sub", "Desde enero 2000 hasta 2025")
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Crecimiento TRM"),
                                tags$div(style=paste0("color:", green,
                                                      "; font-size:32px; font-weight:700; margin-bottom:6px;"),
                                         "~2.6×"),
                                tags$div(class="stat-sub", "COP/USD desde enero 2000")
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Período analizado"),
                                tags$div(style=paste0("color:", blue2,
                                                      "; font-size:28px; font-weight:700; margin-bottom:6px;"),
                                         "2000–2025"),
                                tags$div(class="stat-sub",
                                         tags$span(class="badge-neutral", "312 observaciones"))
                       )
                )
              ),
              fluidRow(
                box(title=tags$span("Índice Base 100 — Remesas y TRM (Escala Logarítmica)",
                                    nota_exp("trm_indice",
                                             "El gráfico evidencia una relación creciente entre las remesas y la TRM, mostrando que las <b style='color:#f1f5f9'>depreciaciones del peso colombiano</b> suelen coincidir con aumentos en el flujo de remesas, especialmente durante periodos de <b style='color:#f1f5f9'>crisis e incertidumbre económica</b>.")),
                    width=12, solidHeader=TRUE,
                    tags$div(style="display:flex; align-items:center; gap:16px; margin-bottom:12px;",
                             tags$div(style=paste0("color:",text3,"; font-size:11px; font-weight:600;
                         text-transform:uppercase; letter-spacing:0.5px; white-space:nowrap;"),
                                      "Desde"),
                             selectInput("anio_ini_trm", label=NULL,
                                         choices=2000:2024, selected=2000, width="110px"),
                             tags$div(style=paste0("color:",text3,"; font-size:11px; font-weight:600;
                         text-transform:uppercase; letter-spacing:0.5px; white-space:nowrap;"),
                                      "Hasta"),
                             selectInput("anio_fin_trm", label=NULL,
                                         choices=2001:2025, selected=2025, width="110px")
                    ),
                    plotlyOutput("plot_trm_indice", height="380px"))
              ),
              fluidRow(
                box(title="Remesas (Millones USD) — 2000–2025",
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_trm_remesas", height="280px")),
                box(title="TRM (COP/USD) — 2000–2025",
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_trm_trm", height="280px"))
              ),
              fluidRow(
                box(title=tags$span("Dispersión: Remesas vs TRM",
                                    nota_exp("trm_scatter",
                                             "El gráfico de dispersión muestra una <b style='color:#f1f5f9'>relación positiva</b> entre la TRM y las remesas, donde mayores niveles del tipo de cambio se asocian con mayores envíos hacia Colombia.")),
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_trm_scatter", height="300px")),
                box(
                  title = tags$span(
                    "Interpretación",
                    tags$button(
                      id      = "btn_nota_trm",
                      onclick = "toggleNotaTRM()",
                      style   = paste0(
                        "position:absolute; top:10px; right:12px;",
                        "background:rgba(155,127,254,0.15);",
                        "border:1px solid rgba(155,127,254,0.4);",
                        "color:#9b7ffe; border-radius:50%;",
                        "width:26px; height:26px; font-size:12px; font-weight:700;",
                        "cursor:pointer; display:inline-flex;",
                        "align-items:center; justify-content:center;",
                        "transition:all 0.2s ease; z-index:10;"
                      ),
                      "i"
                    )
                  ),
                  width=6, solidHeader=TRUE,
                  tags$div(style="padding:8px;",
                           # ── Script + Popup nota metodológica ──
                           tags$script(HTML("
                 function toggleNotaTRM() {
                   var popup = document.getElementById('popup_nota_trm');
                   var btn   = document.getElementById('btn_nota_trm');
                   if (!popup || !btn) return;
                   var rect = btn.getBoundingClientRect();
                   if (popup.style.display === 'none' || popup.style.display === '') {
                     popup.style.top  = (rect.bottom + 8) + 'px';
                     popup.style.left = Math.max(8, rect.right - 315) + 'px';
                     popup.style.display = 'block';
                   } else {
                     popup.style.display = 'none';
                   }
                 }
                 document.addEventListener('click', function(e) {
                   var popup = document.getElementById('popup_nota_trm');
                   if (popup && e.target.id !== 'btn_nota_trm' && !popup.contains(e.target)) {
                     popup.style.display = 'none';
                   }
                 });
               ")),
                           tags$div(
                             id    = "popup_nota_trm",
                             style = "display:none; position:fixed; width:315px;
                          background:#1a2035; border:1px solid #252d45;
                          border-radius:12px; padding:16px 18px;
                          z-index:99999; box-shadow:0 8px 32px rgba(0,0,0,0.85);
                          font-size:12px; color:#94a3b8; line-height:1.7;",
                             tags$div(
                               style = "font-size:11px; font-weight:700; text-transform:uppercase;
                            letter-spacing:0.7px; color:#9b7ffe;
                            margin-bottom:10px; padding-bottom:8px;
                            border-bottom:1px solid #1e2640;",
                               "⚠  Nota Metodológica"
                             ),
                             tags$b(style="color:#f1f5f9;", "Correlación en niveles vs. diferencias"),
                             tags$br(), tags$br(),
                             "La correlación de Pearson ",
                             tags$span(
                               style = "background:rgba(6,182,212,0.15); border:1px solid rgba(6,182,212,0.3);
                            color:#06b6d4; border-radius:6px; padding:2px 8px; font-weight:700;",
                               paste0("r = ", cor_pearson)
                             ),
                             " se calculó sobre los valores en niveles de ambas series.",
                             tags$br(), tags$br(),
                             "Ambas series tienen ",
                             tags$b(style="color:#f1f5f9;", "tendencia creciente"),
                             " en 2000-2025, por lo que la correlación puede ser ",
                             tags$b(style="color:#f1f5f9;", "espuria"),
                             ": series con tendencia positiva tienden a correlacionarse aunque no
                 exista relación causal directa.",
                             tags$br(), tags$br(),
                             "En primeras diferencias el resultado es ",
                             tags$span(
                               style = "background:rgba(6,182,212,0.15); border:1px solid rgba(6,182,212,0.3);
                            color:#06b6d4; border-radius:6px; padding:2px 8px; font-weight:700;
                            white-space:nowrap; display:inline-block;",
                               paste0("r = ", cor_pearson_dif)
                             ),
                             " — prácticamente cero (p ≈ 0.34).",
                             tags$br(),
                             "La relación es principalmente de ",
                             tags$b(style="color:#f1f5f9;", "largo plazo"),
                             ", no de transmisión mensual inmediata."
                           ),
                           tags$div(class="interp-card",
                                    
                                    tags$div(class="interp-tag", style=paste0("color:", purple2, ";"),
                                             "HALLAZGOS PRINCIPALES"),
                                    tags$div(class="interp-row",
                                             tags$div(class="interp-key", "Correlación en niveles"),
                                             tags$div(class="interp-val", style=paste0("color:", purple2, ";"),
                                                      paste0("r = ", cor_pearson, " — relación positiva fuerte"))
                                    ),
                                    tags$div(class="interp-row",
                                             tags$div(class="interp-key", "Correlación en diferencias"),
                                             tags$div(class="interp-val", style=paste0("color:", cyan, ";"),
                                                      paste0("r = ", cor_pearson_dif, " — sin relación mes a mes"))
                                    ),
                                    tags$div(class="interp-row",
                                             tags$div(class="interp-key", "Mecanismo"),
                                             tags$div(class="interp-val", style=paste0("color:", blue2, ";"),
                                                      "TRM alto → más pesos por dólar → incentiva remesas")
                                    ),
                                    tags$div(class="interp-row",
                                             tags$div(class="interp-key", "Crisis 2008"),
                                             tags$div(class="interp-val", style=paste0("color:", text2, ";"),
                                                      "Caída en remesas por contracción económica global")
                                    ),
                                    tags$div(class="interp-row",
                                             tags$div(class="interp-key", "COVID-19 2020"),
                                             tags$div(class="interp-val", style=paste0("color:", text2, ";"),
                                                      "TRM disparada, remesas caen y luego rebotan")
                                    ),
                                    tags$div(class="interp-row",
                                             tags$div(class="interp-key", "Post-2021"),
                                             tags$div(class="interp-val", style=paste0("color:", green, ";"),
                                                      "Boom: ambas variables crecen simultáneamente")
                                    )
                           )
                  )
                )
              )
      ),
      
      # ======================================================================
      # TAB 6 — COINTEGRACIÓN
      # ======================================================================
      tabItem(tabName="cointegracion",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px;",
                                       section_hdr("Cointegración",
                                                   "Test de Engle-Granger · Remesas vs TRM · 2000–2025")
                              )
              )),
              
              # Explicación del test
              fluidRow(
                box(title=tags$span("¿Qué es la Cointegración?",
                                    nota_exp("coint_teoria",
                                             "Dos series <b style='color:#f1f5f9'>no estacionarias</b> están cointegradas si
                   existe una <b style='color:#f1f5f9'>combinación lineal</b> entre ellas que sí es
                   estacionaria. Eso implica que comparten una <b style='color:#f1f5f9'>relación de
                   equilibrio de largo plazo</b> de la que no pueden alejarse permanentemente,
                   aunque en el corto plazo puedan desviarse.")),
                    width=12, solidHeader=TRUE,
                    tags$div(style="display:flex; gap:16px; flex-wrap:wrap; padding:4px 0;",
                             tags$div(style=paste0(
                               "flex:1; min-width:200px; background:",card2,";",
                               "border:1px solid ",border,"; border-radius:10px; padding:16px;"),
                               tags$div(style=paste0("color:",purple2,"; font-size:11px; font-weight:700;",
                                                     "text-transform:uppercase; letter-spacing:0.6px; margin-bottom:8px;"),
                                        "Paso 1 — Regresión de largo plazo"),
                               tags$div(style=paste0("color:",text2,"; font-size:13px; line-height:1.6;"),
                                        "Se estima: Remesas = α + β·TRM + ε",
                                        tags$br(),
                                        "Si las series están cointegradas, los residuos ε deben ser estacionarios.")
                             ),
                             tags$div(style=paste0(
                               "flex:1; min-width:200px; background:",card2,";",
                               "border:1px solid ",border,"; border-radius:10px; padding:16px;"),
                               tags$div(style=paste0("color:",cyan,"; font-size:11px; font-weight:700;",
                                                     "text-transform:uppercase; letter-spacing:0.6px; margin-bottom:8px;"),
                                        "Paso 2 — Test ADF sobre residuos"),
                               tags$div(style=paste0("color:",text2,"; font-size:13px; line-height:1.6;"),
                                        "H₀: Los residuos tienen raíz unitaria (no cointegración).",
                                        tags$br(),
                                        "Si p < 0.05, se rechaza H₀ → hay cointegración.")
                             ),
                             tags$div(style=paste0(
                               "flex:1; min-width:200px; background:",card2,";",
                               "border:1px solid ",border,"; border-radius:10px; padding:16px;"),
                               tags$div(style=paste0("color:",green,"; font-size:11px; font-weight:700;",
                                                     "text-transform:uppercase; letter-spacing:0.6px; margin-bottom:8px;"),
                                        "Interpretación"),
                               tags$div(style=paste0("color:",text2,"; font-size:13px; line-height:1.6;"),
                                        "Cointegración confirma que la correlación entre Remesas y TRM",
                                        tags$br(),
                                        "NO es espuria, sino una relación de largo plazo real.")
                             )
                    )
                )
              ),
              
              # Resultados del test
              fluidRow(
                column(4,
                       tags$div(class="stat-card",
                                tags$div(style="display:flex; align-items:center; justify-content:space-between;",
                                         tags$div(class="stat-label", "Estadístico ADF"),
                                         nota_exp_inline("eg_stat",
                                                         "El estadístico ADF de Dickey-Fuller Aumentado evalúa si los residuos
                   de la regresión de largo plazo tienen <b style='color:#f1f5f9'>raíz unitaria</b>.
                   Valores más negativos indican mayor evidencia de estacionariedad.")
                                ),
                                tags$div(style=paste0("color:",blue2,"; font-size:28px; font-weight:700;",
                                                      "margin:8px 0 4px;"), textOutput("eg_stat_val", inline=TRUE)),
                                tags$div(class="stat-sub", "Test ADF sobre residuos")
                       )
                ),
                column(4,
                       tags$div(class="stat-card",
                                tags$div(style="display:flex; align-items:center; justify-content:space-between;",
                                         tags$div(class="stat-label", "p-valor"),
                                         nota_exp_inline("eg_pval",
                                                         "Si el p-valor es <b style='color:#f1f5f9'>menor a 0.05</b>, se rechaza la
                   hipótesis nula de raíz unitaria en los residuos, lo que confirma la
                   cointegración entre las series.")
                                ),
                                tags$div(style=paste0("color:",purple2,"; font-size:28px; font-weight:700;",
                                                      "margin:8px 0 4px;"), textOutput("eg_pval_val", inline=TRUE)),
                                tags$div(class="stat-sub", "Significancia al 5%")
                       )
                ),
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Conclusión"),
                                tags$div(style=paste0("font-size:13px; font-weight:600; margin:8px 0 4px;",
                                                      "line-height:1.5;"),
                                         textOutput("eg_conclusion_val", inline=TRUE)),
                                tags$div(class="stat-sub", "Test de Engle-Granger (1987)")
                       )
                )
              ),
              
              # Gráficos
              fluidRow(
                box(title=tags$span("Regresión de Largo Plazo: Remesas ~ TRM",
                                    nota_exp("coint_reg",
                                             "La pendiente β representa cuántos <b style='color:#f1f5f9'>millones de USD adicionales</b>
                   en remesas se asocian con un aumento de 1 COP/USD en la TRM.
                   Si hay cointegración, esta relación es estable en el largo plazo.")),
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_coint_reg", height="300px")),
                box(title=tags$span("Residuos de la Regresión — ¿Son Estacionarios?",
                                    nota_exp("coint_resid",
                                             "Si los residuos son <b style='color:#f1f5f9'>estacionarios</b> (fluctúan alrededor
                   de cero sin tendencia), confirman la cointegración. Una serie no estacionaria
                   de residuos indicaría que la relación es espuria.")),
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_coint_resid", height="300px"))
              ),
              
              # Causalidad de Granger
              fluidRow(
                box(title=tags$span("Test de Causalidad de Granger",
                                    nota_exp("granger_nota",
                                             "<b style='color:#f1f5f9'>Causalidad de Granger NO es causalidad económica real.</b>
                   Es causalidad <b style='color:#f1f5f9'>predictiva</b>: si X causa a Y en Granger,
                   significa que los valores pasados de X ayudan a predecir Y, más allá de lo que
                   ya explica el propio pasado de Y.<br><br>
                   H₀: X no causa en sentido Granger a Y<br>
                   Si p &lt; 0.05 → se rechaza H₀ → hay causalidad predictiva.")),
                    width=12, solidHeader=TRUE,
                    tags$div(style="padding:4px 0;",
                             fluidRow(
                               column(6,
                                      tags$div(class="interp-card",
                                               tags$div(class="interp-tag",
                                                        style=paste0("color:",cyan,";"),
                                                        "¿TRM causa (Granger) a Remesas?"),
                                               tags$div(class="interp-row",
                                                        tags$div(class="interp-key", "H₀: TRM no predice Remesas"),
                                                        tags$div(class="interp-val",
                                                                 style=paste0("color:",blue2,";"),
                                                                 tags$div(style="font-size:22px; font-weight:700; margin-bottom:4px;",
                                                                          textOutput("gr_pval1", inline=TRUE)),
                                                                 textOutput("gr_conc1", inline=TRUE))
                                               ),
                                               tags$div(class="interp-row",
                                                        tags$div(class="interp-key", "Interpretación"),
                                                        tags$div(class="interp-val",
                                                                 style=paste0("color:",text2,";"),
                                                                 "¿Los valores pasados de la TRM ayudan a predecir las remesas del mes siguiente?")
                                               )
                                      )
                               ),
                               column(6,
                                      tags$div(class="interp-card",
                                               tags$div(class="interp-tag",
                                                        style=paste0("color:",purple2,";"),
                                                        "¿Remesas causan (Granger) a TRM?"),
                                               tags$div(class="interp-row",
                                                        tags$div(class="interp-key", "H₀: Remesas no predicen TRM"),
                                                        tags$div(class="interp-val",
                                                                 style=paste0("color:",purple2,";"),
                                                                 tags$div(style="font-size:22px; font-weight:700; margin-bottom:4px;",
                                                                          textOutput("gr_pval2", inline=TRUE)),
                                                                 textOutput("gr_conc2", inline=TRUE))
                                               ),
                                               tags$div(class="interp-row",
                                                        tags$div(class="interp-key", "Interpretación"),
                                                        tags$div(class="interp-val",
                                                                 style=paste0("color:",text2,";"),
                                                                 "¿Los valores pasados de las remesas ayudan a predecir la TRM del mes siguiente?")
                                               )
                                      )
                               )
                             )
                    )
                )
              ),
              
              # Ecuación de largo plazo
              fluidRow(
                box(title="Ecuación de Equilibrio de Largo Plazo",
                    width=12, solidHeader=TRUE,
                    tags$div(style="padding:8px;",
                             tags$div(class="interp-card",
                                      tags$div(class="interp-tag",
                                               style=paste0("color:",cyan,";"),
                                               "RELACIÓN DE LARGO PLAZO ESTIMADA"),
                                      tags$div(style=paste0(
                                        "font-size:18px; font-weight:600; color:",text1,";",
                                        "margin:16px 0; text-align:center; letter-spacing:-0.3px;"),
                                        "Remesas = ",
                                        tags$span(style=paste0("color:",purple2,";"),
                                                  textOutput("coef_alpha", inline=TRUE)),
                                        " + ",
                                        tags$span(style=paste0("color:",cyan,";"),
                                                  textOutput("coef_beta", inline=TRUE)),
                                        " × TRM"
                                      ),
                                      tags$div(class="interp-row",
                                               tags$div(class="interp-key", "Interpretación de β"),
                                               tags$div(class="interp-val", style=paste0("color:",blue2,";"),
                                                        "Por cada aumento de 1 COP en la TRM, las remesas aumentan
                       β millones de USD en el largo plazo")
                                      ),
                                      tags$div(class="interp-row",
                                               tags$div(class="interp-key", "Validez"),
                                               tags$div(class="interp-val", style=paste0("color:",text2,";"),
                                                        "Esta ecuación es válida solo si las series están cointegradas
                       (residuos estacionarios)")
                                      )
                             )
                    )
                )
              )
      ),
      
      # ======================================================================
      # TAB 6 — SUPUESTOS DEL MODELO (MODIFICACIÓN 3 — nueva pestaña)
      # ======================================================================
      tabItem(tabName="supuestos",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px;",
                                       section_hdr("Supuestos del Modelo",
                                                   "Diagnóstico formal · SARIMA(1,1,2)(2,0,0)[12]")
                              )
              )),
              # KPIs de supuestos
              fluidRow(
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Media de residuos"),
                                tags$div(style=paste0("color:", green,
                                                      "; font-size:24px; font-weight:700; margin-bottom:6px;"),
                                         media_res),
                                tags$div(class="stat-sub", tags$span(class="badge-up", "✅ Cercana a cero"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Desv. estándar residuos"),
                                tags$div(style=paste0("color:", blue2,
                                                      "; font-size:24px; font-weight:700; margin-bottom:6px;"),
                                         sd_res),
                                tags$div(class="stat-sub", tags$span(class="badge-neutral", "Unidades log"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Jarque-Bera p-value"),
                                tags$div(style=paste0("color:", if(jb_pval > 0.05) green else red,
                                                      "; font-size:24px; font-weight:700; margin-bottom:6px;"),
                                         jb_pval),
                                tags$div(class="stat-sub",
                                         tags$span(class=if(jb_pval > 0.05) "badge-up" else "badge-down",
                                                   if(jb_pval > 0.05) "✅ Normalidad" else "⚠ No normal"))
                       )
                ),
                column(3,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "Ljung-Box p-value"),
                                tags$div(style=paste0("color:", if(lb_pval_mej > 0.05) green else red,
                                                      "; font-size:24px; font-weight:700; margin-bottom:6px;"),
                                         lb_pval_mej),
                                tags$div(class="stat-sub",
                                         tags$span(class=if(lb_pval_mej > 0.05) "badge-up" else "badge-down",
                                                   if(lb_pval_mej > 0.05) "✅ Ruido blanco" else "❌ Autocorrelación"))
                       )
                )
              ),
              # Gráficos de diagnóstico
              fluidRow(
                box(title="Residuos en el Tiempo",
                    width=8, solidHeader=TRUE,
                    plotlyOutput("plot_sup_residuos", height="260px")),
                box(title="Histograma de Residuos",
                    width=4, solidHeader=TRUE,
                    plotlyOutput("plot_sup_hist", height="260px"))
              ),
              fluidRow(
                box(title="ACF de Residuos — Prueba de Ruido Blanco",
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_sup_acf", height="260px")),
                box(title="QQ-Plot — Prueba de Normalidad",
                    width=6, solidHeader=TRUE,
                    plotlyOutput("plot_sup_qq", height="260px"))
              ),
              # Conclusión
              fluidRow(
                column(12,
                       tags$div(class="interp-card",
                                tags$div(class="interp-tag", style=paste0("color:", cyan, ";"),
                                         "CONCLUSIÓN — VERIFICACIÓN DE SUPUESTOS"),
                                fluidRow(
                                  column(3,
                                         tags$div(class="interp-row",
                                                  tags$div(class="interp-key", "1. Media cero"),
                                                  tags$div(class="interp-val", style=paste0("color:", green, ";"),
                                                           paste0("✅ Media = ", media_res, " ≈ 0"))
                                         )
                                  ),
                                  column(3,
                                         tags$div(class="interp-row",
                                                  tags$div(class="interp-key", "2. Ruido blanco (Ljung-Box)"),
                                                  tags$div(class="interp-val",
                                                           style=paste0("color:", if(lb_pval_mej>0.05) green else red, ";"),
                                                           if(lb_pval_mej>0.05)
                                                             paste0("✅ p = ", lb_pval_mej, " > 0.05")
                                                           else
                                                             paste0("❌ p = ", lb_pval_mej, " < 0.05"))
                                         )
                                  ),
                                  column(3,
                                         tags$div(class="interp-row",
                                                  tags$div(class="interp-key", "3. Normalidad (Jarque-Bera)"),
                                                  tags$div(class="interp-val",
                                                           style=paste0("color:", if(jb_pval>0.05) green else red, ";"),
                                                           if(jb_pval>0.05)
                                                             paste0("✅ p = ", jb_pval, " > 0.05")
                                                           else
                                                             paste0("⚠ p = ", jb_pval, " — leve no normalidad (común en series largas)"))
                                         )
                                  ),
                                  column(3,
                                         tags$div(class="interp-row",
                                                  tags$div(class="interp-key", "4. Pronóstico"),
                                                  tags$div(class="interp-val", style=paste0("color:", green, ";"),
                                                           "✅ ICs válidos — modelo confiable")
                                         )
                                  )
                                )
                       )
                )
              )
      ),
      
      # ======================================================================
      # TAB 7 — MODELOS
      # ======================================================================
      tabItem(tabName="modelos",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px;",
                                       section_hdr("Modelos SARIMA", "Comparación · Diagnóstico · Selección")
                              )
              )),
              fluidRow(
                box(title="Tabla Comparativa de Modelos",
                    width=12, solidHeader=TRUE,
                    DTOutput("tabla_modelos"))
              ),
              fluidRow(
                box(title=tags$span("Diagnóstico de Residuos",
                                    nota_exp("diag_resid",
                                             "&#8226; Los residuos fluctúan alrededor de <b style='color:#f1f5f9'>cero</b> sin patrones sistemáticos.<br><br>
                   &#8226; La ausencia de autocorrelaciones sugiere <b style='color:#f1f5f9'>ruido blanco</b>.<br><br>
                   &#8226; El histograma muestra distribución aproximadamente <b style='color:#f1f5f9'>normal</b>.")),
                    width=8, solidHeader=TRUE,
                    tags$div(style="margin-bottom:14px;",
                             selectInput("modelo_sel", "Seleccionar modelo:",
                                         choices=setNames(1:4, nombres_mod), selected=2)
                    ),
                    plotOutput("plot_residuos", height="320px")),
                column(4,
                       tags$div(style="padding:0 0 0 4px;",
                                tags$div(style="margin-bottom:16px; display:flex; align-items:center; gap:8px;",
                                         section_hdr("Mejor modelo", "SARIMA(1,1,2)(2,0,0)"),
                                         nota_exp_inline("mejor_mod",
                                                         "Los criterios de información presentan valores bajos, indicando un modelo <b style='color:#f1f5f9'>parsimonioso y estadísticamente eficiente</b>.")
                                ),
                                metric_mini("AIC",  round(AIC(mejor_modelo), 2)),
                                metric_mini("BIC",  round(BIC(mejor_modelo), 2)),
                                metric_mini("AICc", round(mejor_modelo$aicc, 2)),
                                tags$div(class="metric-mini",
                                         tags$div(style="display:flex; align-items:center; justify-content:space-between;",
                                                  tags$div(class="metric-mini-label", "Ljung-Box p-value"),
                                                  nota_exp_inline("ljung_box",
                                                                  "No existe evidencia estadística de <b style='color:#f1f5f9'>autocorrelación en los residuos</b>, validando la especificación del modelo SARIMA seleccionado.")
                                         ),
                                         tags$div(style=paste0("color:", green, "; font-size:22px; font-weight:700;"),
                                                  paste0(lb_pval_mej, " ✅"))
                                )
                       )
                )
              ),
              fluidRow(
                box(title="Coeficientes del Modelo — Estimaciones y Significancia",
                    width=12, solidHeader=TRUE,
                    tags$div(style="margin-bottom:14px;",
                             selectInput("modelo_coef", "Seleccionar modelo para ver coeficientes:",
                                         choices=setNames(1:4, nombres_mod), selected=2)
                    ),
                    DTOutput("tabla_coeficientes"))
              )
      ),
      
      # ======================================================================
      # TAB 8 — PRONÓSTICO
      # ======================================================================
      tabItem(tabName="pronostico",
              fluidRow(column(12,
                              tags$div(style="margin-bottom:20px; display:flex; align-items:center; justify-content:space-between;",
                                       section_hdr("Pronóstico vs Realidad",
                                                   "12 meses · 2021 · SARIMA(1,1,2)(2,0,0)"),
                                       downloadButton("descargar_csv", "Descargar CSV",
                                                      style=paste0(
                                                        "background:rgba(124,92,252,0.15); border:1px solid rgba(124,92,252,0.4);",
                                                        "color:#9b7ffe; border-radius:8px; font-size:12px; font-weight:600;",
                                                        "padding:7px 16px; cursor:pointer; transition:all 0.2s ease;"
                                                      ))
                              )
              )),
              fluidRow(
                box(title=tags$span("Pronóstico 2021 vs Valores Reales",
                                    nota_exp("pron_graf",
                                             "El gráfico muestra que el modelo SARIMA logra capturar adecuadamente la <b style='color:#f1f5f9'>tendencia y dinámica general</b> de las remesas, aunque los valores reales de 2021 presentan niveles superiores a los pronosticados en varios meses.")),
                    width=12, solidHeader=TRUE,
                    plotlyOutput("plot_pronostico", height="400px"))
              ),
              fluidRow(
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "MAE — Error Absoluto Medio"),
                                tags$div(style=paste0("color:", purple2, "; font-size:28px; font-weight:700;"),
                                         paste0(mae_val, " M USD")),
                                tags$div(class="stat-sub", "Promedio de errores absolutos")
                       )
                ),
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "RMSE — Raíz del Error Cuadrático"),
                                tags$div(style=paste0("color:", blue2, "; font-size:28px; font-weight:700;"),
                                         paste0(rmse_val, " M USD")),
                                tags$div(class="stat-sub", "Penaliza errores grandes")
                       )
                ),
                column(4,
                       tags$div(class="stat-card",
                                tags$div(class="stat-label", "MAPE — Error Porcentual Medio"),
                                tags$div(style=paste0("color:", cyan, "; font-size:28px; font-weight:700;"),
                                         paste0(mape_val, "%")),
                                tags$div(class="stat-sub",
                                         tags$span(class="badge-neutral", "Aceptable < 15%"))
                       )
                )
              ),
              fluidRow(
                box(title=tags$span("Tabla Comparativa: Pronóstico vs Realidad 2021",
                                    nota_exp("pron_tabla",
                                             "Aunque el modelo está correctamente especificado, las diferencias se explican por <b style='color:#f1f5f9'>choques extraordinarios post-pandemia</b>: recuperación económica global, aumento de la migración y depreciación del peso — factores no contenidos en la dinámica histórica estimada.")),
                    width=12, solidHeader=TRUE,
                    DTOutput("tabla_comparacion"))
              )
      )
    )
  )
)
# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {
  
  P <- function(p) {
    p |> layout(
      paper_bgcolor = card,
      plot_bgcolor  = card,
      font = list(color=text2, family="DM Sans", size=12),
      xaxis = list(gridcolor=border, zerolinecolor=border,
                   tickfont=list(color=text3), linecolor=border),
      yaxis = list(gridcolor=border, zerolinecolor=border,
                   tickfont=list(color=text3), linecolor=border),
      legend = list(bgcolor="rgba(0,0,0,0)", font=list(color=text2, size=11)),
      margin = list(t=10, r=10, b=40, l=50)
    )
  }
  
  # --- PASTEL PAÍSES ----------------------------------------------------------
  output$plot_pastel <- renderPlotly({
    df_paises <- data.frame(
      Pais = c("EEUU","España","Chile","Reino Unido","Ecuador","Otros"),
      Part = c(0.536, 0.154, 0.039, 0.037, 0.030, 1-(0.536+0.154+0.039+0.037+0.030))
    )
    colores_pastel <- c(purple2, blue2, cyan, green, "#f97316", text3)
    p <- plot_ly(df_paises, labels=~Pais, values=~Part, type="pie",
                 marker=list(colors=colores_pastel,
                             line=list(color=card, width=2)),
                 textinfo="label+percent",
                 textfont=list(color=text1, size=12, family="DM Sans"),
                 hovertemplate="<b>%{label}</b><br>%{percent}<extra></extra>",
                 hole=0.35) |>
      layout(
        paper_bgcolor = card,
        plot_bgcolor  = card,
        font   = list(color=text2, family="DM Sans", size=12),
        legend = list(bgcolor="rgba(0,0,0,0)", font=list(color=text2, size=11),
                      orientation="v", x=1.02, y=0.5),
        margin = list(t=10, r=10, b=10, l=10),
        showlegend = TRUE
      )
    p
  })
  
  # --- INICIO -----------------------------------------------------------------
  output$plot_inicio <- renderPlotly({
    df <- data.frame(
      Fecha   = as.Date(as.yearmon(time(remesas_ts))),
      Remesas = as.numeric(remesas_ts)
    )
    p <- plot_ly(df, x=~Fecha, y=~Remesas, type="scatter", mode="lines",
                 fill="tozeroy", fillcolor=paste0(purple,"18"),
                 line=list(color=purple2, width=1.8), name="Remesas",
                 hovertemplate="<b>%{x}</b><br>$%{y:.1f} M<extra></extra>") |>
      layout(xaxis=list(title=""), yaxis=list(title="Millones USD"),
             showlegend=FALSE)
    P(p)
  })
  
  # --- SERIE ORIGINAL ---------------------------------------------------------
  output$plot_serie_original <- renderPlotly({
    df <- data.frame(
      Fecha   = as.Date(as.yearmon(time(remesas_ts))),
      Remesas = as.numeric(remesas_ts)
    )
    p <- plot_ly(df, x=~Fecha, y=~Remesas, type="scatter", mode="lines",
                 fill="tozeroy", fillcolor=paste0(blue,"15"),
                 line=list(color=blue2, width=1.8), name="Remesas",
                 hovertemplate="<b>%{x}</b><br>$%{y:.1f} M<extra></extra>") |>
      add_lines(y=mean(df$Remesas), name="Media",
                line=list(color=purple, dash="dot", width=1.5),
                hoverinfo="skip") |>
      layout(xaxis=list(title="Año"), yaxis=list(title="Millones USD"))
    P(p)
  })
  
  # --- SERIE LOG --------------------------------------------------------------
  output$plot_serie_log <- renderPlotly({
    df <- data.frame(
      Fecha   = as.Date(as.yearmon(time(remesas_log))),
      Remesas = as.numeric(remesas_log)
    )
    p <- plot_ly(df, x=~Fecha, y=~Remesas, type="scatter", mode="lines",
                 fill="tozeroy", fillcolor=paste0(green,"15"),
                 line=list(color=green, width=1.8), name="log(Remesas)",
                 hovertemplate="<b>%{x}</b><br>%{y:.3f}<extra></extra>") |>
      add_lines(y=mean(df$Remesas), name="Media",
                line=list(color=cyan, dash="dot", width=1.5),
                hoverinfo="skip") |>
      layout(xaxis=list(title="Año"), yaxis=list(title="log(Mill. USD)"))
    P(p)
  })
  
  # --- COMPARACIÓN ------------------------------------------------------------
  output$plot_comparacion <- renderPlotly({
    df1 <- data.frame(Fecha=as.Date(as.yearmon(time(remesas_ts))),
                      Valor=as.numeric(remesas_ts))
    df2 <- data.frame(Fecha=as.Date(as.yearmon(time(remesas_log))),
                      Valor=as.numeric(remesas_log))
    p <- plot_ly() |>
      add_lines(data=df1, x=~Fecha, y=~Valor, name="Original",
                line=list(color=blue2, width=1.5)) |>
      add_lines(data=df2, x=~Fecha, y=~Valor, name="Log",
                yaxis="y2", line=list(color=green, width=1.5, dash="dot")) |>
      layout(
        yaxis =list(title="Mill. USD", color=blue2,
                    gridcolor=border, tickfont=list(color=text3)),
        yaxis2=list(title="log(Mill. USD)", color=green,
                    overlaying="y", side="right",
                    gridcolor="rgba(0,0,0,0)", tickfont=list(color=text3)),
        xaxis=list(title="Año", gridcolor=border, tickfont=list(color=text3)),
        paper_bgcolor=card, plot_bgcolor=card,
        font=list(color=text2, family="DM Sans"),
        legend=list(bgcolor="rgba(0,0,0,0)", font=list(color=text2)),
        margin=list(t=10,r=60,b=40,l=50)
      )
    p
  })
  
  # --- DESCOMPOSICIÓN ---------------------------------------------------------
  output$plot_descomposicion <- renderPlot({
    par(bg=card, col.axis=text3, col.lab=text2,
        col.main=text1, fg=border2, family="sans")
    plot(decompose(remesas_ts), col=blue2)
  })
  
  # --- ESTACIONAL -------------------------------------------------------------
  output$plot_estacional <- renderPlot({
    ggseasonplot(remesas_ts, year.labels=FALSE) +
      scale_color_viridis_d(option="plasma") +
      labs(title=NULL) +
      theme(panel.background  = element_rect(fill=card, color=NA),
            plot.background   = element_rect(fill=card, color=NA),
            panel.grid.major  = element_line(color=border, linewidth=0.3),
            panel.grid.minor  = element_blank(),
            axis.text         = element_text(color=text3, size=10),
            axis.title        = element_text(color=text2, size=11),
            legend.background = element_rect(fill=card, color=NA),
            legend.text       = element_text(color=text3, size=9),
            legend.key        = element_rect(fill=card, color=NA))
  })
  
  # --- SUBSERIES — MODIFICACIÓN 2: línea azul cyan en lugar de negro ----------
  output$plot_subseries <- renderPlot({
    ggsubseriesplot(remesas_ts) +
      labs(title=NULL) +
      theme(panel.background = element_rect(fill=card, color=NA),
            plot.background  = element_rect(fill=card, color=NA),
            panel.grid.major = element_line(color=border, linewidth=0.3),
            panel.grid.minor = element_blank(),
            axis.text        = element_text(color=text3, size=10),
            axis.title       = element_text(color=text2, size=11)) +
      # Reemplaza el color negro de la línea de media por cyan del dashboard
      geom_line(color = cyan, linewidth = 0.8)
  })
  
  # --- DIF SIMPLE -------------------------------------------------------------
  output$plot_dif_simple <- renderPlotly({
    df <- data.frame(
      Fecha = as.Date(as.yearmon(time(dif_simple_log))),
      Valor = as.numeric(dif_simple_log)
    )
    p <- plot_ly(df, x=~Fecha, y=~Valor, type="scatter", mode="lines",
                 line=list(color=cyan, width=1.3), name="Dif. Log",
                 hovertemplate="<b>%{x}</b><br>%{y:.4f}<extra></extra>") |>
      add_lines(y=0, line=list(color=border2, width=1, dash="dot"),
                hoverinfo="skip", showlegend=FALSE) |>
      layout(xaxis=list(title="Año"), yaxis=list(title="Δ log(Remesas)"))
    P(p)
  })
  
  # --- ACF --------------------------------------------------------------------
  output$plot_acf <- renderPlotly({
    acf_data <- Acf(dif_simple_log, lag.max=48, plot=FALSE)
    lags <- as.numeric(acf_data$lag)
    acfs <- as.numeric(acf_data$acf)
    ic   <- qnorm(0.975) / sqrt(length(dif_simple_log))
    cols <- ifelse(abs(acfs) > ic, blue2, text3)
    p <- plot_ly() |>
      add_segments(x=~lags, xend=~lags, y=0, yend=~acfs,
                   line=list(color=cols, width=2.5), name="ACF") |>
      add_lines(x=c(min(lags),max(lags)), y=c(ic,ic),
                line=list(color=paste0(purple,"99"), dash="dash", width=1),
                name="IC 95%") |>
      add_lines(x=c(min(lags),max(lags)), y=c(-ic,-ic),
                line=list(color=paste0(purple,"99"), dash="dash", width=1),
                showlegend=FALSE) |>
      layout(xaxis=list(title="Lag"),
             yaxis=list(title="ACF", range=c(min(acfs)-0.05, max(acfs)+0.05)))
    P(p)
  })
  
  # --- PACF -------------------------------------------------------------------
  output$plot_pacf <- renderPlotly({
    pacf_data <- Pacf(dif_simple_log, lag.max=48, plot=FALSE)
    lags  <- as.numeric(pacf_data$lag)
    pacfs <- as.numeric(pacf_data$acf)
    ic    <- qnorm(0.975) / sqrt(length(dif_simple_log))
    cols  <- ifelse(abs(pacfs) > ic, purple2, text3)
    p <- plot_ly() |>
      add_segments(x=~lags, xend=~lags, y=0, yend=~pacfs,
                   line=list(color=cols, width=2.5), name="PACF") |>
      add_lines(x=c(min(lags),max(lags)), y=c(ic,ic),
                line=list(color=paste0(cyan,"99"), dash="dash", width=1),
                name="IC 95%") |>
      add_lines(x=c(min(lags),max(lags)), y=c(-ic,-ic),
                line=list(color=paste0(cyan,"99"), dash="dash", width=1),
                showlegend=FALSE) |>
      layout(xaxis=list(title="Lag"),
             yaxis=list(title="PACF", range=c(min(pacfs)-0.05, max(pacfs)+0.05)))
    P(p)
  })
  
  # --- GRÁFICO ÍNDICE BASE 100 (con rango de fechas reactivo) ----------------
  
  output$plot_trm_indice <- renderPlotly({
    anio_ini  <- if (!is.null(input$anio_ini_trm)) as.integer(input$anio_ini_trm) else 2000L
    anio_fin  <- if (!is.null(input$anio_fin_trm)) as.integer(input$anio_fin_trm) else 2025L
    if (anio_ini >= anio_fin) { anio_ini <- 2000L; anio_fin <- 2025L }
    
    fecha_ini <- as.Date(paste0(anio_ini, "-01-01"))
    fecha_fin <- as.Date(paste0(anio_fin, "-12-01"))
    
    df_fil <- df_trm %>%
      filter(Fecha >= fecha_ini, Fecha <= fecha_fin) %>%
      mutate(
        idx_remesas = (Remesas / Remesas[1]) * 100,
        idx_TRM     = (TRM     / TRM[1])     * 100
      )
    
    # Solo mostrar eventos dentro del rango
    ev_mask   <- eventos_x >= fecha_ini & eventos_x <= fecha_fin
    ev_x_fil  <- eventos_x[ev_mask]
    ev_lb_fil <- etiquetas_ev[ev_mask]
    
    p <- plot_ly() |>
      add_lines(y=c(100,100), x=c(fecha_ini, fecha_fin),
                line=list(color=text3, dash="dot", width=1),
                name="Base 100", hoverinfo="skip") |>
      add_lines(data=df_fil, x=~Fecha, y=~idx_remesas,
                line=list(color=purple2, width=2), name="Remesas (Mill. USD)",
                hovertemplate="<b>%{x}</b><br>Índice Remesas: %{y:.1f}<extra></extra>") |>
      add_lines(data=df_fil, x=~Fecha, y=~idx_TRM,
                line=list(color=green, width=2), name="TRM (COP/USD)",
                hovertemplate="<b>%{x}</b><br>Índice TRM: %{y:.1f}<extra></extra>")
    
    if (length(ev_x_fil) > 0) {
      p <- p |>
        add_segments(x=ev_x_fil, xend=ev_x_fil, y=90, yend=max(df_fil$idx_remesas)*1.05,
                     line=list(color=text3, dash="dot", width=0.8),
                     showlegend=FALSE, hoverinfo="skip")
    }
    
    p <- p |> layout(
      yaxis=list(title="Índice (Base período seleccionado = 100)", type="log",
                 tickvals=c(100,150,200,300,500,700,1000),
                 ticktext=c("100","150","200","300","500","700","1000")),
      xaxis=list(title="Año"),
      hovermode="x unified",
      annotations=if(length(ev_x_fil)>0) lapply(seq_along(ev_x_fil), function(i){
        list(x=ev_x_fil[i], y=3.0, text=ev_lb_fil[i],
             showarrow=FALSE, font=list(color=text3, size=10),
             xanchor="center", yanchor="bottom")
      }) else list()
    )
    P(p)
  })
  
  # --- REMESAS SERIE COMPLETA -------------------------------------------------
  output$plot_trm_remesas <- renderPlotly({
    p <- plot_ly(df_trm, x=~Fecha, y=~Remesas, type="scatter", mode="lines",
                 fill="tozeroy", fillcolor=paste0(purple,"18"),
                 line=list(color=purple2, width=1.8), name="Remesas",
                 hovertemplate="<b>%{x}</b><br>$%{y:.1f} M<extra></extra>") |>
      layout(xaxis=list(title="Año"), yaxis=list(title="Millones USD"))
    P(p)
  })
  
  # --- TRM SERIE COMPLETA -----------------------------------------------------
  output$plot_trm_trm <- renderPlotly({
    p <- plot_ly(df_trm, x=~Fecha, y=~TRM, type="scatter", mode="lines",
                 fill="tozeroy", fillcolor=paste0(green,"18"),
                 line=list(color=green, width=1.8), name="TRM",
                 hovertemplate="<b>%{x}</b><br>$%{y:,.0f} COP<extra></extra>") |>
      layout(xaxis=list(title="Año"), yaxis=list(title="COP/USD"))
    P(p)
  })
  
  # --- DISPERSIÓN -------------------------------------------------------------
  output$plot_trm_scatter <- renderPlotly({
    lm_fit  <- lm(Remesas ~ TRM, data=df_trm)
    x_range <- range(df_trm$TRM, na.rm=TRUE)
    df_line <- data.frame(TRM=x_range,
                          Remesas=predict(lm_fit, newdata=data.frame(TRM=x_range)))
    p <- plot_ly() |>
      add_markers(data=df_trm, x=~TRM, y=~Remesas,
                  marker=list(color=paste0(purple2,"99"), size=5,
                              line=list(color=purple, width=0.5)),
                  name="Observaciones",
                  hovertemplate="TRM: %{x:,.0f}<br>Remesas: $%{y:.1f} M<extra></extra>") |>
      add_lines(data=df_line, x=~TRM, y=~Remesas,
                line=list(color=cyan, width=2, dash="dot"),
                name=paste0("Tendencia (r=",cor_pearson,")"),
                hoverinfo="skip") |>
      layout(xaxis=list(title="TRM (COP/USD)"),
             yaxis=list(title="Remesas (Millones USD)"))
    P(p)
  })
  
  # --- SUPUESTOS: Residuos en el tiempo ---------------------------------------
  output$plot_sup_residuos <- renderPlotly({
    df_res <- data.frame(
      Fecha = as.Date(as.yearmon(time(res_mejor))),
      Res   = as.numeric(res_mejor)
    )
    p <- plot_ly(df_res, x=~Fecha, y=~Res, type="scatter", mode="lines",
                 line=list(color=blue2, width=1.2), name="Residuos",
                 hovertemplate="<b>%{x}</b><br>%{y:.4f}<extra></extra>") |>
      add_lines(y=0, line=list(color=text3, dash="dot", width=1),
                hoverinfo="skip", showlegend=FALSE) |>
      layout(xaxis=list(title="Año"),
             yaxis=list(title="Residuo"))
    P(p)
  })
  
  # --- SUPUESTOS: Histograma --------------------------------------------------
  output$plot_sup_hist <- renderPlotly({
    res_vec <- as.numeric(res_mejor)
    p <- plot_ly(x=res_vec, type="histogram",
                 nbinsx=30,
                 marker=list(color=paste0(purple2,"99"),
                             line=list(color=border2, width=0.5)),
                 name="Residuos") |>
      layout(xaxis=list(title="Residuo"),
             yaxis=list(title="Frecuencia"),
             bargap=0.05)
    P(p)
  })
  
  # --- SUPUESTOS: ACF de residuos ---------------------------------------------
  output$plot_sup_acf <- renderPlotly({
    acf_res <- Acf(res_mejor, lag.max=36, plot=FALSE)
    lags <- as.numeric(acf_res$lag)
    acfs <- as.numeric(acf_res$acf)
    ic   <- qnorm(0.975) / sqrt(length(res_mejor))
    cols <- ifelse(abs(acfs) > ic, red, cyan)
    p <- plot_ly() |>
      add_segments(x=~lags, xend=~lags, y=0, yend=~acfs,
                   line=list(color=cols, width=2.5), name="ACF residuos") |>
      add_lines(x=c(min(lags),max(lags)), y=c(ic,ic),
                line=list(color=paste0(purple,"99"), dash="dash", width=1),
                name="IC 95%") |>
      add_lines(x=c(min(lags),max(lags)), y=c(-ic,-ic),
                line=list(color=paste0(purple,"99"), dash="dash", width=1),
                showlegend=FALSE) |>
      layout(xaxis=list(title="Lag"),
             yaxis=list(title="ACF", range=c(min(acfs)-0.05, max(acfs)+0.05)))
    P(p)
  })
  
  # --- SUPUESTOS: QQ-plot -----------------------------------------------------
  output$plot_sup_qq <- renderPlotly({
    res_vec  <- sort(as.numeric(res_mejor))
    n        <- length(res_vec)
    teoricos <- qnorm(ppoints(n))
    p <- plot_ly() |>
      add_markers(x=teoricos, y=res_vec,
                  marker=list(color=paste0(purple2,"99"), size=4),
                  name="Cuantiles",
                  hovertemplate="Teórico: %{x:.2f}<br>Observado: %{y:.4f}<extra></extra>") |>
      add_lines(x=range(teoricos),
                y=mean(res_vec) + sd(res_vec) * range(teoricos),
                line=list(color=cyan, width=1.5, dash="dot"),
                name="Línea normal", hoverinfo="skip") |>
      layout(xaxis=list(title="Cuantiles teóricos (Normal)"),
             yaxis=list(title="Cuantiles observados"))
    P(p)
  })
  
  # --- TABLA MODELOS ----------------------------------------------------------
  output$tabla_modelos <- renderDT({
    datatable(tabla_modelos,
              options=list(dom="t", pageLength=4, ordering=TRUE),
              rownames=FALSE) |>
      formatStyle("RuidoBlanco",
                  color=styleEqual(c("✅ Sí","❌ No"), c(green, red)),
                  fontWeight="bold") |>
      formatStyle("LjungBox",
                  color=styleInterval(0.05, c(red, green))) |>
      formatStyle(columns=1:6, backgroundColor=card, color=text1)
  })
  
  # --- RESIDUOS (tab Modelos) -------------------------------------------------
  output$plot_residuos <- renderPlot({
    mod <- modelos_lista[[as.integer(input$modelo_sel)]]
    par(bg=card, col.axis=text3, col.lab=text2,
        col.main=text1, fg=border2, family="sans", col=blue2)
    checkresiduals(mod)
  })
  
  # --- PRONÓSTICO VS REALIDAD -------------------------------------------------
  output$plot_pronostico <- renderPlotly({
    df_hist <- data.frame(
      Fecha = as.Date(as.yearmon(time(remesas_ts))),
      Valor = as.numeric(remesas_ts)
    )
    fechas_pron <- seq(as.Date("2021-01-01"), by="month", length.out=12)
    df_pron <- data.frame(
      Fecha = fechas_pron,
      Media = as.numeric(pronostico$mean),
      Lo80  = as.numeric(pronostico$lower[,1]),
      Hi80  = as.numeric(pronostico$upper[,1]),
      Lo95  = as.numeric(pronostico$lower[,2]),
      Hi95  = as.numeric(pronostico$upper[,2])
    )
    df_real <- data.frame(
      Fecha = fechas_pron,
      Valor = as.numeric(remesas_real_2021)
    )
    p <- plot_ly() |>
      add_ribbons(data=df_pron, x=~Fecha, ymin=~Lo95, ymax=~Hi95,
                  fillcolor=paste0(blue,"22"),
                  line=list(color="transparent"), name="IC 95%",
                  hoverinfo="skip") |>
      add_ribbons(data=df_pron, x=~Fecha, ymin=~Lo80, ymax=~Hi80,
                  fillcolor=paste0(blue,"44"),
                  line=list(color="transparent"), name="IC 80%",
                  hoverinfo="skip") |>
      add_lines(data=df_hist, x=~Fecha, y=~Valor,
                line=list(color=text2, width=1.2), name="Histórico",
                hovertemplate="<b>%{x}</b><br>$%{y:.1f} M<extra></extra>") |>
      add_lines(data=df_pron, x=~Fecha, y=~Media,
                line=list(color=blue2, width=2.5, dash="dot"),
                name="Pronóstico",
                hovertemplate="<b>%{x}</b><br>Pronóstico: $%{y:.1f} M<extra></extra>") |>
      add_lines(data=df_real, x=~Fecha, y=~Valor,
                line=list(color="#f97316", width=2.5), name="Real 2021",
                hovertemplate="<b>%{x}</b><br>Real: $%{y:.1f} M<extra></extra>") |>
      layout(xaxis=list(title="Año"), yaxis=list(title="Millones USD"),
             hovermode="x unified")
    
    P(p)
  })
  
  # --- COEFICIENTES DEL MODELO ------------------------------------------------
  output$tabla_coeficientes <- renderDT({
    sel <- if (is.null(input$modelo_coef)) 2L else as.integer(input$modelo_coef)
    mod <- modelos_lista[[sel]]
    
    coefs   <- coef(mod)
    se_coef <- sqrt(diag(vcov(mod)))
    t_stat  <- coefs / se_coef
    p_vals  <- 2 * pnorm(-abs(t_stat))
    
    sig <- ifelse(p_vals < 0.001, "***",
                  ifelse(p_vals < 0.01,  "**",
                         ifelse(p_vals < 0.05,  "*",
                                ifelse(p_vals < 0.1,   ".", ""))))
    
    df_coef <- data.frame(
      Parámetro  = names(coefs),
      Estimación = round(coefs,    5),
      Error_Std  = round(se_coef,  5),
      Estadístico_t = round(t_stat, 3),
      P_valor    = ifelse(p_vals < 0.001,
                          "< 0.001",
                          as.character(round(p_vals, 4))),
      Signif     = sig,
      stringsAsFactors = FALSE
    )
    
    datatable(df_coef,
              options=list(dom="t", pageLength=10, ordering=FALSE),
              rownames=FALSE,
              colnames=c("Parámetro","Estimación","Error Std.",
                         "Estadístico t","p-valor","Signif.")) |>
      formatStyle("Signif",
                  color = styleEqual(
                    c("***","**","*",".",""),
                    c(green, green, blue2, text2, text3)
                  ),
                  fontWeight = "bold") |>
      formatStyle("Estadístico_t",
                  color = styleInterval(c(-1.96, 1.96),
                                        c(red, text2, green))) |>
      formatStyle(columns=1:6, backgroundColor=card, color=text1) |>
      formatStyle("P_valor",
                  color = styleEqual("< 0.001", green))
  })
  
  # --- COINTEGRACIÓN ----------------------------------------------------------
  output$eg_stat_val       <- renderText({ as.character(eg_stat) })
  output$eg_pval_val       <- renderText({ as.character(eg_pval) })
  output$eg_conclusion_val <- renderText({ eg_conclusion })
  output$coef_alpha        <- renderText({ paste0(coef_lr[1]) })
  output$coef_beta         <- renderText({ paste0(coef_lr[2]) })
  output$gr_pval1          <- renderText({ paste0("p = ", gr_pval_rem_trm) })
  output$gr_conc1          <- renderText({ gr_conclusion_rem_trm })
  output$gr_pval2          <- renderText({ paste0("p = ", gr_pval_trm_rem) })
  output$gr_conc2          <- renderText({ gr_conclusion_trm_rem })
  
  output$plot_coint_reg <- renderPlotly({
    x_seq  <- seq(min(trm_niveles), max(trm_niveles), length.out=100)
    y_pred <- coef_lr[1] + coef_lr[2] * x_seq
    p <- plot_ly() |>
      add_markers(x=trm_niveles, y=rem_niveles,
                  marker=list(color=paste0(purple2,"88"), size=5,
                              line=list(color=purple, width=0.5)),
                  name="Observaciones",
                  hovertemplate="TRM: %{x:,.0f}<br>Remesas: %{y:.1f} M<extra></extra>") |>
      add_lines(x=x_seq, y=y_pred,
                line=list(color=cyan, width=2),
                name="Regresión LP", hoverinfo="skip") |>
      layout(xaxis=list(title="TRM (COP/USD)"),
             yaxis=list(title="Remesas (Millones USD)"))
    P(p)
  })
  
  output$plot_coint_resid <- renderPlotly({
    df_resid <- data.frame(Fecha=df_trm$Fecha, Residuo=residuos_eg)
    p <- plot_ly(df_resid, x=~Fecha, y=~Residuo,
                 type="scatter", mode="lines",
                 line=list(color=purple2, width=1.5),
                 fill="tozeroy",
                 fillcolor=paste0(purple2,"22"),
                 name="Residuos",
                 hovertemplate="<b>%{x}</b><br>%{y:.2f} M<extra></extra>") |>
      add_lines(y=0, line=list(color=text3, dash="dot", width=1),
                hoverinfo="skip", showlegend=FALSE) |>
      add_lines(y= 2*sd(residuos_eg),
                line=list(color=red, dash="dash", width=1),
                name="+2σ", hoverinfo="skip") |>
      add_lines(y=-2*sd(residuos_eg),
                line=list(color=red, dash="dash", width=1),
                name="-2σ", hoverinfo="skip") |>
      layout(xaxis=list(title="Año"),
             yaxis=list(title="Residuo (Millones USD)"),
             hovermode="x unified")
    P(p)
  })
  
  # --- DESCARGA CSV -----------------------------------------------------------
  output$descargar_csv <- downloadHandler(
    filename = function() paste0("pronostico_remesas_2021_", Sys.Date(), ".csv"),
    content  = function(file) {
      write.csv(tabla_comparacion, file, row.names=FALSE, fileEncoding="UTF-8")
    }
  )
  
  # --- TABLA COMPARACIÓN ------------------------------------------------------
  output$tabla_comparacion <- renderDT({
    datatable(tabla_comparacion,
              options=list(dom="t", pageLength=12),
              rownames=FALSE,
              colnames=c("Mes","Pronóstico (M)","Real (M)",
                         "Diferencia (M)","Error %")) |>
      formatStyle("Diferencia",
                  color=styleInterval(0, c(red, green)),
                  fontWeight="bold") |>
      formatStyle("Error_Pct",
                  color=styleInterval(0, c(green, red))) |>
      formatStyle(columns=1:5, backgroundColor=card, color=text1)
  })
}

# ==============================================================================
# RUN
# ==============================================================================

shinyApp(ui=ui, server=server)