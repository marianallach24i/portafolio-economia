# Factores asociados al número de hijos en jóvenes de 18 a 25 años

Trabajo de Econometría I elaborado en coautoría con Rodrigo Rivera de la Asunción, sobre los factores asociados al número de hijos esperados por jóvenes colombianos entre 18 y 25 años.

📄 **Póster completo:** ver archivo PDF en esta carpeta.

## Contexto

La disminución de la natalidad en Colombia tiene un trasfondo socioeconómico y psicológico que incide en las expectativas de los jóvenes sobre el número de hijos que desean tener.

## Datos y metodología

- Encuesta propia aplicada a una muestra de 98 jóvenes.
- Modelo de regresión **quasi-Poisson** para estimar el número de hijos esperados como variable dependiente de conteo.
- Variables explicativas: salario esperado, cursos aprobados, estrato socioeconómico e incentivos económicos para no tener hijos.
- Diagnósticos de multicolinealidad (Factores de Inflación de la Varianza, VIF) y de dispersión (residuos de Pearson) para validar la robustez del modelo.

## Resultados principales

- Bajo nivel de multicolinealidad entre las variables (todos los VIF menores a 2).
- El modelo mostró subdispersión leve, corregida adecuadamente mediante el ajuste quasi-Poisson.
- Los incentivos económicos que desincentivan la natalidad resultaron ser el factor más relevante: a mayor exposición a estos incentivos, menor el número de hijos esperados.
- El salario esperado y el número de cursos aprobados no mostraron una relación estadísticamente clara con la expectativa de tener hijos.

## Herramientas

R (modelo de regresión quasi-Poisson y diagnósticos estadísticos).
