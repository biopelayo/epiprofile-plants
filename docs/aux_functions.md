# Funciones auxiliares y de infraestructura de EpiProfile-PLANTS

Listado completo de las **55 funciones auxiliares** que componen la infraestructura de `EpiProfile-PLANTS`, con líneas de código y función asignada. Este fichero se referencia desde el Apéndice A2 del Capítulo 2 de la tesis doctoral de Pelayo González de Lena («Cuantificación reproducible de modificaciones postraduccionales de histonas en *Arabidopsis thaliana* a lo largo del desarrollo», Universidad de Oviedo).

Total: **55 archivos · 19 794 líneas de código** (funciones puras, sin contar módulos de cuantificación por péptido ni orquestador).

## Tabla

| Archivo | Líneas | Función |
|----------------|---------------|----------------------------------------|
| EpiProfile.m | 110 | Punto de entrada principal; GUI |
| DrawISOProfile0.m | 89 | Orquestador original (organismo mamífero) |
| DrawISOProfile1.m | 331 | Orquestador adaptado (*A. thaliana*) |
| ReadInput.m | 75 | Lectura de parámetros de entrada |
| Raw2MS.m | 55 | Interfaz de conversión raw-to-mzML |
| GetIPV.m | 12.003 | Tabla de patrones isotópicos (11.999 filas de datos) |
| GetMS1ScanNo.m | 227 | Extracción de *scan numbers* MS1 |
| GetMS2ScanNo.m | 310 | Extracción de *scan numbers* MS2 |
| GetPSM.m | 232 | Asignación péptido-espectro (*peptide-spectrum match*) |
| GetPSM2.m | 242 | Asignación PSM extendida |
| GetProfiles.m | 93 | Extracción de perfiles cromatográficos |
| GetTopBottom.m | 102 | Detección de picos con *splitting* (umbrales 1,5 %/20 %/5 %) |
| GetTopBottom11.m | 104 | Detección de picos sin *splitting* |
| GetLocal.m | 69 | Búsqueda de máximos/mínimos locales |
| GetBenchmark.m | 23 | Evaluación de rendimiento |
| GetMods.m | 184 | Gestión de modificaciones PTM |
| Getaamass.m | 30 | Tabla de masas de aminoácidos |
| JudgeLocalmaxmin.m | 99 | Clasificación de extremos locales |
| OutputFigures.m | 402 | Generación de figuras de diagnóstico |
| OutputSinglePTMs.m | 229 | Exportación de PTMs individuales |
| OutputTogether.m | 164 | Exportación consolidada |
| calculate_pepmz.m | 46 | Cálculo de *m/z* teórico de péptidos |
| check_layout.m | 19 | Verificación de disposición de datos |
| check_otherparas.m | 21 | Verificación de parámetros adicionales |
| check_ref.m | 26 | Verificación de RT de referencia |
| draw_layout.m | 235 | Renderizado de disposición gráfica |
| find_pair.m | 35 | Deconvolución de pares coeluyentes |
| find_pair_new.m | 43 | Deconvolución de pares (versión mejorada) |
| find_triple.m | 104 | Deconvolución de tripletes coeluyentes (contribución original) |
| get_area.m | 212 | Integración de áreas cromatográficas (incluye `judgeOverlap1`) |
| get_histone0.m | 195 | Inicialización de estructura de histona (nivel 0) |
| get_histone1.m | 62 | Inicialización de histona (nivel 1) |
| get_histone2.m | 383 | Inicialización de histona (nivel 2) |
| get_histone4.m | 431 | Inicialización de histona (nivel 4) |
| get_histone6.m | 680 | Inicialización de histona (nivel 6) |
| get_histone10.m | 62 | Inicialización de histona (nivel 10) |
| get_histone11.m | 62 | Inicialización de histona (nivel 11) |
| get_histone12.m | 62 | Inicialización de histona (nivel 12) |
| get_histone13.m | 62 | Inicialización de histona (nivel 13) |
| get_histone22.m | 369 | Inicialización de histona (nivel 22) |
| get_key_ions.m | 63 | Extracción de iones fragmento diagnósticos |
| get_key_ions1.m | 50 | Iones fragmento (variante 1) |
| get_key_ions2.m | 65 | Iones fragmento (variante 2) |
| get_main_ch.m | 14 | Determinación del estado de carga principal |
| get_mod_mass.m | 30 | Masas de modificaciones PTM |
| get_mod_postype.m | 28 | Tipo y posición de modificaciones |
| get_rts.m | 49 | Extracción de tiempos de retención |
| get_rts0.m | 51 | Extracción de RT (variante inicial) |
| get_rts2.m | 275 | Validación cruzada MS1-MS2 por similitud coseno |
| get_rts22.m | 275 | Validación cruzada (variante) |
| get_theo_mz.m | 62 | Cálculo de *m/z* teórico |
| init_histone0.m | 245 | Catálogo de referencia (18 péptidos, \~66 hPF) |
| output_histone.m | 100 | Exportación de resultados por histona |
| output_histone2.m | 129 | Exportación de resultados (formato extendido) |
| relocateD.m | 18 | Utilidad de reubicación de directorios |


## Cómo citar

Este catálogo forma parte del bundle EpiProfile-PLANTS ([`biopelayo/epiprofile-plants`](https://github.com/biopelayo/epiprofile-plants)). Ver el CITATION.cff en la raíz del repositorio para el formato de cita canónico.

## Contexto

Los módulos de cuantificación por péptido (11 H3 core, 8 controles H3, 9 H4, 11 H2A, 8 H2B, 4 H1) se documentan en la Tabla A2.1 del Capítulo 2 de la tesis. Los módulos snapshot y legacy se documentan en las Tablas A2.1.G y A2.1.H respectivamente.
