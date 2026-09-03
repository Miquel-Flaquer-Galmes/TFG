# Grafs metabòlics i cadenes de Markov a partir d'abundàncies de KOs

El projecte construeix **grafs dirigits de reaccions metabòliques** a partir de matrius d'abundància d'ortòlegs KEGG (KOs), els analitza com a **cadenes de Markov** i en compara l'estructura entre mostres mitjançant un **kernel de propagació/transició**, resumit en dendrogrames i heatmaps.

L'anàlisi s'aplica a tres cohorts independents:

| Cohort | Grups | Matriu d'abundància | Metadades |
|--------|-------|---------------------|-----------|
| Grups d'edat | Adults / Infants / Elders | `GutMicrobiota.tsv` | (deduïdes de l'ID de mostra) |
| Dietes | Korean / Western | `KoreanWestern.csv` | `Results.csv` |
| Regnes taxonòmics | Animal, Planta, Fong, Arqueu, Bacteri… | `matriu_KO.csv` | `metadata.csv` |

## Idea general

A partir d'una matriu de KOs (files = KOs, columnes = mostres) el flux és:

1. **KO → reaccions.** Cada KO es mapeja als seus enzims (EC) i cada EC a les reaccions que catalitza, fent servir els diccionaris JSON.
2. **Graf de reaccions.** Es posa una aresta `R1 → R2` quan un producte de `R1` és substrat de `R2`, excloent compostos ubics (aigua, ATP, ADP, NAD(P)(H), fosfat…). Les reaccions reversibles es desdoblen en un node `_rev`.
3. **Cadena de Markov.** El graf es converteix en una matriu de transició estocàstica per files (ponderada per l'abundància del node de destí) i s'analitza: classificació d'estats, ergodicitat, distribucions estacionàries i anàlisi d'absorció.
4. **Comparació entre mostres.** Un kernel de propagació calcula la similitud entre els grafs de mostres diferents; la matriu de similitud resultant alimenta dendrogrames i heatmaps.

## Estructura del repositori

### Scripts (`.qmd`)

| Fitxer | Funció |
|--------|--------|
| `Construcció_K0s_web.qmd` | Descarrega en línia els KOs de cada organisme des de l'API de KEGG i construeix la matriu i les metadades de la cohort de regnes. |
| `Lectura_dades_possibilitat_excloure.qmd` | Construeix els grafs de reaccions per mostra a partir d'una matriu d'abundància, amb filtre de prevalença i exclusió de compostos ubics. |
| `Cadena_Markov_unificat.qmd` | Flux unificat: estimació estadística de l'esperança per KO (Zero-Inflated Gamma / Gamma), construcció del graf, matriu de transició i anàlisi completa de la cadena de Markov. |
| `kernel_transicio.qmd` | Kernel de propagació/transició entre grafs metabòlics i heatmaps de similitud (per a la cohort d'edat). |
| `kernel_trans_Benja.qmd` | Dendrogrames circulars (fan) i heatmaps del kernel per grups d'edat. |
| `dendrograma_KW.qmd` | Pipeline de dendrogrames de la cohort de dietes (Korean/Western). |
| `dendrograma_regnes.qmd` | Pipeline de dendrogrames de la cohort de regnes taxonòmics. |

### Dades d'entrada

| Carpeta | Contingut |
|---------|-----------|
| `Dades/` | Matrius d'abundància de KOs de les tres cohorts i fitxers de metadades (`metadata.csv`, `Results.csv`). |
| `Diccionaris/` | Diccionaris de mapeig `KO → EC → reaccions`: `enzymes_KO.json`, `enzymes_lite_new.json`, `reactions_lite_pw.json`. |

### Grafs metabòlics (`.graphml`)

Els grafs de reaccions de totes les mostres, en format `.graphml`, estan disponibles en línia (fitxers massa voluminosos per al repositori):

**Enllaç a Drive:** (https://drive.google.com/drive/folders/1UTT0xBHlQ-aPiDn_2mDtLfrqzlHIwQzo?usp=sharing)

### Resultats

| Carpeta / fitxer | Contingut |
|------------------|-----------|
| `GutMicrobiota/dendrogrames/` | Dendrogrames circulars del kernel de transició per grups d'edat (5 iteracions i max_diam × Adults/Infants/Elders). |
| `Dendogrames_KW/` | Dendrogrames de la cohort de dietes: 30 rèpliques × 3 configuracions d'iteracions, distàncies entre arbres i diferència simètrica. |
| `Dendogrames_regnes/` | El mateix per a la cohort de regnes taxonòmics. |
| `Heatmaps/` | Heatmaps de similitud |


## El kernel de propagació

La similitud entre dos grafs es calcula amb `Prop_att_ker()`: es parteix del vector d'abundàncies de reaccions de cada graf i, a cada iteració, es propaga pel graf multiplicant per la matriu de propagació (adjacència ponderada per abundància i normalitzada per files). La similitud és la **suma dels cosinus** entre els dos vectors al llarg de les iteracions, i després es normalitza. Es proven tres configuracions d'iteracions: **1**, **5** i **max_diam** (el diàmetre màxim del conjunt de grafs).

Per avaluar l'estabilitat, cada cohort genera **30 rèpliques** simulant les abundàncies positives amb una Gamma per KO, i es comparen els arbres resultants amb distàncies de Robinson-Foulds (RF, wRF), *path difference* i nodal, així com la diferència simètrica respecte a un arbre de referència (abundàncies uniformes).

## Requisits

Projecte en **R** (obre `TFG_rec.Rproj`). Paquets principals:

```r
install.packages(c(
  "tidyverse", "igraph", "Matrix", "jsonlite",
  "fitdistrplus", "markovchain", "goftest",
  "lsa", "dendextend", "ape", "phangorn",
  "circlize", "pheatmap"
))
# ComplexHeatmap (Bioconductor):
# install.packages("BiocManager"); BiocManager::install("ComplexHeatmap")
```

`Construcció_K0s_web.qmd` requereix connexió a internet (API de KEGG, `https://rest.kegg.jp`).

## Ús

L'ordre habitual d'execució és:

1. Obtenir les matrius d'abundància (ja incloses a `Dades/`, o regenerar la de regnes amb `Construcció_K0s_web.qmd`).
2. Construir els grafs metabòlics amb `Lectura_dades_possibilitat_excloure.qmd`.
3. Executar l'anàlisi desitjada: cadena de Markov (`Cadena_Markov_unificat.qmd`) o comparació entre mostres (`kernel_transicio.qmd`, `kernel_trans_Benja.qmd`, `dendrograma_KW.qmd`, `dendrograma_regnes.qmd`).

> **Nota sobre rutes.** Els scripts fan referència a rutes com `Mostres/`, `Regnes/` i `Dades_Pere/`. Segons on hagis desat els fitxers d'entrada d'aquest repositori (`Dades/`, `Diccionaris/`), ajusta les variables de ruta i nom del principi de cada `.qmd` (`ABUND_FILE`, `META_FILE`, `KO_EC_FILE`, `EC_RXN_FILE`, `RXN_FILE`) perquè apuntin a les carpetes correctes.
