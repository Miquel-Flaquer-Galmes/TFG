```{r}
## ============================================================
##  build_metabolic_graphs.R
##  Builds per-patient directed metabolic reaction graphs from
##  KO abundance data using proportional distribution to avoid
##  functional inflation.
## ============================================================

suppressPackageStartupMessages({
  library(igraph)
  library(jsonlite)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
})

## ── 0. CONFIGURATION ─────────────────────────────────────────
ABUND_FILE   <- "data/Taula_abundàncies_KOs.tsv" # definim la ruta de la taula d'abundàncies
KO_EC_FILE   <- "Dades_Pere/enzymes_KO.json" # definim la ruta que relaciona K00 i enzims
EC_RXN_FILE  <- "Dades_Pere/enzymes_lite.json" # definim la ruta que relaciona reaccions i enzims
RXN_FILE     <- "Dades_Pere/reactions_lite_pw.json" # definim la ruta que ens dona reaccions reversibles
OUT_DIR      <- "graphml_output"
dir.create(OUT_DIR, showWarnings = FALSE)

## ── 1. LOAD DATA ──────────────────────────────────────────────
message("Loading data...")

abund_raw <- read.delim(ABUND_FILE, row.names = 1, check.names = FALSE)
# files = KOs, Columnes = pacients

ko_ec_raw  <- fromJSON(KO_EC_FILE,  simplifyVector = FALSE) # de K0 a enzims
ec_rxn_raw <- fromJSON(EC_RXN_FILE, simplifyVector = FALSE) # d'enzims a reaccions
rxn_raw    <- fromJSON(RXN_FILE,    simplifyVector = FALSE) # informació de reaccions

## ── 2. BUILD KO → REACTIONS MAPPING ──────────────────────────
message("Building KO → reactions mapping...")


# EC → reaction IDs
ec_to_rxns <- lapply(
  setNames(ec_rxn_raw, sapply(ec_rxn_raw, function(x) x[["enzyme"]])), 
  function(ec_entry) {
    rxns <- unlist(ec_entry[["reaction"]])
  }
)


# KO → unique reaction IDs  (via EC intermediary)
ko_to_rxns <- lapply(ko_ec_raw, function(ec_list) {
  rxns <- unlist(lapply(unlist(ec_list), function(ec) ec_to_rxns[[ec]]))
})
names(ko_to_rxns) <- names(ko_ec_raw)

# Drop KOs that map to no reactions
ko_to_rxns <- Filter(function(x) length(x) > 0, ko_to_rxns)

## ── 3. BUILD REACTION METADATA TABLE ─────────────────────────
message("Building reaction metadata...")

rxn_meta <- tibble(
  rxn_id     = sapply(rxn_raw, `[[`, "reaction"),           # adjust key as needed
  reversible = sapply(rxn_raw, function(r) isTRUE(r[["reversible"]])),
  substrates = lapply(rxn_raw, function(r) unlist(r[["substrate"]])),
  products   = lapply(rxn_raw, function(r) unlist(r[["product"]]))
)

## ── 4. HELPER: BUILD ONE PATIENT GRAPH ───────────────────────

build_patient_graph <- function(patient_id, abund_vec, ko_to_rxns, rxn_meta) {
  
  # -- 4a. Proportional abundance per reaction -------------------
  active_kos <- names(abund_vec)[abund_vec > 0]
  active_kos <- intersect(active_kos, names(ko_to_rxns))
  
  if (length(active_kos) == 0) {
    message("  No active KOs found for patient: ", patient_id)
    return(NULL)
  }
  
  # Build a long data frame: each row = (rxn, fractional_abundance)
  rxn_abund_long <- lapply(active_kos, function(ko) {
    rxns <- ko_to_rxns[[ko]]
    n    <- length(rxns)
    tibble(rxn_id = rxns, frac_abund = abund_vec[[ko]] / n)
  }) |>
    bind_rows() |>
    summarise(abundance = sum(frac_abund), .by = rxn_id)
  
  # Keep only reactions with abundance > 0 and known metadata
  active_rxns <- rxn_abund_long |>
    filter(abundance > 0) |>
    inner_join(rxn_meta, by = "rxn_id")
  
  if (nrow(active_rxns) == 0) {
    message("  No active reactions for patient: ", patient_id)
    return(NULL)
  }
  
  message("  Patient ", patient_id, ": ", nrow(active_rxns), " active reactions")
  
  # -- 4b. Vectorised edge finding --------------------------------
  #
  # For each reaction we pre-compute:
  #   produced_compounds  = products  ∪ (substrates if reversible)
  #   consumed_compounds  = substrates ∪ (products  if reversible)
  #
  # Then we use a compound → reactions index to find R1→R2 edges
  # without a nested loop.
  
  n_rxn <- nrow(active_rxns)
  
  produced <- mapply(function(prod, subs, rev) {
    if (rev) union(prod, subs) else prod
  }, active_rxns$products, active_rxns$substrates, active_rxns$reversible,
  SIMPLIFY = FALSE)
  
  consumed <- mapply(function(subs, prod, rev) {
    if (rev) union(subs, prod) else subs
  }, active_rxns$substrates, active_rxns$products, active_rxns$reversible,
  SIMPLIFY = FALSE)
  
  # Index: compound → indices of reactions that CONSUME it
  all_consumed_cpds <- unique(unlist(consumed))
  cpd_to_consumer_idx <- setNames(
    lapply(all_consumed_cpds, function(cpd) {
      which(vapply(consumed, function(x) cpd %in% x, logical(1)))
    }),
    all_consumed_cpds
  )
  
  # For every reaction R1, find all R2s that consume something R1 produces
  edge_list <- lapply(seq_len(n_rxn), function(i) {
    cpds_produced <- produced[[i]]
    if (length(cpds_produced) == 0) return(NULL)
    
    # Union of all consumer indices for all produced compounds, excluding self
    targets <- unlist(cpd_to_consumer_idx[cpds_produced], use.names = FALSE)
    targets <- targets[targets != i]
    
    if (length(targets) == 0) return(NULL)
    tibble(from = active_rxns$rxn_id[i],
           to   = active_rxns$rxn_id[targets])
  }) |> bind_rows()
  
  # -- 4c. Build igraph object -----------------------------------
  nodes_df <- active_rxns |>
    select(rxn_id, abundance, reversible) |>
    rename(name = rxn_id)
  
  if (nrow(edge_list) == 0) {
    g <- graph_from_data_frame(
      d        = data.frame(from = character(0), to = character(0)),
      directed = TRUE,
      vertices = nodes_df
    )
  } else {
    g <- graph_from_data_frame(
      d        = edge_list,
      directed = TRUE,
      vertices = nodes_df
    )
  }
  
  graph_attr(g, "patient_id") <- patient_id
  return(g)
}

## ── 5. LOOP OVER ALL PATIENTS ─────────────────────────────────
message("\nProcessing patients...")

patient_ids <- colnames(abund_raw)
graphs      <- list()

for (pid in patient_ids) {
  message("Processing: ", pid)
  abund_vec <- setNames(abund_raw[[pid]], rownames(abund_raw))
  
  g <- build_patient_graph(
    patient_id = pid,
    abund_vec  = abund_vec,
    ko_to_rxns = ko_to_rxns,
    rxn_meta   = rxn_meta
  )
  
  if (!is.null(g)) {
    graphs[[pid]] <- g
    out_path <- file.path(OUT_DIR, paste0(pid, "_metabolic_graph.graphml"))
    write_graph(g, out_path, format = "graphml")
    message("  Saved → ", out_path)
  }
}

message("\nDone. ", length(graphs), " graphs exported to '", OUT_DIR, "/'")

## ── 6. OPTIONAL: QUICK SUMMARY STATS ─────────────────────────
summary_df <- tibble(
  patient_id   = names(graphs),
  n_reactions  = sapply(graphs, vcount),
  n_edges      = sapply(graphs, ecount),
  density      = sapply(graphs, graph.density),
  n_components = sapply(graphs, function(g) components(g)$no)
)

print(summary_df)
write.csv(summary_df, file.path(OUT_DIR, "graph_summary.csv"), row.names = FALSE)
```