#--------------------------------------------
  # Install Required Packages (if needed)
#----------------------------------------------
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#----------------------------------------------
 # Enrichment Function Module
#----------------------------------------------

run_enrichment <- function(deg_annotated, organism_db = org.Hs.eg.db, save_results = FALSE, output_dir = "results") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  enrichment_results <- list()

  for (contrast_name in names(deg_annotated)) {
    message(paste0("🔬 Enrichment for contrast: ", contrast_name))

    # Extract Entrez IDs
    entrez_up <- deg_annotated[[contrast_name]]$annotated_up$entrezid
    entrez_down <- deg_annotated[[contrast_name]]$annotated_down$entrezid

    enrich_up <- enrich_down <- list()

    # Function to run all 3 GO types
    run_all_go <- function(genes) {
      list(
        BP = enrichGO(gene = genes, OrgDb = organism_db, ont = "BP", readable = TRUE),
        MF = enrichGO(gene = genes, OrgDb = organism_db, ont = "MF", readable = TRUE),
        CC = enrichGO(gene = genes, OrgDb = organism_db, ont = "CC", readable = TRUE)
      )
    }

    if (!is.null(entrez_up) && length(entrez_up) > 10) {
      enrich_up$GO <- run_all_go(entrez_up)
      enrich_up$KEGG <- enrichKEGG(gene = entrez_up, organism = "hsa", pvalueCutoff = 0.05)
    }

    if (!is.null(entrez_down) && length(entrez_down) > 10) {
      enrich_down$GO <- run_all_go(entrez_down)
      enrich_down$KEGG <- enrichKEGG(gene = entrez_down, organism = "hsa", pvalueCutoff = 0.05)
    }

    enrichment_results[[contrast_name]] <- list(
      up = enrich_up,
      down = enrich_down
    )

    # Optional: Save results
    if (save_results) {
      save_go_results <- function(go_list, prefix) {
        for (ont in names(go_list)) {
          write.csv(
            as.data.frame(go_list[[ont]]),
            file.path(output_dir, paste0(contrast_name, "_", prefix, "_GO_", ont, ".csv")),
            row.names = FALSE
          )
        }
      }
      if (!is.null(enrich_up$GO)) save_go_results(enrich_up$GO, "UP")
      if (!is.null(enrich_down$GO)) save_go_results(enrich_down$GO, "DOWN")
    }
  }

  return(enrichment_results)
}
#----------------------------------------------
 # Enrichment Plotting Function Module
#----------------------------------------------

# Function to generate enrichment plots using entrez IDs

run_hub_enrichment <- function(hub_annotated, organism_db = org.Hs.eg.db,
                               save_results = FALSE, output_dir = "results") {


  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    BiocManager::install("clusterProfiler")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    BiocManager::install("org.Hs.eg.db")
  }


  library(clusterProfiler)
  library(org.Hs.eg.db)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Check input type
  if (is.data.frame(hub_annotated)) {
    hub_list <- list("Hub_Genes" = hub_annotated)
  } else {
    hub_list <- hub_annotated
  }

  enrichment_results <- list()

  for (hub_set_name in names(hub_list)) {
    message(paste0("🔬 Enrichment for hub set: ", hub_set_name))

    # Extract Entrez IDs
    hub_df <- hub_list[[hub_set_name]]
    entrez_col <- if ("entrezgene_id" %in% colnames(hub_df)) "entrezgene_id" else "entrezid"
    entrez_ids <- na.omit(unique(hub_df[[entrez_col]]))

    if (length(entrez_ids) < 5) {
      message("⚠️ Only ", length(entrez_ids), " valid Entrez IDs found - skipping enrichment")
      next
    }

    # Run enrichments with error handling
    enrichment_results[[hub_set_name]] <- list()

    tryCatch({
      # GO Enrichment
      enrichment_results[[hub_set_name]]$GO <- list(
        BP = enrichGO(gene = entrez_ids, OrgDb = organism_db, ont = "BP",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE),
        MF = enrichGO(gene = entrez_ids, OrgDb = organism_db, ont = "MF",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE),
        CC = enrichGO(gene = entrez_ids, OrgDb = organism_db, ont = "CC",
                      pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
      )
    }, error = function(e) {
      message("GO enrichment failed: ", e$message)
    })

    tryCatch({
      # KEGG Enrichment
      enrichment_results[[hub_set_name]]$KEGG <- enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    }, error = function(e) {
      message("KEGG enrichment failed: ", e$message)
    })

    # Save results if requested
    if (save_results) {
      save_enrichment_results(enrichment_results[[hub_set_name]],
                              hub_set_name,
                              output_dir)
    }
  }

  return(enrichment_results)
}

# Helper function to save results
save_enrichment_results <- function(enrichment, hub_set_name, output_dir) {
  hub_dir <- file.path(output_dir, "hub_enrichments")
  if (!dir.exists(hub_dir)) dir.create(hub_dir)

  # Save GO results
  if (!is.null(enrichment$GO)) {
    for (ont in names(enrichment$GO)) {
      if (!is.null(enrichment$GO[[ont]])) {
        write.csv(
          as.data.frame(enrichment$GO[[ont]]),
          file.path(hub_dir, paste0(hub_set_name, "_GO_", ont, ".csv")),
          row.names = FALSE
        )
      }
    }
  }

  # Save KEGG results
  if (!is.null(enrichment$KEGG)) {
    write.csv(
      as.data.frame(enrichment$KEGG),
      file.path(hub_dir, paste0(hub_set_name, "_KEGG.csv")),
      row.names = FALSE
    )
  }

  # Save Reactome results
  if (!is.null(enrichment$Reactome)) {
    write.csv(
      as.data.frame(enrichment$Reactome),
      file.path(hub_dir, paste0(hub_set_name, "_Reactome.csv")),
      row.names = FALSE
    )
  }

  # Generate plots
  pdf(file.path(hub_dir, paste0(hub_set_name, "_enrichment_plots.pdf")))

  # GO plots
  if (!is.null(enrichment$GO)) {
    for (ont in names(enrichment$GO)) {
      if (!is.null(enrichment$GO[[ont]]) && nrow(enrichment$GO[[ont]]) > 0) {
        print(clusterProfiler::dotplot(enrichment$GO[[ont]],
                                       showCategory = 15,
                                       title = paste("GO", ont, "Enrichment")))
      }
    }
  }

  # KEGG plot
  if (!is.null(enrichment$KEGG) && nrow(enrichment$KEGG) > 0) {
    print(clusterProfiler::dotplot(enrichment$KEGG,
                                   showCategory = 15,
                                   title = "KEGG Pathway Enrichment"))
  }

  # Reactome plot
  if (!is.null(enrichment$Reactome) && nrow(enrichment$Reactome) > 0) {
    print(clusterProfiler::dotplot(enrichment$Reactome,
                                   showCategory = 15,
                                   title = "Reactome Pathway Enrichment"))
  }

  dev.off()
}


