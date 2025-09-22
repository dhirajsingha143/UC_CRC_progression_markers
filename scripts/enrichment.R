run_enrichment <- function(
    gene_list,
    fc_values = NULL,
    prefix,
    outdir = "results",
    save_tables = TRUE,
    showCategory = 15
) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  enrichments <- list(
    BP   = enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE),
    MF   = enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE),
    CC   = enrichGO(gene = gene_list, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE),
    KEGG = {
      k <- enrichKEGG(gene = gene_list, organism = "hsa", pvalueCutoff = 0.05)
      setReadable(k, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    },
    Reactome = enrichPathway(gene = gene_list, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
  )
  
  results <- list()
  
  for (nm in names(enrichments)) {
    enrich_res <- enrichments[[nm]]
    if (is.null(enrich_res) || nrow(as.data.frame(enrich_res)) == 0) {
      message("No enrichment results for ", nm)
      next
    }
    
    # save table if requested
    if (save_tables) {
      write.csv(as.data.frame(enrich_res),
                file.path(outdir, paste0(prefix, "_", nm, ".csv")),
                row.names = FALSE)
    }
    
    # ---- save all plots into a single PDF ----
    pdf(file.path(outdir, paste0(prefix, "_", nm, ".pdf")), width = 12, height = 10)
    
    # standard plots
    print(barplot(enrich_res, showCategory = showCategory) + ggtitle(paste(nm, "Barplot")))
    print(dotplot(enrich_res, showCategory = showCategory) + ggtitle(paste(nm, "Dotplot")))
    
    if (nm %in% c("BP", "MF", "CC")) {
      print(goplot(enrich_res, showCategory = showCategory) + 
              ggtitle(paste(nm, "Goplot")) + 
              theme(plot.title = element_text(hjust = 0.5)))
    }
    
    # cnetplot with fold changes if provided
    if (!is.null(fc_values)) {
      print(
        cnetplot(enrich_res,
                 showCategory = min(showCategory, 10),
                 foldChange = fc_values,
                 layout = "circle",
                 node_label = "all") +
          ggtitle(paste(nm, "Cnetplot")) + 
          theme(plot.title = element_text(hjust = 0.5))
      )
    }
    
    # term similarity (emap & upset)
    sim <- tryCatch(pairwise_termsim(enrich_res), error = function(e) NULL)
    if (!is.null(sim)) {
      print(emapplot(sim, showCategory = showCategory, layout = "fr", color = "p.adjust") +
              ggtitle(paste("Enriched", nm, "emapplot")) + theme(plot.title = element_text(hjust = 0.5)))
    }
    
    print(upsetplot(enrich_res, n = 10, order.by = "freq") +
            ggtitle(paste("UpSet Plot of", nm)))
    
    dev.off()
    
    results[[nm]] <- enrich_res
  }
  
  return(results)
}

