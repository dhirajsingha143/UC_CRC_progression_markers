# Required libraries
library(limma)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Automated Contrast generation based on disease column
  generate_pairwise_contrasts <- function(conditions) {
      levels <- unique(as.character(conditions))
      pairs <- combn(levels, 2, simplify = FALSE)
      contrast_strings <- sapply(pairs, function(pair) {
        paste0(pair[1], "_vs_", pair[2], " = ", pair[1], " - ", pair[2])
      })
      return(contrast_strings)
    }

# DEG analysis function
  run_deg_analysis <- function(expr_data, meta_data, condition_col = "disease", logfc_threshold = 2, pval_threshold = 0.05, contrasts_list = NULL) {

  
  message("\n🔍 Running DEG Analysis ...\n_🛠SET Working directory to save results/DEG")

  # Step 1: Log transformation
  #log_expr_data <- log2(expr_data + 1)

  # Step 2: Quantile Normalization
  #norm_expr_data <- normalizeBetweenArrays(log_expr_data, method = "quantile")

  # Step 3: Design Matrix
    meta_data[[condition_col]] <- factor(meta_data[[condition_col]])

    design <- model.matrix(~ 0 + meta_data[[condition_col]])
    colnames(design) <- levels(meta_data[[condition_col]])


  # Step 4: Contrast Matrix
  if (is.null(contrasts_list)) {
    stop("❌ You must provide a named list of contrasts.")
  }

  contrast_matrix <- makeContrasts(contrasts = contrasts_list, levels = design)

  # Step 5: Linear Model and Fitting
  fit <- lmFit(expr_data, design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit3 <- eBayes(fit2)

  # Step 6: Collect DEG Results per Contrast
  deg_results <- list()

  for (contrast_name in colnames(contrast_matrix)) {

    cat(paste0("\n🧬 Processing contrast: ", contrast_name))

    res <- topTable(fit3, coef = contrast_name, number = Inf) %>%
      rownames_to_column("Gene") %>%
      mutate(sig = case_when(
        logFC > logfc_threshold & adj.P.Val < pval_threshold ~ "Up",
        logFC < -logfc_threshold & adj.P.Val < pval_threshold ~ "Down",
        TRUE ~ "Not Sig"
      ))

    up <- res %>% filter(sig == "Up")
    down <- res %>% filter(sig == "Down")

    deg_results[[contrast_name]] <- list(
      full_result = res,
      up = up,
      down = down,
      top100_up = up %>% arrange(desc(logFC)) %>% slice_head(n = 100),
      top100_down = down %>% arrange(logFC) %>% slice_head(n = 100)
    )

    # Optional: Save CSV (commented)
    write.csv(res, paste0(contrast_name, "_DEG_results.csv"), row.names = FALSE)
    message("\n 📍 Currently ACTIVATED saving results.CSV for testing purpose")

    # Optional: Volcano Plot
    p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
      geom_point(aes(color = sig), alpha = 0.7) +
      scale_color_manual(values = c("Down" = "#1b7837", "Up" = "red", "Not Sig" = "grey")) +
      geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
      geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed") +
      xlim(min(res$logFC), max(res$logFC)) +
      labs(title = paste("Volcano Plot:", contrast_name),
           x = "Log2 Fold Change",
           y = "-log10 Adjusted P-value") +
      theme_minimal()
    
    print(p)
  }

  return(deg_results)
}
