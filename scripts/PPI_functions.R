# Install if not already
if (!requireNamespace("STRINGdb", quietly = TRUE)) {
  BiocManager::install("STRINGdb")
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

# ================================================
# COMPLETE PPI NETWORK ANALYSIS PIPELINE (igraph only)
# ================================================

# Load required packages
library(STRINGdb)
library(igraph)
library(dplyr)
library(scales)

# -------------------------------------------------
# 1. Build PPI Network from DEGs
# -------------------------------------------------
build_ppi_network <- function(deg_out,
                              condition = "HC_vs_RA = HC - RA",
                              regulation = "up",
                              top_n = 100,
                              string_version = "11.5",
                              species_id = 9606,
                              score_thresh = 400,
                              save = TRUE,
                              out_dir = "PPI_network",
                              return_plot = TRUE) {

  # Create output directory if needed
  if (save && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("Created output directory: ", out_dir)
  }

  message("\nBuilding PPI network for condition: ", condition)
  message("----------------------------------------")

  # Extract DEG data
  deg_df <- deg_out[[condition]]
  if (is.null(deg_df)) stop("Condition not found in DEG results")

  # Select regulated genes
  if ("up" %in% regulation) {
    deg_df <- deg_df[[paste0("top", top_n, "_up")]]
  } else if ("down" %in% regulation) {
    deg_df <- deg_df[[paste0("top", top_n, "_down")]]
  } else {
    stop("regulation must be 'up', 'down', or both")
  }

  if (is.null(deg_df)) stop("No DEG data found for specified regulation")

  message("Using ", nrow(deg_df), " ", regulation, "-regulated genes")

  # Connect to STRINGdb
  message("Connecting to STRINGdb...")
  string_db <- STRINGdb$new(
    version = string_version,
    species = species_id,
    score_threshold = score_thresh,
    input_directory = ""
  )

  # Map genes to STRING IDs
  message("Mapping gene symbols to STRING IDs...")
  mapped_df <- string_db$map(deg_df, "Gene", removeUnmappedRows = TRUE)
  mapped_df <- mapped_df %>% filter(!is.na(STRING_id) & STRING_id != "")
  message(nrow(mapped_df), " genes successfully mapped")

  # Get interactions
  message("Retrieving protein interactions...")
  interactions <- string_db$get_interactions(mapped_df$STRING_id)
  interactions_clean <- interactions %>%
    filter(!is.na(from) & !is.na(to)) %>%
    distinct(from, to, .keep_all = TRUE)
  message(nrow(interactions_clean), " interactions retrieved")

  # Build graph object
  message("Constructing network graph...")
  g <- graph_from_data_frame(interactions_clean, directed = FALSE)

  # Only simplify if there are multiple edges
  if (any(which_multiple(g))) {
    g <- simplify(g,
                  remove.multiple = TRUE,
                  remove.loops = TRUE,
                  edge.attr.comb = list(combined_score = "mean"))
  }

  # Add vertex attributes
  V(g)$STRING_id <- V(g)$name
  V(g)$label <- mapped_df$Gene[match(V(g)$name, mapped_df$STRING_id)]
  V(g)$degree <- degree(g)
  V(g)$size <- scales::rescale(V(g)$degree, to = c(3, 10))

  # Set default visual attributes
  V(g)$color <- "skyblue"
  V(g)$frame.color <- "grey20"
  E(g)$width <- scales::rescale(E(g)$combined_score, to = c(0.5, 3))
  E(g)$color <- "gray60"

  # Calculate layout
  layout <- layout_with_kk(g)

  # Create plot function
  create_network_plot <- function() {
    par(mar = c(0,0,2,0)) # Reduce margins
    plot(
      g,
      layout = layout,
      vertex.label =  V(g)$label, #ifelse(V(g)$degree > quantile(V(g)$degree, 0.75), V(g)$label, NA), #V(g)$label for all names
      vertex.label.cex = 0.7,
      vertex.label.color = "black",
      vertex.label.family = "Helvetica",
      vertex.label.dist = 0,         # offset from center
      vertex.label.degree = -pi/4,   # angle offset
      main = paste("PPI Network:", condition, regulation),
      margin = -0.1
    )
  }

  # Save plot if requested
  if (save) {
    png(file.path(out_dir, paste0("ppi_network_", condition, "_", regulation, ".png")),
        width = 1200, height = 1200, res = 300)
    create_network_plot()
    dev.off()
    message("Saved network plot to ", file.path(out_dir, paste0("ppi_network_", condition, "_", regulation, ".png")))
  }

  # Return results
  result <- list(
    graph = g,
    plot = if (return_plot) create_network_plot else NULL,
    mapped_genes = mapped_df,
    interactions = interactions_clean,
    layout = layout
  )

  message("PPI network construction complete\n")
  return(result)
}

# -------------------------------------------------
# 2. Identify Hub Genes (igraph only)
# -------------------------------------------------
identify_hub_genes <- function(g,
                               condition = "HC_vs_RA = HC - RA",
                               regulation = "up",
                               quantile_cutoff = 0.75,
                               save = TRUE,
                               out_dir = "PPI_network",
                               return_plot = TRUE) {

  message("\nIdentifying hub genes for condition: ", condition)
  message("----------------------------------------")

  # Calculate network metrics
  deg <- degree(g)
  bet <- betweenness(g)
  clos <- closeness(g)

  V(g)$degree <- deg
  V(g)$betweenness <- bet
  V(g)$closeness <- clos

  # Identify hub genes
  hub_thresh <- quantile(deg, quantile_cutoff)
  hub_nodes <- V(g)[deg > hub_thresh]
  message(length(hub_nodes), " hub genes identified (degree > ", round(hub_thresh, 1), ")")

  # Create hub gene data frame
  hub_df <- data.frame(
    STRING_ID = hub_nodes$name,
    Gene_Symbol = hub_nodes$label,
    Degree = hub_nodes$degree,
    Betweenness = hub_nodes$betweenness,
    Closeness = hub_nodes$closeness,
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Degree))

  # Save hub gene information
  if (save) {
    hub_file <- file.path(out_dir, paste0("hub_genes_", condition, "_", regulation, ".csv"))
    write.csv(hub_df, hub_file, row.names = FALSE)
    message("Saved hub gene data to ", hub_file)
  }

  # Create subgraph with hubs and neighbors
  neighbors_list <- ego(g, order = 1, nodes = hub_nodes, mode = "all")
  hub_plus_neighbors <- unique(unlist(neighbors_list))
  g_hub <- induced_subgraph(g, vids = hub_plus_neighbors)

  # Style the subgraph
  V(g_hub)$color <- ifelse(V(g_hub)$name %in% hub_nodes$name, "red", "orange")
  V(g_hub)$size <- scales::rescale(V(g_hub)$degree, to = c(5, 15))
  V(g_hub)$label.color <- "black"
  V(g_hub)$label.cex <- 0.7

  # Create plot function
  create_hub_plot <- function() {
    par(mar = c(0,0,2,0))
    plot(
      g_hub,
      layout = layout_with_kk(g_hub),
      main = paste("Hub Genes and First Neighbors:", condition, regulation)
    )
  }

  # Save plot if requested
  if (save) {
    png(file.path(out_dir, paste0("hub_network_", condition, "_", regulation, ".png")),
        width = 1200, height = 1200, res = 200)
    create_hub_plot()
    dev.off()
    message("Saved hub network plot to ", file.path(out_dir, paste0("hub_network_", condition, "_", regulation, ".png")))
  }

  # Return results
  result <- list(
    hub_df = hub_df,
    hub_graph = g_hub,
    plot = if (return_plot) create_hub_plot else NULL
  )

  message("Hub identification complete\n")
  return(result)
}

# -------------------------------------------------
# 3. Detect PPI Modules (igraph only)
# -------------------------------------------------
detect_ppi_modules <- function(g,
                               condition = "HC_vs_RA = HC - RA",
                               regulation = "up",
                               algorithm = "walktrap",
                               save = TRUE,
                               out_dir = "PPI_network",
                               return_plot = TRUE) {

  message("\nDetecting modules for condition: ", condition)
  message("----------------------------------------")

  # Perform clustering
  if (algorithm == "walktrap") {
    clust <- cluster_walktrap(g)
  } else if (algorithm == "louvain") {
    clust <- cluster_louvain(g)
  } else {
    stop("Unsupported algorithm. Choose 'walktrap' or 'louvain'")
  }

  modules <- communities(clust)
  message(length(modules), " modules detected")

  # Calculate network metrics
  V(g)$degree <- degree(g)
  V(g)$clustering_coeff <- transitivity(g, type = "local", isolates = "zero")

  # Create module data frames and plots
  module_info_list <- list()
  plot_list <- list()

  for (i in seq_along(modules)) {
    nodes <- modules[[i]]
    subg <- induced_subgraph(g, vids = nodes)

    # Create module data frame
    df <- data.frame(
      gene = V(subg)$label,
      degree = degree(subg),
      clustering_coeff = transitivity(subg, type = "local", isolates = "zero"),
      module = paste0("Module_", i),
      stringsAsFactors = FALSE
    ) %>% arrange(desc(degree))

    module_info_list[[i]] <- df

    # Create module plot function with proper environment capture
    plot_module <- local({
      module_num <- i  # Capture current module number
      module_nodes <- nodes  # Capture current nodes
      module_subg <- subg  # Capture current subgraph
      module_color <- rainbow(length(modules))[module_num]  # Capture color

      function() {
        par(mar = c(0,0,2,0))
        plot(
          module_subg,
          layout = layout_with_fr(module_subg),
          vertex.size = 6,
          vertex.label.cex = 0.6,
          vertex.color = module_color,
          main = paste("Module", module_num, ":", condition, regulation)
        )
      }
    })

    plot_list[[paste0("Module_", i)]] <- plot_module

    # Save module-specific files if requested
    if (save) {
      # Save module data
      module_file <- file.path(out_dir, paste0("module_", i, "_", condition, "_", regulation, ".csv"))
      write.csv(df, module_file, row.names = FALSE)

      # Save module plot
      png(file.path(out_dir, paste0("module_", i, "_", condition, "_", regulation, ".png")),
          width = 2000, height = 1800, res = 300)
      plot_module()
      dev.off()
    }
  }

  # Combine all module stats
  all_stats <- do.call(rbind, module_info_list)

  # Save combined stats if requested
  if (save) {
    stats_file <- file.path(out_dir, paste0("all_modules_", condition, "_", regulation, ".csv"))
    write.csv(all_stats, stats_file, row.names = FALSE)
    message("Saved module data to ", stats_file)
  }

  # Create summary plot function
  plot_all_modules <- function() {
    par(mar = c(0,0,2,0))
    plot(
      clust,
      g,
      vertex.size = 5,
      vertex.label.cex = 0.5,
      main = paste("All Modules:", condition, regulation)
    )
  }

  # Save summary plot if requested
  if (save) {
    png(file.path(out_dir, paste0("all_modules_", condition, "_", regulation, ".png")),
        width = 2000, height = 1800, res = 300)
    plot_all_modules()
    dev.off()
    message("Saved module plots to ", out_dir)
  }

  # Return results
  result <- list(
    module_stats = all_stats,
    cluster_obj = clust,
    module_plots = if (return_plot) plot_list else NULL,
    all_modules_plot = if (return_plot) plot_all_modules else NULL
  )

  message("Module detection complete\n")
  return(result)
}

# -------------------------------------------------
# 4. Run Complete PPI Analysis
# -------------------------------------------------
run_ppi_analysis <- function(deg_out_annotated,
                             conditions = "HC_vs_RA = HC - RA",
                             regulations = c("up", "down"),
                             top_n = 100,
                             save = TRUE,
                             out_dir = "PPI_results") {

  # Create output directory if needed
  if (save && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("Created output directory: ", out_dir)
  }

  all_results <- list()

  for (cond in conditions) {
    for (reg in regulations) {
      message("\nStarting PPI analysis for: ", cond, " (", reg, ")")
      message("========================================")

      # Build PPI network
      network <- build_ppi_network(
        deg_out_annotated = deg_out_annotated,
        condition = cond,
        regulation = reg,
        top_n = top_n,
        save = save,
        out_dir = out_dir,
        return_plot = TRUE
      )

      # Identify hub genes
      hubs <- identify_hub_genes(
        g = network$graph,
        condition = cond,
        regulation = reg,
        save = save,
        out_dir = out_dir,
        return_plot = TRUE
      )

      # Detect modules
      modules <- detect_ppi_modules(
        g = network$graph,
        condition = cond,
        regulation = reg,
        save = save,
        out_dir = out_dir,
        return_plot = TRUE
      )

      # Store results
      all_results[[paste0(cond, "_", reg)]] <- list(
        network = network,
        hubs = hubs,
        modules = modules
      )

      message("Completed PPI analysis for: ", cond, " (", reg, ")\n")
    }
  }

  return(all_results)
}





