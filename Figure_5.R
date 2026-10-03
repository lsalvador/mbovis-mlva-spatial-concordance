################################################################################
# Code to generate Figure 5 of the manuscript:
# "Scale- and genotype-dependent spatial concordance of Mycobacterium bovis
#  in cattle and badgers"
#
# Authors: Assel Akhmetova, Liliana C. M. Salvador, Katarina Oravcova,
# Eleanor Presho, Carl McCormick, Suzan Thompson, Lorraine Wright,
# Fraser Menzies, Nigel Trimble, Roland Harwood, Rowland R. Kao,
# Adrian Allen, Theo Pepler, Robin Skuce
#
# Code developers: Assel Akhmetova, Liliana C. M. Salvador.
#
# Description: Code to generat Figure 5. Clonal complex diagram of M. bovis 
# MLVA types. The minimum spanning tree for 36 MLVA types using PHYLOViZ. MLVA006
# is shown as the proposed founder of the clonal group; MLVA001 and MLVA004 are 
# predicted subgroup founders of linked clusters.  Red labels between MLVA types
# (nodes) on the graph demonstrate links between profiles according to tiebreak 
# rules implemented in the goeBurst algorithm implemented in PHYLOViZ.     
#
# Note: working directory setup needs to be set up. This code is divided into three
# stages. Stage 2 needs manual intevention, so the code should be run independently
# for each stage
#
# Input:  input_files_phyloviz/isolate_data.txt
#         input_files_phyloviz/mst_edges_phyloviz.csv (based on phyloviz-goeburst analysis)
#
# Output: Stage 1: 
#            edges.rds
#            nodes.rds
#         Stage 2:
#            positions.json
#         Stage 3:
#            Figure_5.pdf
#
# Requirements: R 4.4.2; packages: igraph, visNetwork, shiny, jsonlite, 
# webshot2, magick)
#
# Contact: lilianasalvador@arizona.edu
# Last updated: October 1, 2026
##########################################################################


install.packages(c("igraph", "visNetwork", "shiny", "jsonlite", "webshot2", "magick"))

library(igraph)
library(visNetwork)
library(shiny)
library(jsonlite)
library(webshot2)
library(magick)


#####################
# Stage 1
#####################
################################################################################
# MAIN SCRIPT — LOAD DATA, BUILD GRAPH, CONSTRUCT NODES/EDGES
#
# Run this first, top to bottom. It builds `nodes` and `edges` and saves
# them to disk so 01_arrange_positions_shiny.R and 02_build_static_export.R
# can pick them up.
################################################################################

# All working/output files live here. Change this one line if you move
# the project folder — must match project_dir in the other two scripts.
# project dir needs to be defined
#project_dir <-

################################################################################
# 1. LOAD DATA
################################################################################
cat("Loading data...\n")

freq_data <- read.csv(file.path(project_dir, "input_files/isolate_data.txt"), sep = "\t", stringsAsFactors = FALSE)
names(freq_data) <- c("ST", "Freq")

mst_edges <- read.csv(file.path(project_dir, "input_filmst_edges_phyloviz.csv"), stringsAsFactors = FALSE)

################################################################################
# 2. BUILD GRAPH
################################################################################
g_mst <- graph_from_data_frame(
  d = mst_edges,
  directed = FALSE,
  vertices = data.frame(
    name = as.character(freq_data$ST),
    Freq = freq_data$Freq,
    stringsAsFactors = FALSE
  )
)

################################################################################
# 3. PREPARE NODE DATA FOR visNetwork
################################################################################
cat("Preparing node data...\n")

freq_in_mst <- V(g_mst)$Freq
min_freq <- min(freq_in_mst)
max_freq <- max(freq_in_mst)

size_min <- 80
size_max <- 200
node_sizes <- size_min +
  (freq_in_mst - min_freq) / (max_freq - min_freq) * (size_max - size_min)

# Color palette - greys, teals, blues, greens
color_palette <- c(
  "#2F4F4F", "#4A4A4A", "#696969", "#808080", "#A9A9A9", "#C0C0C0", "#D3D3D3",
  "#008B8B", "#20B2AA", "#48D1CC", "#5FD3D3", "#7FFFD4", "#AFEEEE",
  "#000080", "#00008B", "#0000CD", "#0000FF", "#4169E1", "#1E90FF", "#6495ED", "#87CEEB", "#ADD8E6",
  "#006400", "#228B22", "#2E8B57", "#3CB371", "#66CDAA", "#90EE90", "#98FB98", "#B0E0E6",
  "#4B0082", "#9370DB", "#DA70D6"
)

st_names <- V(g_mst)$name
n_sts <- length(st_names)
node_colors <- color_palette[((seq_len(n_sts) - 1) %% length(color_palette)) + 1]

nodes <- data.frame(
  id = st_names,
  label = st_names,
  size = node_sizes,
  color = node_colors,
  title = paste0("ST ", st_names, "<br>Frequency: ", V(g_mst)$Freq),
  stringsAsFactors = FALSE
)

################################################################################
# 4. PREPARE EDGE DATA FOR visNetwork
################################################################################
cat("Preparing edge data...\n")

# Custom layout weights for ST-6 neighbors (used as "length", the actual
# vis.js property for per-edge spring length override — NOT "weight",
# which vis.js edges don't recognize as a layout property)
st6_long_neighbors <- c("421", "297", "464", "423", "158", "266")
st6_neighbors <- neighbors(g_mst, "6")

custom_lengths <- E(g_mst)$weight
for (i in seq_along(E(g_mst))) {
  endpoints <- ends(g_mst, i)
  
  if ("6" %in% endpoints && (endpoints[1] %in% st6_long_neighbors | endpoints[2] %in% st6_long_neighbors)) {
    custom_lengths[i] <- 20 / custom_lengths[i]
  } else if ("6" %in% endpoints) {
    custom_lengths[i] <- 15 / custom_lengths[i]
  } else if (endpoints[1] %in% st6_neighbors | endpoints[2] %in% st6_neighbors) {
    custom_lengths[i] <- 10 / custom_lengths[i]
  } else {
    custom_lengths[i] <- 5 / custom_lengths[i]
  }
}

edges_data <- igraph::as_data_frame(g_mst, what = "edges")

edges <- data.frame(
  from = edges_data$from,
  to = edges_data$to,
  label = edges_data$weight,        # actual loci-distance, used as the displayed edge label
  length = custom_lengths,          # layout-only spring length override
  weight = edges_data$weight,       # kept for reference / tooltip use
  title = paste0("Distance: ", edges_data$weight, " loci"),
  color = "grey50",
  width = 2,
  stringsAsFactors = FALSE
)

cat("✓ Built", nrow(nodes), "nodes and", nrow(edges), "edges\n")

################################################################################
# 5. SAVE nodes/edges FOR THE NEXT STEPS
################################################################################
saveRDS(nodes, file.path(project_dir, "nodes.rds"))
saveRDS(edges, file.path(project_dir, "edges.rds"))

cat("\nSaved nodes.rds and edges.rds to:\n  ", project_dir, "\n")
cat("\nNext: run Stage 2 below \n")



##############################
# Stage 2
#############################

################################################################################
# STAGE 2 (Shiny version) — ARRANGE NODES AND SAVE POSITIONS
#
# IMPORTANT: after building `nodes` and `edges` in your main script, save
# them to disk first. This Shiny app loads them from that fixed folder rather than relying on
# the working directory or shared session state — RStudio's "Run App"
# button launches the app in a separate background process that does NOT
# have access to objects in your console, and the working directory it
# starts in may not match your console's
#
# Run this stage 2 script. A browser/viewer window opens with your network. 
# Drag nodes into place, then click "Save Positions" — this writes positions.json 
# to the same folder, and you can close the app.
################################################################################

# All working files (nodes.rds, edges.rds, positions.json) live here.
# Change this one line if you move the project folder.

nodes <- readRDS(file.path(project_dir, "nodes.rds"))
edges <- readRDS(file.path(project_dir, "edges.rds"))

ui <- fluidPage(
  actionButton("save_btn", "Save Positions", class = "btn-warning"),
  br(), br(),
  visNetworkOutput("network", height = "1200px")
)

server <- function(input, output, session) {
  
  output$network <- renderVisNetwork({
    
    # Correct node label font the same way: the original construction used
    # `family = "Arial"` (not a real vis.js option — only `face` controls
    # font family) and `face = "bold"` (which tries to use "bold" itself as
    # a font-family name). Override explicitly with the right property names.
    nodes_display <- nodes
    if (!"font.size" %in% names(nodes_display)) nodes_display$font.size <- 70
    nodes_display$font.face <- "Arial"
    nodes_display$font.bold <- TRUE
    if (!"font.color" %in% names(nodes_display)) nodes_display$font.color <- "black"
    
    # Auto-contrast label color per node: black text disappears against the
    # darker entries in the palette (e.g. "#2F4F4F", "#000080"), so pick
    # white or black per node based on background luminance instead of one
    # fixed color for everyone.
    nodes_display$font.color <- sapply(nodes_display$color, function(col) {
      rgb_vals <- grDevices::col2rgb(col)
      luminance <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
      if (luminance < 140) "white" else "black"
    })
    
    nodes_display$shape <- "circle"
    nodes_display$widthConstraint <- 2 * nodes_display$size
    
    # Ensure edge weight labels exist AND are character type (vis.js's label
    # field is documented as String — a numeric label, even if populated,
    # can silently fail to render).
    edges_display <- edges
    if (!"label" %in% names(edges_display) && "weight" %in% names(edges_display)) {
      edges_display$label <- edges_display$weight
    }
    edges_display$label <- as.character(edges_display$label)
    
    # Per-edge font.* columns (baked in earlier in the pipeline) take
    # precedence over visEdges()'s global font default, so the global
    # size = 64 below gets silently overridden unless we set it here too.
    edges_display$font.size <- 64
    edges_display$font.color <- "red"
    edges_display$font.bold <- TRUE
    edges_display$font.face <- "Arial"
    edges_display$font.strokeWidth <- 6
    edges_display$font.strokeColor <- "white"
    edges_display$font.align <- "horizontal"
    
    visNetwork(nodes_display, edges_display, height = "1200px", width = "100%") %>%
      visPhysics(
        enabled = TRUE,
        stabilization = list(iterations = 200),
        barnesHut = list(
          gravitationalConstant = -80000,
          centralGravity = 0.3,
          springLength = 300,
          springConstant = 0.05
        )
      ) %>%
      visInteraction(
        navigationButtons = TRUE,
        keyboard = TRUE,
        zoomView = TRUE,
        dragView = TRUE
      ) %>%
      visLayout(randomSeed = 42) %>%
      visEdges(smooth = FALSE, font = list(size = 64, color = "red", strokeWidth = 6, strokeColor = "white", align = "horizontal")) %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 2),
        nodesIdSelection = TRUE
      ) %>%
      visEvents(stabilizationIterationsDone = "function () { this.setOptions({physics:false}); }")
  })
  
  observeEvent(input$save_btn, {
    visNetworkProxy("network") %>% visGetPositions()
  })
  
  observeEvent(input$network_positions, {
    positions <- input$network_positions
    
    positions_df <- data.frame(
      id = names(positions),
      x = vapply(positions, function(p) p$x, numeric(1)),
      y = vapply(positions, function(p) p$y, numeric(1)),
      stringsAsFactors = FALSE
    )
    
    write_json(positions_df, file.path(project_dir, "positions.json"), auto_unbox = TRUE)
    showNotification(paste0("Positions saved to ", file.path(project_dir, "positions.json"), " — you can close this app."),
                     type = "message", duration = NULL)
  })
}

shinyApp(ui, server)


######################
# STAGE 3
######################
################################################################################
# STAGE 3 — REBUILD WITH HAND-ARRANGED POSITIONS, EXPORT HIGH-RES FIGURE
# Run this after you've dragged nodes in the Shiny app and saved
# positions.json. Loads `nodes` and `edges` from the same RDS files 
# the Shiny app used, and all output files (HTML, PNG, TIFF, PDF) are written to
# the same project folder.
################################################################################

nodes <- readRDS(file.path(project_dir, "nodes.rds"))
edges <- readRDS(file.path(project_dir, "edges.rds"))

################################################################################
# 1. READ BACK YOUR HAND-ARRANGED POSITIONS
################################################################################
cat("Reading saved node positions...\n")

# The Shiny app writes positions.json as a
# flat array of {id, x, y} objects, so fromJSON gives us a data frame directly.
positions_df <- jsonlite::fromJSON(file.path(project_dir, "positions.json"))
positions_df$id <- as.character(positions_df$id)

# Merge into nodes, preserving original node order/attributes
nodes_fixed <- merge(nodes, positions_df, by = "id", all.x = TRUE, sort = FALSE)
nodes_fixed$fixed <- TRUE   # lock these nodes in place — no physics jitter on export

# Correct node label font: the original construction used `family = "Arial"`
# (not a real vis.js option — only `face` controls font family) and
# `face = "bold"` (which tries to use "bold" itself as a font-family name).
# Override explicitly with the right property names.
if (!"font.size" %in% names(nodes_fixed)) nodes_fixed$font.size <- 70
nodes_fixed$font.face <- "Arial"
nodes_fixed$font.bold <- TRUE
if (!"font.color" %in% names(nodes_fixed)) nodes_fixed$font.color <- "black"

# Auto-contrast label color per node (see explanation in Shiny script).
nodes_fixed$font.color <- sapply(nodes_fixed$color, function(col) {
  rgb_vals <- grDevices::col2rgb(col)
  luminance <- 0.299 * rgb_vals[1] + 0.587 * rgb_vals[2] + 0.114 * rgb_vals[3]
  if (luminance < 140) "white" else "black"
})

# Center labels inside nodes (see explanation in Shiny script).
nodes_fixed$shape <- "circle"
nodes_fixed$widthConstraint <- 2 * nodes_fixed$size

cat("Loaded positions for", nrow(positions_df), "nodes\n")

################################################################################
# 2. STATIC/CLEAN VERSION (uses your hand-arranged layout, physics off)
################################################################################
cat("Creating clean version for figure export...\n")

# Ensure edge weight labels exist and are visibly styled (same defensive
# handling as the Shiny arrangement step).
edges_display <- edges
if (!"label" %in% names(edges_display) && "weight" %in% names(edges_display)) {
  edges_display$label <- edges_display$weight
}
edges_display$label <- as.character(edges_display$label)

# Per-edge font.* columns (baked in earlier in the pipeline) take precedence
# over visEdges()'s global font default — override them directly.
edges_display$font.size <- 64
edges_display$font.color <- "red"
edges_display$font.bold <- TRUE
edges_display$font.face <- "Arial"
edges_display$font.strokeWidth <- 6
edges_display$font.strokeColor <- "white"
edges_display$font.align <- "horizontal"

network_static <- visNetwork(nodes_fixed, edges_display, height = "7200px", width = "9600px") %>%
  visPhysics(enabled = FALSE) %>%   # positions are fixed — no simulation needed
  visEdges(smooth = FALSE, font = list(size = 64, color = "red", strokeWidth = 6, strokeColor = "white", align = "horizontal")) %>%
  visInteraction(
    navigationButtons = FALSE,
    dragView = FALSE,
    zoomView = FALSE,
    keyboard = FALSE
  ) %>%
  visOptions(
    highlightNearest = FALSE,
    nodesIdSelection = FALSE
  ) %>%
  # Explicitly fit the view to the larger canvas — don't rely on vis.js's
  # default auto-fit behavior, since that's normally tied to the
  # stabilization event, which never fires with physics disabled.
  htmlwidgets::onRender("function(el, x) { this.fit(); }")

visSave(network_static, file = file.path(project_dir, "Figure_5.html"))
cat("Saved: Figure_5.html (clean, hand-arranged layout)\n")

################################################################################
# 3. HIGH-RESOLUTION CAPTURE
################################################################################
cat("Capturing high-resolution screenshot...\n")

# If it can't find one (e.g. "Google Chrome was not found"), point it at
# any other Chromium-based browser you have installed instead, such as Edge:
Sys.setenv(CHROMOTE_CHROME = "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge")
# Or check what chromote can find automatically with: chromote::find_chrome()

webshot2::webshot(
  url     = file.path(project_dir, "MST_network_for_export.html"),
  file    = file.path(project_dir, "MST_network_final.png"),
  vwidth  = 9600,    # MUST match the visNetwork height/width above —
  vheight = 7200,    # this is what actually drives resolution, not zoom
  zoom    = 1,       # zoom only upscales the already-rendered canvas bitmap
  # (causes pixelation) — real resolution comes from the
  # native container size set above
  delay   = 4        # bigger canvas still — give it a bit more time to render
)
cat("Captured: Figure_5.png\n")

################################################################################
# 4. TRIM WHITESPACE + CONVERT TO TIFF/PDF
################################################################################
cat("Processing final outputs...\n")

img <- image_read(file.path(project_dir, "Figure_5.png"))
print(image_info(img))           # sanity check on pixel dimensions

img_trimmed <- image_trim(img)   # auto-crop surrounding white space
print(image_info(img_trimmed))   

target_dpi <- 600   

image_write(img_trimmed, file.path(project_dir, "Figure_5.png"),
            format = "png")
image_write(img_trimmed, file.path(project_dir, "Figure_5.tiff"),
            format = "tiff", density = paste0(target_dpi, "x", target_dpi))
image_write(img_trimmed, file.path(project_dir, "Figure_5.pdf"),
            format = "pdf", density = paste0(target_dpi, "x", target_dpi))

cat("\n Saved: Figure_5.png  (trimmed, high-res)\n")
cat("Saved: Figure_5.tiff (", target_dpi, "dpi)\n")
cat("Saved: Figure_5.pdf  (", target_dpi, "dpi)\n")

