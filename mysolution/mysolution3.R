library(shiny)
library(bslib)
library(igraph)

# ================================================================
# 1. Wczytywanie oraz przygotowanie grafu
# ================================================================

build_graph <- function() {
  message("Downloading dataset...")
  
  raw <- read.csv2(
    "https://bergplace.org/share/out.radoslaw_email_email",
    skip = 2, sep = " ", header = FALSE
  )[, 1:2]
  colnames(raw) <- c("sender", "receiver")
  
  g <- graph_from_data_frame(raw, directed = TRUE)
  g <- simplify(g, remove.loops = TRUE)
  
  message("Nodes: ", vcount(g), "  Edges: ", ecount(g))
  message("Computing edge probabilities...")
  
  # Liczność połączeń
  edge_freq <- as.data.frame(table(raw$sender, raw$receiver))
  names(edge_freq) <- c("sender", "receiver", "cnt")
  edge_freq <- subset(edge_freq, cnt > 0)
  
  total_sent <- aggregate(cnt ~ sender, data = edge_freq, sum)
  colnames(total_sent) <- c("sender", "total")
  
  edge_freq <- merge(edge_freq, total_sent, by = "sender")
  edge_freq$prob <- edge_freq$cnt / edge_freq$total
  
  # Dopasowanie wag do krawędzi grafu
  edge_list <- as_edgelist(g, names = TRUE)
  keys_graph <- paste(edge_list[,1], edge_list[,2], sep = "_")
  keys_table <- paste(edge_freq$sender, edge_freq$receiver, sep = "_")
  
  prob_map <- setNames(edge_freq$prob, keys_table)
  E(g)$weight <- prob_map[keys_graph]
  E(g)$weight[is.na(E(g)$weight)] <- 0.001
  
  message("Caching adjacency lists...")
  
  n <- vcount(g)
  edg <- as_edgelist(g, names = FALSE)
  wgt <- E(g)$weight
  
  neigh_cache <- vector("list", n)
  weight_cache <- vector("list", n)
  
  for (i in seq_len(n)) {
    idx <- which(edg[,1] == i)
    neigh_cache[[i]] <- edg[idx, 2]
    weight_cache[[i]] <- wgt[idx]
  }
  
  g$neighbors <- neigh_cache
  g$weights <- weight_cache
  
  message("Graph ready.")
  return(g)
}


# ================================================================
# 2. Symulacja modelu niezależnej kaskady
# ================================================================

cascade <- function(g, start_nodes, multiplier = 1, max_iter = 50) {
  n <- vcount(g)
  
  active <- rep(FALSE, n)
  active[start_nodes] <- TRUE
  
  current <- start_nodes
  counts <- c(length(start_nodes))
  
  attempted <- new.env(hash = TRUE)
  
  for (step in seq_len(max_iter)) {
    if (length(current) == 0) break
    
    next_round <- integer(0)
    
    for (node in current) {
      neigh <- g$neighbors[[node]]
      w <- g$weights[[node]]
      
      if (length(neigh) == 0) next
      
      for (k in seq_along(neigh)) {
        target <- neigh[k]
        if (active[target]) next
        
        key <- paste0(node, "_", target)
        if (exists(key, envir = attempted)) next
        assign(key, TRUE, envir = attempted)
        
        p <- min(w[k] * multiplier, 1)
        
        if (runif(1) <= p) {
          active[target] <- TRUE
          next_round <- c(next_round, target)
        }
      }
    }
    
    next_round <- unique(next_round)
    counts <- c(counts, length(next_round))
    current <- next_round
  }
  
  return(counts)
}

# ================================================================
# 3. Wybór węzłów startowych
# ================================================================

choose_start_nodes <- function(g, method, pct = 0.05) {
  k <- max(1, round(vcount(g) * pct))
  
  if (method == "outdegree") {
    return(order(degree(g, mode="out"), decreasing=TRUE)[1:k])
  }
  if (method == "betweenness") {
    return(order(betweenness(g, directed=TRUE), decreasing=TRUE)[1:k])
  }
  if (method == "closeness") {
    cl <- closeness(g, mode = "out")
    cl[is.infinite(cl)] <- 0
    return(order(cl, decreasing=TRUE)[1:k])
  }
  if (method == "random") {
    return(sample(vcount(g), k))
  }
  if (method == "pagerank") {
    pr <- page.rank(g)$vector
    return(order(pr, decreasing=TRUE)[1:k])
  }
}

# ================================================================
# 4. Uruchamianie wielu prób
# ================================================================

experiment <- function(g, method, multiplier, max_iter, runs = 100) {
  mat <- matrix(0, nrow = runs, ncol = max_iter+1)
  
  for (i in seq_len(runs)) {
    init <- choose_start_nodes(g, method)
    res <- cascade(g, init, multiplier, max_iter)
    len <- min(length(res), max_iter+1)
    mat[i, 1:len] <- res[1:len]
  }
  
  colMeans(mat)
}

ui <- navbarPage(
  title = "Diffusion Simulation",
  theme = bs_theme(
    version = 5,
    bootswatch = "minty",     # Całkowita zmiana stylu
    primary = "#2c7fb8",
    base_font = font_google("Inter")
  ),
  
  tabPanel("Simulation",
           
           fluidPage(
             br(),
             
             fluidRow(
               column(
                 width = 4,
                 
                 # ------ PANEL STEROWANIA W KARTACH ------
                 card(
                   full_screen = FALSE,
                   card_header("Parameters"),
                   card_body(
                     sliderInput(
                       "mult", "Activation probability multiplier:",
                       min = 10, max = 200, step = 10, value = 100, post = "%"
                     ),
                     
                     sliderInput(
                       "iters", "Maximum iterations:",
                       min = 1, max = 50, value = 10
                     ),
                     
                     br(),
                     actionButton("run", "Run simulation",
                                  class = "btn btn-success w-100"),
                     br(), br()
                   )
                 ),
                 
                 card(
                   card_header("Strategies description"),
                   card_body(
                     tags$ul(
                       tags$li("Outdegree –> nodes with many outgoing edges"),
                       tags$li("Betweenness –> bridge nodes"),
                       tags$li("Closeness –> central nodes"),
                       tags$li("Random –> baseline"),
                       tags$li("PageRank –> importance measure")
                     ),
                   )
                 )
               ),
               
               # ------ WYKRES W OSOBNEJ KARCIE ------
               column(
                 width = 8,
                 card(
                   full_screen = TRUE,
                   card_header("Diffusion over iterations"),
                   card_body(
                     plotOutput("plot", height = "650px")
                   )
                 )
               )
             )
           )
  )
)

# ================================================================
# 6. Server logic (unchanged)
# ================================================================

server <- function(input, output) {
  
  g <- build_graph()
  message("System initialized.")
  
  sim <- eventReactive(input$run, {
    mult <- input$mult / 100
    iters <- input$iters
    
    strategies <- c("outdegree", "betweenness", "closeness", "random", "pagerank")
    results <- vector("list", length(strategies))
    
    withProgress(message = "Running experiments...", value = 0, {
      t0 <- Sys.time()
      
      for (i in seq_along(strategies)) {
        incProgress(1/length(strategies),
                    detail = paste("Strategy:", strategies[i]))
        results[[i]] <- experiment(g, strategies[i], mult, iters, runs = 100)
      }
      
      dt <- difftime(Sys.time(), t0, units = "secs")
      
      list(results = results, mult = mult, iters = iters, time = dt)
    })
  })
  
  output$plot <- renderPlot({
    s <- sim()
    
    methods <- c("Out-degree", "Betweenness", "Closeness", "Random", "PageRank")
    cols <- c("#2c7fb8", "#d95f02", "#1b9e77", "#7570b3", "#e7298a")
    
    plot(0:s$iters, s$results[[1]][1:(s$iters+1)], type = "l", lwd = 3,
         col = cols[1], ylim = c(0, max(unlist(s$results))),
         xlab = "Iteration", ylab = "New activations",
         main = paste0("Diffusion (×", s$mult, ") – ", round(s$time,1), "s"))
    
    for (i in 2:5) {
      lines(0:s$iters, s$results[[i]][1:(s$iters+1)], lwd = 3, col = cols[i])
    }
    
    legend("topright", legend = methods, col = cols, lwd = 3)
    grid()
  })
}

shinyApp(ui, server)

