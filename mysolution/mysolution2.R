# mysolution/mysolution2.R

library(igraph)
set.seed(42)

# 1. Graf Barabási–Albert: 1000 węzłów, m = 1
g_ba <- barabasi.game(
  n = 1000,
  m = 1,
  directed = FALSE
)

# 2. Wizualizacja (duży graf → małe węzły)
par(mar = c(1, 1, 1, 1))
cat("\n--- Wizualizacja grafu ---\n")

plot(
  g_ba,
  layout = layout_with_fr,
  vertex.size = 1.5,
  vertex.label = NA,
  edge.arrow.size = 0.1,
  main = "Graf Barabási–Albert (1000 węzłów)"
)

# 3. Betweenness centrality – najbardziej centralny węzeł
cat("\n--- Obliczanie betweenness ---\n")
bw <- betweenness(g_ba)

max_bw_node <- which.max(bw)

cat("Najbardziej centralny węzeł:", max_bw_node, "\n")
cat("Wartość betweenness:", max(bw), "\n")

# 4. Średnica grafu
cat("\n--- Średnica grafu ---\n")
diam <- diameter(g_ba)
cat("Średnica:", diam, "\n")

# 5. Różnica między modelami BA i ER:
#
# Model Erdős–Rényi (ER):
# - Graf powstaje przez losowe łączenie par węzłów z ustalonym prawdopodobieństwem p.
# - Stopnie węzłów mają rozkład zbliżony do Poissona.
# - W grafie nie pojawiają się wyraźne “huby” – większość węzłów ma podobny stopień.
#
# Model Barabási–Albert (BA):
# - Graf rośnie w czasie — nowe węzły dołączają do już istniejących.
# - Preferencyjne przyłączanie: węzły o dużym stopniu są częściej wybierane.
# - Stopnie mają rozkład potęgowy: wiele węzłów o małym stopniu, kilka ogromnych hubów.
# - Model lepiej odzwierciedla rzeczywiste sieci, np. społeczne czy internetowe.
