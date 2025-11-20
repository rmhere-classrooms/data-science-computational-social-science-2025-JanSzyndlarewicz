# mysolution/mysolution1.R

library(igraph)

set.seed(123)

# 1. Sieć Erdős–Rényi: 100 wierzchołków, prawdopodobieństwo krawędzi 0.05
g <- erdos.renyi.game(n = 100, p = 0.05)

cat("\n--- Podsumowanie grafu (przed dodaniem wag) ---\n")
print(summary(g))

# 2. Wypisanie wierzchołków i krawędzi
cat("\n--- Wierzchołki ---\n")
print(V(g))

cat("\n--- Krawędzie ---\n")
print(E(g))

# 3. Losowe wagi krawędzi z zakresu 0.01–1.00
E(g)$weight <- runif(ecount(g), min = 0.01, max = 1)

cat("\n--- Podsumowanie grafu (po dodaniu wag) ---\n")
print(summary(g))

# 4. Stopnie węzłów + histogram
node_degrees <- degree(g)

cat("\n--- Stopnie węzłów ---\n")
print(node_degrees)

hist(
  node_degrees,
  main = "Histogram stopni węzłów",
  xlab = "Stopień",
  ylab = "Liczba węzłów",
  col = "lightblue",
  border = "black"
)

# 5. Liczba składowych spójnych
comp <- components(g)

cat("\n--- Liczba składowych spójnych ---\n")
print(comp$no)

# 6. Wizualizacja grafu (rozmiar = PageRank)
cat("\n--- Wizualizacja grafu (PageRank) ---\n")

pr <- page.rank(g)$vector

# Dopasowanie wielkości do wartości PageRank
V(g)$size <- pr * 200

plot(
  g,
  vertex.label = NA,
  edge.arrow.size = 0.2,
  layout = layout_with_fr,
  main = "Graf Erdős–Rényi (rozmiar węzła = PageRank)"
)