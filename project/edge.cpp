#include <algorithm>

struct Edge {
  size_t v;
  size_t u;
  double weight;

  Edge(size_t vertice1, size_t vertice2) {
    v = std::min(vertice1, vertice2);
    u = std::max(vertice1, vertice2);
  }
  Edge(size_t vertice1, size_t vertice2, double w) {
    v = std::min(vertice1, vertice2);
    u = std::max(vertice1, vertice2);
    weight = w;
  }
};

bool operator<(Edge e1, Edge e2) {
  return e1.v < e2.v || e1.v == e2.v && e1.u < e2.u;
}