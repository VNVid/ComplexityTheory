#include <algorithm>
#include <cstdint>
#include <queue>
#include <vector>

#include "edge.cpp"

class DinitzAlgorithm {
 private:
  std::vector<std::vector<double>> weight;
  std::vector<int> distance;
  size_t s;
  size_t t;
  const double MAXFLOW = 9 * 1e18;

  /**
   *Returns true if t-vertice is reachable from s-vertice, otherwise returns
   *false. And counts distances from the start vertice s to all others (puts
   *result in 'distance' vector)
   **/
  bool BFS();

  double PushFlowDFS(size_t v, double flow);

 public:
  DinitzAlgorithm(std::vector<std::vector<double>>& weight, unsigned int s,
                  unsigned int t)
      : weight(weight), s(s), t(t) {}

  std::pair<std::vector<Edge>, double> GetCut() {
    double cut_weight = 0;
    while (BFS()) {
      while (double flow = PushFlowDFS(s, MAXFLOW)) {
        cut_weight += flow;
      }
    }

    std::vector<Edge> cut;
    std::vector<size_t> reachable;
    std::vector<size_t> unreacheable;
    for (size_t i = 0; i < weight.size(); i++) {
      if (distance[i] == -1) {
        unreacheable.push_back(i);
      } else {
        reachable.push_back(i);
      }
    }
    for (size_t v : reachable) {
      for (size_t u : unreacheable) {
        if (weight[v][u] >= 0) {
          cut.push_back(Edge(v, u));
        }
      }
    }
    return std::make_pair(cut, cut_weight);
  }
};

bool DinitzAlgorithm::BFS() {
  std::queue<size_t> q;
  distance = std::vector<int>(weight.size(), -1);

  distance[s] = 0;
  q.push(s);
  while (!q.empty()) {
    unsigned int from = q.front();
    q.pop();
    for (size_t to = 0; to < weight.size(); to++) {
      if (distance[to] == -1 && weight[from][to] > 0) {
        distance[to] = distance[from] + 1;
      }
    }
  }

  return distance[t] != -1;
}

double DinitzAlgorithm::PushFlowDFS(size_t v, double flow) {
  if (v == t) {
    return flow;
  }
  if (flow == 0) {
    return 0;
  }

  for (size_t to = 0; to < weight.size(); to++) {
    if (distance[to] == distance[v] + 1 && weight[v][to] > 0) {
      double new_flow = PushFlowDFS(to, std::min(flow, weight[v][to]));
      weight[v][to] -= new_flow;
      weight[to][v] -= new_flow;
      return new_flow;
    }
  }
  return 0;
}