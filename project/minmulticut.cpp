#include <algorithm>
#include <map>
#include <set>
#include <vector>

#include "dinitz.cpp"

class AlmostMinMulticutSolver {
 private:
  std::vector<std::vector<double>> weight;

  // Almost minimum multicut will be found for vertices in
  std::vector<size_t> vertices;

  /**
   * Builds a new graph where all terminal vertices except one are united.
   * Returns new indexes of chosen terminal vertice (s) and new united vertice
   *(t) (we will find (s-t)-cut for these ones further)
   **/
  std::pair<size_t, size_t> BuildNewGraph(
      std::vector<std::vector<double>>& new_weight, size_t cur_vertice,
      std::map<size_t, size_t>& new_index);

 public:
  AlmostMinMulticutSolver(std::vector<std::vector<double>>& weight,
                          std::vector<size_t>& vertices)
      : weight(weight), vertices(vertices) {}

  std::pair<std::vector<Edge>, double> GetSolution() {
    size_t N = weight.size();
    size_t K = vertices.size();

    std::vector<double> cut_sum(K);
    std::vector<std::vector<Edge>> cuts(K);

    for (size_t i = 0; i < K; i++) {
      size_t cur_vertice = vertices[i];
      // Building new graph
      std::vector<std::vector<double>> new_weight(N - K + 2);
      for (size_t j = 0; j < new_weight.size(); j++) {
        new_weight[j].resize(N - K + 2, -1);
      }
      std::map<size_t, size_t> new_index;
      std::pair<size_t, size_t> s_t =
          BuildNewGraph(new_weight, cur_vertice, new_index);
      // Running maxflow algorithm
      DinitzAlgorithm dinitz(new_weight, s_t.first, s_t.second);
      std::pair<std::vector<Edge>, double> result = dinitz.GetCut();
      cut_sum[i] = result.second;

      // Finding min (s-t)-cut for new graph
      std::map<size_t, std::vector<size_t>> old_index;
      for (auto j : new_index) {
        old_index[j.second].push_back(j.first);
      }
      std::vector<Edge> cur_cut;
      for (auto edge : result.first) {
        for (size_t v : old_index[edge.v]) {
          for (size_t u : old_index[edge.u]) {
            if (weight[u][v] >= 0) {
              cur_cut.push_back(Edge(v, u));
            }
          }
        }
      }
      cuts[i] = cur_cut;
    }

    // Choosing the 'heaviest' cut and getting final result
    size_t heaviest_cut = 0;
    for (size_t i = 1; i < K; i++) {
      if (cut_sum[i] > cut_sum[heaviest_cut]) {
        heaviest_cut = i;
      }
    }
    std::set<Edge> final_cut;
    for (size_t i = 0; i < K; i++) {
      if (i == heaviest_cut) {
        continue;
      }
      for (Edge edge : cuts[i]) {
        final_cut.insert(edge);
      }
    }
    std::pair<std::vector<Edge>, double> answer;
    answer.second = 0;
    for (Edge edge : final_cut) {
      answer.first.push_back(edge);
      answer.second += weight[edge.v][edge.u];
    }
    return answer;
  }
};

std::pair<size_t, size_t> AlmostMinMulticutSolver::BuildNewGraph(
    std::vector<std::vector<double>>& new_weight, size_t cur_vertice,
    std::map<size_t, size_t>& new_index) {
  std::set<size_t> vert(vertices.begin(), vertices.end());
  vert.erase(cur_vertice);
  int ind = 1;
  for (size_t i = 0; i < weight.size(); i++) {
    if (vert.find(i) != vert.end()) {
      new_index[i] = 0;
    } else {
      new_index[i] = ind;
      ++ind;
    }
  }

  for (size_t i = 0; i < weight.size(); i++) {
    for (size_t j = 0; j < weight.size(); j++) {
      if (weight[i][j] == -1) {
        continue;
      }
      if (new_weight[new_index[i]][new_index[j]] < 0) {  // == -1
        new_weight[new_index[i]][new_index[j]] += 1 + weight[i][j];
      } else {
        new_weight[new_index[i]][new_index[j]] += weight[i][j];
      }
    }
  }

  for (size_t i = 0; i < new_weight.size(); i++) {
    new_weight[i][i] = -1;
  }

  return std::make_pair(new_index[cur_vertice], 0);
}
