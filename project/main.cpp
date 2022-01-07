#include <sys/time.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <vector>

#include "minmulticut.cpp"

// for time measuring (in milliseconds)
long GetTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  long mt = (long)t.tv_sec * 1000 + t.tv_usec / 1000;
  return mt;
}

////////////////////////////////////////////////////////////////////////////////
// Test with cycle and terminals attached to it
std::vector<Edge> GenerateGraphCycleTerminals(size_t num_vertices) {
  std::vector<Edge> graph;
  for (size_t i = num_vertices; i < 2 * num_vertices - 1; i++) {
    graph.push_back(Edge(i, i + 1, 1));
    graph.push_back(Edge(i, i - num_vertices, 1.9999));
  }
  graph.push_back(Edge(2 * num_vertices - 1, num_vertices, 1));
  graph.push_back(Edge(2 * num_vertices - 1, num_vertices - 1, 1.9999));
  return graph;
}

void TestCycleTerminals(size_t num_vert, std::ofstream& fout) {
  std::vector<Edge> edges = GenerateGraphCycleTerminals(num_vert);
  std::vector<std::vector<double>> weight(2 * num_vert);
  std::vector<size_t> vertices;
  for (size_t i = 0; i < 2 * num_vert; i++) {
    weight[i].resize(2 * num_vert, -1);
  }
  for (size_t i = 0; i < num_vert; i++) {
    vertices.push_back(i);
  }
  for (Edge e : edges) {
    weight[e.v][e.u] = e.weight;
    weight[e.u][e.v] = e.weight;
  }
  AlmostMinMulticutSolver solver(weight, vertices);
  long time = GetTime();
  std::pair<std::vector<Edge>, double> ans = solver.GetSolution();
  time = GetTime() - time;
  fout << num_vert << "," << ans.second << "," << (ans.second / num_vert) << ","
       << time << "\n";
}

void TestCycleTerminalsStatistics() {
  std::ofstream fout("CycleTerminals.csv");
  fout << "N,Answer,Accuracy,Time(ms)\n";
  for (size_t num_vert = 3; num_vert < 100; num_vert++) {
    TestCycleTerminals(num_vert, fout);
  }
  for (size_t num_vert = 100; num_vert < 800; num_vert += 100) {
    TestCycleTerminals(num_vert, fout);
  }
}

/////////////////////////////////////////////////////////////////////////////
// Test with cycle
std::vector<Edge> GenerateGraphCycle(size_t num_vertices) {
  std::random_device random_device;
  std::mt19937 generator(random_device());
  std::uniform_int_distribution<> distribution(0, 1e5);
  std::vector<Edge> graph;
  for (size_t i = 0; i < num_vertices - 1; i++) {
    graph.push_back(Edge(i, i + 1, (double)distribution(generator) / 1e3));
  }
  graph.push_back(
      Edge(0, num_vertices - 1, (double)distribution(generator) / 1e3));
  return graph;
}

void TestCycle(size_t num_vert, size_t num_term, std::ofstream& fout) {
  std::random_device random_device;
  std::mt19937 generator(random_device());
  std::uniform_int_distribution<> distribution(0, num_term - 1);
  std::vector<Edge> edges = GenerateGraphCycle(num_vert);
  std::vector<std::vector<double>> weight(num_vert);
  std::vector<size_t> vertices;
  for (size_t i = 0; i < num_vert; i++) {
    weight[i].resize(num_vert, -1);
  }
  std::set<size_t> vset;
  while (vset.size() < num_term) {
    vset.insert((int)distribution(generator));
  }
  for (auto v : vset) {
    vertices.push_back(v);
  }
  for (Edge e : edges) {
    weight[e.v][e.u] = e.weight;
    weight[e.u][e.v] = e.weight;
  }
  AlmostMinMulticutSolver solver(weight, vertices);
  long time = GetTime();
  std::pair<std::vector<Edge>, double> ans = solver.GetSolution();
  time = GetTime() - time;

  double correct_ans = 0;
  std::sort(vertices.begin(), vertices.end());
  for (int i = 0; i < num_term - 1; i++) {
    int start = vertices[i];
    int end = vertices[(i + 1) % num_term];
    double cur_min = 1e10;
    for (int delta = 0; delta < end - start; delta++) {
      if (weight[start + delta][start + delta + 1] < cur_min) {
        cur_min = weight[start + delta][start + delta + 1];
      }
    }
    correct_ans += cur_min;
  }
  int start = vertices[num_term - 1];
  int end = vertices[0];
  double cur_min = 1e10;
  for (int delta = 0; delta < end + num_term - start; delta++) {
    if (weight[(start + delta) % num_vert][(start + delta + 1) % num_vert] <
        cur_min) {
      cur_min =
          weight[(start + delta) % num_vert][(start + delta + 1) % num_vert];
    }
  }

  fout << num_vert << "," << num_term << "," << ans.second << ","
       << (correct_ans / ans.second) << "," << time << "\n";
}

void TestCycleStatistics() {
  std::ofstream fout("Cycle.csv");
  fout << "Vertices,Terminals,Answer,Accuracy,Time(ms)\n";
  for (size_t num_vert = 3; num_vert < 10; num_vert++) {
    for (size_t num_term = 3; num_term < num_vert + 1; num_term++) {
      TestCycle(num_vert, num_term, fout);
    }
  }
  for (size_t num_vert = 100; num_vert < 700; num_vert += 100) {
    for (size_t num_term = 10; num_term < num_vert + 1; num_term += 60) {
      TestCycle(num_vert, num_term, fout);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Test with small random graphs
std::vector<Edge> GenerateGraph(size_t num_vertices) {
  std::random_device random_device;
  std::mt19937 generator(random_device());
  std::uniform_int_distribution<> vertice(0, num_vertices - 1);
  std::uniform_int_distribution<> weight(-1e2, 1e5);
  std::vector<Edge> graph;
  for (size_t i = 0; i < (size_t)num_vertices * sqrt(num_vertices); i++) {
    int v1 = vertice(generator);
    int v2 = vertice(generator);
    if (v1 == v2) {
      continue;
    }
    double w = (double)weight(generator) / 1e3;
    if (w < 0) {
      w = -1;
    }
    graph.push_back(Edge(v1, v2, w));
  }
  return graph;
}

bool ReacheableBFS(std::vector<std::vector<double>>& weight,
                   std::vector<size_t>& vertices, size_t s) {
  std::vector<int> distance;
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
  bool reachable = false;
  for (auto v : vertices) {
    if (v == s) {
      continue;
    }
    reachable |= (distance[v] != -1);
  }

  return reachable;
}

double NaiveAlgo(size_t num_vert, std::vector<size_t>& vertices,
                 std::vector<Edge>& edges, std::vector<size_t>& ind,
                 double cur_min) {
  if (ind.size() == edges.size()) {
    std::vector<std::vector<double>> weight(num_vert);
    for (size_t i = 0; i < num_vert; i++) {
      weight[i].resize(num_vert, -1);
    }
    for (size_t i = 0; i < edges.size(); i++) {
      if (ind[i]) {
        Edge e = edges[i];
        weight[e.v][e.u] = e.weight;
        weight[e.u][e.v] = e.weight;
      }
    }

    bool reachable = false;
    for (auto v : vertices) {
      reachable |= ReacheableBFS(weight, vertices, v);
    }
    if (!reachable) {
      double new_min = 0;
      for (size_t i = 0; i < edges.size(); i++) {
        if (!ind[i]) {
          new_min += edges[i].weight;
        }
      }
      cur_min = std::min(cur_min, new_min);
    }
  } else {
    ind.push_back(0);
    cur_min = NaiveAlgo(num_vert, vertices, edges, ind, cur_min);
    ind.pop_back();
    ind.push_back(1);
    cur_min = NaiveAlgo(num_vert, vertices, edges, ind, cur_min);
    ind.pop_back();
  }
  return cur_min;
}

void TestSmallGraph(size_t num_vert, size_t num_term, std::ofstream& fout) {
  std::random_device random_device;
  std::mt19937 generator(random_device());
  std::uniform_int_distribution<> distribution(0, num_term - 1);
  std::vector<Edge> edges = GenerateGraph(num_vert);
  std::vector<std::vector<double>> weight(num_vert);
  std::vector<size_t> vertices;
  for (size_t i = 0; i < num_vert; i++) {
    weight[i].resize(num_vert, -1);
  }
  std::set<size_t> vset;
  while (vset.size() < num_term) {
    vset.insert((int)distribution(generator));
  }
  for (auto v : vset) {
    vertices.push_back(v);
  }
  for (Edge e : edges) {
    weight[e.v][e.u] = e.weight;
    weight[e.u][e.v] = e.weight;
  }
  AlmostMinMulticutSolver solver(weight, vertices);
  long time = GetTime();
  std::pair<std::vector<Edge>, double> ans = solver.GetSolution();
  time = GetTime() - time;

  edges.clear();
  for (size_t i = 0; i < num_vert; i++) {
    for (size_t j = i + 1; j < num_vert; j++) {
      if (weight[i][j] > -1e-3) {
        edges.push_back(Edge(i, j, weight[i][j]));
      }
    }
  }
  std::vector<size_t> ind;
  long naive_time = GetTime();
  double correct_ans = NaiveAlgo(num_vert, vertices, edges, ind, 1e10);
  naive_time = GetTime() - naive_time;

  fout << num_vert << "," << num_term << "," << ans.second << ","
       << (correct_ans / ans.second) << "," << time << "," << naive_time
       << "\n";
}

void TestSmallGraphStatistics() {
  std::ofstream fout("SmallGraphs.csv");
  fout << "Vertices,Terminals,Answer,Accuracy,Time(ms),NaiveAlgo Time(ms)\n";

  for (size_t num_vert = 3; num_vert < 9; num_vert++) {
    for (size_t num_term = 3; num_term < num_vert + 1; num_term++) {
      for (size_t k = 0; k < 5; k++) {
        TestSmallGraph(num_vert, num_term, fout);
      }
    }
  }
}
int main() {
  // TestCycleTerminalsStatistics();
  // TestCycleStatistics();
  // TestSmallGraphStatistics();
}