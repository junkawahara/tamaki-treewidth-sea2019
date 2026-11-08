// Greedy upper bound (min-fill elimination) and lower bound (MMD+ / min-d).
#pragma once
#include <algorithm>
#include <vector>

#include "graph.hpp"

// Tree decomposition from a min-fill elimination ordering.
inline TreeDecomp minFillTD(const Graph& G) {
  int n = G.n;
  if (n == 0) return TreeDecomp{};
  std::vector<Bits> H = G.adj;
  Bits alive = G.all();
  std::vector<int> order, pos(n);
  std::vector<Bits> bagOf(n);
  for (int step = 0; step < n; ++step) {
    int best = -1;
    long long bestFill = 0;
    int bestDeg = 0;
    alive.forEach([&](int v) {
      Bits Nv = H[v] & alive;
      long long fill = 0;
      Nv.forEach([&](int u) { fill += (Nv - H[u]).count() - 1; });
      int deg = Nv.count();
      if (best < 0 || fill < bestFill || (fill == bestFill && deg < bestDeg)) {
        best = v;
        bestFill = fill;
        bestDeg = deg;
      }
    });
    int v = best;
    Bits Nv = H[v] & alive;
    Nv.forEach([&](int u) { H[u] |= Nv; H[u].reset(u); });
    bagOf[v] = Nv;
    bagOf[v].set(v);
    pos[v] = step;
    order.push_back(v);
    alive.reset(v);
  }
  TreeDecomp td;
  for (int i = 0; i < n; ++i) td.bags.push_back(bagOf[order[i]].toVector());
  for (int i = 0; i < n; ++i) {
    int v = order[i];
    int parent = -1;
    bagOf[v].forEach([&](int u) {
      if (u != v && (parent < 0 || pos[u] < parent)) parent = pos[u];
    });
    if (parent < 0 && i + 1 < n) parent = i + 1;  // different component: link anywhere
    if (parent >= 0) td.edges.push_back({i, parent});
  }
  return compressTD(td, n);
}

// MMD+ (min-d heuristic): contract a minimum degree vertex into its neighbor of minimum degree.
inline int mmdPlusLowerBound(const Graph& G) {
  int n = G.n;
  std::vector<Bits> H = G.adj;
  Bits alive = G.all();
  int lb = 0, left = n;
  while (left >= 2) {
    int v = -1, dv = 0;
    alive.forEach([&](int x) {
      int d = (H[x] & alive).count();
      if (v < 0 || d < dv) {
        v = x;
        dv = d;
      }
    });
    lb = std::max(lb, dv);
    if (dv > 0) {
      int u = -1, du = 0;
      (H[v] & alive).forEach([&](int x) {
        int d = (H[x] & alive).count();
        if (u < 0 || d < du) {
          u = x;
          du = d;
        }
      });
      (H[v] & alive).forEach([&](int x) {
        if (x == u) return;
        H[x].set(u);
        H[u].set(x);
      });
    }
    alive.reset(v);
    --left;
  }
  return lb;
}
