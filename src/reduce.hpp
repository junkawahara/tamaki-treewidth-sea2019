// Safe reduction rules (not part of the slides; standard preprocessing,
// Bodlaender & Koster): with a lower bound `low` on tw(G),
//   * simplicial v (N(v) a clique):            tw(G) = max(deg v, tw(G - v))
//   * almost simplicial v with deg v <= low:   tw(G) = max(deg v, tw(G'))
//     (N(v) \ {u} is a clique for some u; G' = G - v with N(v) made a clique)
// Eliminated vertices are re-inserted into the decomposition of the reduced
// graph as bags N[v] attached to a bag containing the clique N(v).
#pragma once
#include <algorithm>
#include <vector>

#include "graph.hpp"

struct Reduction {
  Graph reduced;                    // graph on the remaining vertices
  std::vector<int> keep;            // reduced vertex i == original vertex keep[i]
  int low = 0;                      // lower bound on tw(original)
  std::vector<int> elimOrder;       // eliminated vertices in order
  std::vector<std::vector<int>> elimNbrs;  // N(v) at elimination time (original ids)
};

inline Reduction reduceGraph(const Graph& G, int low) {
  Reduction R;
  R.low = low;
  std::vector<Bits> H = G.adj;
  Bits alive = G.all();
  auto isClique = [&](const Bits& S) {
    bool ok = true;
    S.forEach([&](int x) {
      if (!ok) return;
      Bits rest = S - H[x];
      rest.reset(x);
      if (rest.any()) ok = false;
    });
    return ok;
  };
  auto eliminate = [&](int v, const Bits& N) {
    N.forEach([&](int x) {
      H[x] |= N;
      H[x].reset(x);
      H[x].reset(v);
    });
    alive.reset(v);
    R.elimOrder.push_back(v);
    R.elimNbrs.push_back(N.toVector());
  };
  bool changed = true;
  while (changed) {
    changed = false;
    for (int v = 0; v < G.n; ++v) {
      if (!alive.test(v)) continue;
      Bits N = H[v] & alive;
      int d = N.count();
      if (isClique(N)) {
        R.low = std::max(R.low, d);
        eliminate(v, N);
        changed = true;
        continue;
      }
      if (d > R.low) continue;
      // almost simplicial: some u ∈ N with N \ {u} a clique
      bool done = false;
      N.forEach([&](int u) {
        if (done) return;
        Bits M = N;
        M.reset(u);
        if (isClique(M)) done = true;
      });
      if (done) {
        eliminate(v, N);
        changed = true;
      }
    }
  }
  R.keep = alive.toVector();
  Graph Hg(static_cast<int>(R.keep.size()));
  std::vector<int> pos(G.n, -1);
  for (size_t i = 0; i < R.keep.size(); ++i) pos[R.keep[i]] = static_cast<int>(i);
  for (size_t i = 0; i < R.keep.size(); ++i)
    (H[R.keep[i]] & alive).forEach([&](int u) { Hg.adj[i].set(pos[u]); });
  R.reduced = std::move(Hg);
  return R;
}

// Lifts a decomposition of R.reduced to one of the original graph (n vertices).
inline TreeDecomp liftDecomposition(const Reduction& R, const TreeDecomp& tdReduced, int n) {
  TreeDecomp td = tdReduced;
  for (auto& b : td.bags)
    for (int& x : b) x = R.keep[x];
  std::vector<Bits> bag;
  for (const auto& b : td.bags) {
    Bits B(n);
    for (int x : b) B.set(x);
    bag.push_back(B);
  }
  for (int i = static_cast<int>(R.elimOrder.size()) - 1; i >= 0; --i) {
    int v = R.elimOrder[i];
    Bits N(n);
    for (int x : R.elimNbrs[i]) N.set(x);
    int host = -1;
    for (size_t j = 0; j < bag.size() && host < 0; ++j)
      if (N.subsetOf(bag[j])) host = static_cast<int>(j);
    Bits B = N;
    B.set(v);
    int id = static_cast<int>(td.bags.size());
    td.bags.push_back(B.toVector());
    bag.push_back(B);
    if (host >= 0)
      td.edges.push_back({host, id});
    else if (id > 0)
      td.edges.push_back({0, id});  // N(v) empty here only if v was isolated
  }
  return compressTD(td, n);
}
