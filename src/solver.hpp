// Ascend / Descend / Alternate (slide 4) and the driver over connected components.
#pragma once
#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "graph.hpp"
#include "heuristics.hpp"
#include "minseps.hpp"
#include "msdp.hpp"
#include "reduce.hpp"

struct Options {
  std::string mode = "ascend";  // ascend | descend | alternate
  int verbose = 0;
  bool reduce = true;  // safe reductions (simplicial / almost simplicial vertices)
};

class Solver {
 public:
  explicit Solver(Options opt) : opt_(std::move(opt)) {}

  TreeDecomp solve(const Graph& G) {
    if (G.n == 0) return TreeDecomp{};
    if (!opt_.reduce) return solveComponents(G, 0);
    Reduction R = reduceGraph(G, mmdPlusLowerBound(G));
    log(1, "reduced n=" + std::to_string(G.n) + " -> " + std::to_string(R.reduced.n) +
               ", low=" + std::to_string(R.low));
    TreeDecomp t = solveComponents(R.reduced, R.low);
    return liftDecomposition(R, t, G.n);
  }

  // Solves each connected component; widths <= cur are acceptable.
  TreeDecomp solveComponents(const Graph& G, int cur) {
    if (G.n == 0) return TreeDecomp{};
    std::vector<Bits> comps = G.components(G.all());
    std::sort(comps.begin(), comps.end(), [](const Bits& a, const Bits& b) { return a.count() > b.count(); });
    TreeDecomp result;
    for (const Bits& comp : comps) {
      std::vector<int> vs = comp.toVector();
      Graph H = G.induced(vs);
      TreeDecomp t = solveConnected(H, cur);
      cur = std::max(cur, t.width());
      int off = static_cast<int>(result.bags.size());
      for (auto& b : t.bags) {
        for (int& x : b) x = vs[x];
        result.bags.push_back(b);
      }
      for (auto [a, b] : t.edges) result.edges.push_back({a + off, b + off});
      if (off > 0) result.edges.push_back({0, off});
    }
    return result;
  }

  // Connected G; any decomposition of width <= floor is acceptable.
  TreeDecomp solveConnected(const Graph& G, int floor) {
    if (G.n <= floor + 1) return trivialTD(G.n);
    int lb = std::max(floor, mmdPlusLowerBound(G));
    log(1, "component n=" + std::to_string(G.n) + " m=" + std::to_string(G.numEdges()) +
               " lb=" + std::to_string(lb));
    if (opt_.mode == "descend") return descend(G, lb);
    if (opt_.mode == "alternate") return alternate(G, lb);
    return ascend(G, lb);
  }

  // Ascend: k = lb, lb+1, ...; stop at the first YES of msDP(G, k, listExact(G, k)).
  TreeDecomp ascend(const Graph& G, int lb) {
    for (int k = lb;; ++k) {
      TreeDecomp td;
      if (exactStep(G, k, td)) return td;
    }
  }

  // Descend: improve a greedy decomposition with heuristic separator lists,
  // then certify optimality with the exact list.
  TreeDecomp descend(const Graph& G, int lb) {
    TreeDecomp T = minFillTD(G);
    int w = T.width();
    log(1, "min-fill width " + std::to_string(w));
    while (w > lb) {
      int k = w - 1;
      TreeDecomp td;
      if (heuristicImprove(G, T, k, td)) {
        T = td;
        w = T.width();
        continue;
      }
      if (exactStep(G, k, td)) {
        T = td;
        w = T.width();
        continue;
      }
      break;  // tw(G) = w
    }
    return T;
  }

  // Alternate: interleave ascend steps (exact, at k = lb) and descend steps
  // (one heuristic list expansion, at k = ub-1), giving each about the same time.
  TreeDecomp alternate(const Graph& G, int lb) {
    using clock = std::chrono::steady_clock;
    TreeDecomp T = minFillTD(G);
    int ub = T.width();
    log(1, "min-fill width " + std::to_string(ub));
    double tAsc = 0, tDesc = 0;
    bool descDone = false;
    std::unique_ptr<HeuristicLister> H;
    bool first = true;
    while (lb < ub) {
      auto t0 = clock::now();
      if (descDone || tAsc <= tDesc) {
        TreeDecomp td;
        if (exactStep(G, lb, td)) return td;
        ++lb;
        tAsc += std::chrono::duration<double>(clock::now() - t0).count();
      } else {
        int k = ub - 1;
        if (!H) {
          H = std::make_unique<HeuristicLister>(G, T, k);
          first = true;
        }
        if (!first && !H->expand()) {
          descDone = true;
        } else {
          first = false;
          MsDP dp(G, k, H->current());
          bool yes = dp.run();
          log(2, "  heuristic k=" + std::to_string(k) + " |Δ_i|=" + std::to_string(H->current().size()) +
                     (yes ? " YES" : " NO"));
          if (yes) {
            T = dp.decomposition();
            ub = T.width();
            H.reset();
            log(1, "improved to width " + std::to_string(ub));
          }
        }
        tDesc += std::chrono::duration<double>(clock::now() - t0).count();
      }
    }
    return T;
  }

 private:
  Options opt_;

  void log(int level, const std::string& s) const {
    if (opt_.verbose >= level) std::cerr << "c " << s << std::endl;
  }

  bool exactStep(const Graph& G, int k, TreeDecomp& out) {
    auto t0 = std::chrono::steady_clock::now();
    ExactLister lister(k);
    std::vector<Bits> delta = lister.list(G);
    auto t1 = std::chrono::steady_clock::now();
    MsDP dp(G, k, delta);
    bool yes = dp.run();
    auto t2 = std::chrono::steady_clock::now();
    log(1, "exact k=" + std::to_string(k) + " |Δ_k|=" + std::to_string(delta.size()) +
               " inbound=" + std::to_string(dp.numInbound()) + " feasible=" + std::to_string(dp.numMarked()) +
               " nodes=" + std::to_string(lister.nodes()) +
               " list=" + std::to_string(std::chrono::duration<double>(t1 - t0).count()) +
               "s dp=" + std::to_string(std::chrono::duration<double>(t2 - t1).count()) + "s -> " +
               (yes ? "YES" : "NO"));
    if (yes) out = dp.decomposition();
    return yes;
  }

  // Tries msDP(G, k, Δ_i) for Δ_0 ⊂ Δ_1 ⊂ ... built from T.
  bool heuristicImprove(const Graph& G, const TreeDecomp& T, int k, TreeDecomp& out) {
    HeuristicLister H(G, T, k);
    for (int i = 0;; ++i) {
      MsDP dp(G, k, H.current());
      bool yes = dp.run();
      log(2, "  heuristic k=" + std::to_string(k) + " i=" + std::to_string(i) +
                 " |Δ_i|=" + std::to_string(H.current().size()) + (yes ? " YES" : " NO"));
      if (yes) {
        out = dp.decomposition();
        log(1, "improved to width " + std::to_string(out.width()));
        return true;
      }
      if (!H.expand()) return false;
    }
  }
};
