// msDP(G, k, Δ): decides whether the connected graph G has a tree decomposition
// of width <= k whose separators are taken from Δ (Tamaki, SEA 2019, slides 20-29).
//
// Fixed vertex order 0 < 1 < ... < n-1; vertex sets are compared by their minimum.
// A connected set C is *inbound* if some full component D of G \ N(C) has
// min D < min C.  The main loop processes every inbound C with N(C) ∈ Δ in
// increasing order of |C| and marks it feasible when ISFEASIBLE(C) holds:
//
//   ISFEASIBLE(C):
//     if |N[C]| <= k+1 return true
//     for all marked D with min(C) ∈ D, D ⊊ C, |N(C) ∪ N(D)| <= k+1:
//        if ALLFEASIBLE(N(C) ∪ N(D), C) return true
//     return ALLFEASIBLE(N(C) ∪ {min C}, C)
//   ALLFEASIBLE(S, C):
//     if |S| > k+1 return false
//     for each component D of G \ S inside C:
//        if N(D) = S:  require ISFEASIBLE(D)          (recursion)
//        else:         require D marked feasible
//     return true
//
// G is feasible iff ISFEASIBLE(V(G)) (N(V) = ∅).  Every successful call records
// the bag S (or the leaf bag N[C]) that witnessed it, from which a tree
// decomposition is rebuilt.
//
// Deviations from the slide pseudocode (all preserve soundness, i.e. a YES
// answer always comes with a valid decomposition of width <= k):
//  * line 2 does not require N(C) ∈ Δ (a single bag N[C] is always valid);
//  * in ALLFEASIBLE a non-full component D with |N[D]| <= k+1 is accepted even
//    if it is not marked (again a single bag suffices);
//  * negative results are memoized until the next set gets marked, positive
//    results are memoized permanently.
#pragma once
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graph.hpp"

class MsDP {
 public:
  MsDP(const Graph& G, int k, const std::vector<Bits>& delta) : G_(G), k_(k), byMin_(G.n) {
    for (const Bits& S : delta) {
      std::vector<Bits> full;
      for (Bits& c : G_.components(G_.all() - S))
        if (G_.nbr(c) == S) full.push_back(std::move(c));
      // components() returns components ordered by min vertex: all full
      // components but the first are inbound.
      for (size_t i = 1; i < full.size(); ++i) L_.push_back(std::move(full[i]));
    }
    std::stable_sort(L_.begin(), L_.end(),
                     [](const Bits& a, const Bits& b) { return a.count() < b.count(); });
  }

  // Runs the DP; returns true iff a decomposition of width <= k was found.
  bool run() {
    for (const Bits& C : L_) {
      if (isFeasible(C) && marked_.insert(C).second) {
        byMin_[C.first()].add(C, G_.nbr(C));
        ++gen_;
      }
    }
    if (G_.n == 0) return feasible_ = true;
    return feasible_ = isFeasible(G_.all());
  }

  TreeDecomp decomposition() const {
    TreeDecomp td;
    if (!feasible_ || G_.n == 0) return td;
    buildRec(G_.all(), td);
    return compressTD(td, G_.n);
  }

  size_t numInbound() const { return L_.size(); }
  size_t numMarked() const { return marked_.size(); }

 private:
  struct Witness {
    bool leaf;
    Bits S;
  };

  const Graph& G_;
  int k_;
  std::vector<Bits> L_;
  // Marked sets D grouped by min D, stored contiguously for a fast linear scan.
  struct Bucket {
    std::vector<uint64_t> words;  // per entry: W words of D, then W words of N(D)
    std::vector<int> size;        // |D|
    void add(const Bits& D, const Bits& ND) {
      words.insert(words.end(), D.w.begin(), D.w.end());
      words.insert(words.end(), ND.w.begin(), ND.w.end());
      size.push_back(D.count());
    }
  };
  std::vector<Bucket> byMin_;
  std::unordered_set<Bits, BitsHash> marked_;
  std::unordered_map<Bits, Witness, BitsHash> pos_;  // proven feasible
  std::unordered_map<Bits, int, BitsHash> neg_;      // proven infeasible at generation
  int gen_ = 0;
  bool feasible_ = false;

  bool isFeasible(const Bits& C) {
    if (pos_.count(C)) return true;
    if (auto it = neg_.find(C); it != neg_.end() && it->second == gen_) return false;
    Bits NC = G_.nbr(C);
    int nc = NC.count(), cc = C.count();
    if (nc <= k_ && cc + nc <= k_ + 1) {
      pos_.emplace(C, Witness{true, Bits()});
      return true;
    }
    if (nc <= k_) {
      int m = C.first();
      const size_t W = C.w.size();
      const Bucket& bk = byMin_[m];
      std::unordered_set<Bits, BitsHash> tried;
      for (size_t i = 0; i < bk.size.size(); ++i) {
        if (bk.size[i] >= cc) continue;
        const uint64_t* d = bk.words.data() + 2 * W * i;
        const uint64_t* nd = d + W;
        bool ok = true;
        int u = 0;
        for (size_t j = 0; j < W && ok; ++j) {
          if (d[j] & ~C.w[j]) ok = false;  // D ⊆ C
          u += std::popcount(nd[j] | NC.w[j]);
        }
        if (!ok || u > k_ + 1) continue;
        Bits S = NC;
        for (size_t j = 0; j < W; ++j) S.w[j] |= nd[j];
        if (!tried.insert(S).second) continue;
        if (allFeasible(S, C)) {
          pos_.emplace(C, Witness{false, std::move(S)});
          return true;
        }
      }
      Bits S = NC;
      S.set(m);
      if (allFeasible(S, C)) {
        pos_.emplace(C, Witness{false, std::move(S)});
        return true;
      }
    }
    neg_[C] = gen_;
    return false;
  }

  bool allFeasible(const Bits& S, const Bits& C) {
    if (S.count() > k_ + 1) return false;
    std::vector<Bits> fulls;
    for (Bits& D : G_.components(C - S)) {
      Bits ND = G_.nbr(D);
      if (ND == S) {
        fulls.push_back(std::move(D));
      } else if (marked_.count(D)) {
        // ok
      } else if (D.count() + ND.count() <= k_ + 1) {
        pos_.try_emplace(D, Witness{true, Bits()});
      } else {
        return false;
      }
    }
    for (const Bits& D : fulls)
      if (!isFeasible(D)) return false;
    return true;
  }

  // Builds a tree decomposition of G[N[C]] whose returned root bag contains N(C).
  int buildRec(const Bits& C, TreeDecomp& td) const {
    const Witness& w = pos_.at(C);
    int id = static_cast<int>(td.bags.size());
    if (w.leaf) {
      td.bags.push_back((C | G_.nbr(C)).toVector());
      return id;
    }
    td.bags.push_back(w.S.toVector());
    for (const Bits& D : G_.components(C - w.S)) {
      int child = buildRec(D, td);
      td.edges.push_back({id, child});
    }
    return id;
  }
};
