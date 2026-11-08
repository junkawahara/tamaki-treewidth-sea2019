// Listing minimal separators.
//
//  * closeSeparators / lambdaSuccessors / allMinimalSeparators:
//      the Λ(G) structure of Berry et al. (slides 30-32).
//  * listExact(G, k): Δ_k(G), all minimal separators of size <= k
//      (slides 33-44), via the recurrence on Δ_ab(A, F) with pruning and
//      divide and conquer on X = N(v) for a min-fill vertex v.
//  * HeuristicLister: Δ0 ⊂ Δ1 ⊂ ... ⊆ Δ_k(G) from a tree decomposition (slide 45).
#pragma once
#include <algorithm>
#include <functional>
#include <unordered_set>
#include <vector>

#include "graph.hpp"

using SepSet = std::unordered_set<Bits, BitsHash>;

// Minimal separators close to the connected set B: N(D) for every component D of G \ N[B].
template <class F>
void closeSeparators(const Graph& G, const Bits& B, F&& emit) {
  for (const Bits& D : G.components(G.all() - G.closedNbr(B))) emit(G.nbr(D));
}

// Successors of S in Λ(G): separators close to C ∪ {v} for a full component C of S and v ∈ S.
template <class F>
void lambdaSuccessors(const Graph& G, const Bits& S, F&& emit) {
  for (const Bits& C : G.components(G.all() - S)) {
    if (G.nbr(C) != S) continue;
    S.forEach([&](int v) {
      Bits B = C;
      B.set(v);
      closeSeparators(G, B, emit);
    });
  }
}

// All minimal separators (Berry, Bordat, Cogis 2000).  Exponential output; for tests.
inline std::vector<Bits> allMinimalSeparators(const Graph& G) {
  SepSet seen;
  std::vector<Bits> list;
  auto add = [&](const Bits& S) {
    if (S.any() && seen.insert(S).second) list.push_back(S);
  };
  for (int v = 0; v < G.n; ++v) {
    Bits B = G.empty();
    B.set(v);
    closeSeparators(G, B, add);
  }
  for (size_t i = 0; i < list.size(); ++i) {
    Bits S = list[i];
    lambdaSuccessors(G, S, add);
  }
  return list;
}

// ---------------------------------------------------------------------------
// Exact listing of Δ_k(G).
//
// For a connected A with a ∈ A, b ∉ N[A] such that N(A) is a minimal a-b
// separator with full components A and C_b(A) (the component of G \ N(A)
// containing b), and F ⊆ N(A), Δ_ab(A, F) is the set of minimal a-b
// separators S with
//    |S| <= k,  F ⊆ S,  A ⊆ C_a,  C_b ∩ N[A] = ∅,  |C_a| <= |C_b|
// (C_a, C_b: the components of G \ S containing a, b).  For crossing
// separators we additionally require a = min(C_a ∩ X), see below.
//
//   Δ_ab(A, F) = Δ_ab(A, F ∪ {v}) ∪ Δ_ab(A', F),   v ∈ N(A) \ F,
// where N(A') is the minimal (A ∪ {v})-b separator close to A ∪ {v}.
//
// Base cases / pruning (all are necessary conditions for non-emptiness):
//    |F| > k;  |A| > |C_b(A)|;  |F| = k and N(A) ≠ F;  N(A) = F → {F};
//    F ∩ A' ≠ ∅ (for the second branch);
//    vertices of N(A) \ F adjacent to b must be in S (forced);
//    disjoint paths: take vertex-disjoint paths starting at the vertices of
//      N(A) \ F through V \ N[A] \ {b}.  At most k - |F| of them can be hit
//      by S \ F; every other path lies entirely in C_a.  So
//      |C_a| >= |A| + (sum of the shortest (#paths - (k-|F|)) path sizes) =: lb,
//      and we need lb <= |C_b(A)| and 2·lb <= n - |F|.
//
// Divide and conquer: v = min-fill vertex, X = N(v).  A minimal separator S
// either crosses X (two vertices a, b ∈ X lie in different full components
// of G \ S; then a, b are non-adjacent), or has a full component disjoint
// from X, hence S = X = N(C) or S is local to a component C of G \ N[v].
// The local ones are minimal separators of H_C = G[N[C]] with N(C) completed
// into a clique (N(C) is a minimal separator since {v} and C are full
// components), which are listed recursively.  Every listed set is finally
// re-checked to be a minimal separator of G.
// ---------------------------------------------------------------------------

class ExactLister {
 public:
  explicit ExactLister(int k) : k_(k) {}

  std::vector<Bits> list(const Graph& G) {
    SepSet out;
    listRec(G, out);
    std::vector<Bits> res;
    for (const Bits& S : out)
      if (S.count() <= k_ && S.any() && G.isMinimalSeparator(S)) res.push_back(S);
    return res;
  }

  long long nodes() const { return nodes_; }

 private:
  int k_;
  long long nodes_ = 0;

  struct Ctx {
    const Graph* G;
    int a, b;
    Bits all;
    Bits lowX;  // vertices of X smaller than a
    SepSet* out;
  };

  static int minFillVertex(const Graph& G) {
    int best = -1;
    long long bestFill = -1;
    int bestDeg = 0;
    for (int v = 0; v < G.n; ++v) {
      long long fill = 0;
      const Bits& Nv = G.adj[v];
      Nv.forEach([&](int u) { fill += (Nv - G.adj[u]).count() - 1; });  // -1: u itself
      fill /= 2;
      int deg = Nv.count();
      if (best < 0 || fill < bestFill || (fill == bestFill && deg < bestDeg)) {
        best = v;
        bestFill = fill;
        bestDeg = deg;
      }
    }
    return best;
  }

  void listRec(const Graph& G, SepSet& out) {
    if (G.n <= 2) return;
    int v = minFillVertex(G);
    Bits X = G.adj[v];
    Bits all = G.all();

    // separators crossing X
    std::vector<int> xs = X.toVector();
    for (int a : xs) {
      Bits lowX = G.empty();
      for (int x : xs)
        if (x < a) lowX.set(x);
      for (int b : xs) {
        if (b == a || G.adjacent(a, b)) continue;
        Ctx c{&G, a, b, all, lowX, &out};
        Bits Na = G.adj[a];
        Na.set(a);
        Bits Cb = G.compOf(b, all - Na);
        Bits S0 = G.nbr(Cb);
        Bits A0 = G.compOf(a, all - S0);
        recAB(c, A0, G.empty());
      }
    }

    // separators local to components of G \ N[v]
    Bits rest = all - X;
    rest.reset(v);
    for (const Bits& C : G.components(rest)) {
      Bits NC = G.nbr(C);
      if (NC.count() <= k_) out.insert(NC);
      std::vector<int> vs = (C | NC).toVector();
      Graph H = G.induced(vs);
      std::vector<int> pos(G.n, -1);
      for (size_t i = 0; i < vs.size(); ++i) pos[vs[i]] = static_cast<int>(i);
      NC.forEach([&](int x) {
        NC.forEach([&](int y) {
          if (x != y) H.adj[pos[x]].set(pos[y]);
        });
      });
      SepSet sub;
      listRec(H, sub);
      for (const Bits& T : sub) {
        Bits S = G.empty();
        T.forEach([&](int i) { S.set(vs[i]); });
        out.insert(std::move(S));
      }
    }
  }

  // Lower bound on the number of vertices that must join the a-side.
  int pathBound(const Graph& G, const Bits& all, const Bits& A, const Bits& NA, const Bits& F, int b,
                int budget) {
    std::vector<int> ends = (NA - F).toVector();
    int m = static_cast<int>(ends.size());
    if (m <= budget) return 0;
    Bits freeV = all - A - NA;
    freeV.reset(b);
    std::vector<int> len(m, 1);
    bool progress = true;
    while (progress) {
      progress = false;
      for (int i = 0; i < m; ++i) {
        int u = G.adj[ends[i]].firstCommon(freeV);
        if (u >= 0) {
          freeV.reset(u);
          ends[i] = u;
          ++len[i];
          progress = true;
        }
      }
    }
    std::sort(len.begin(), len.end());
    int s = 0;
    for (int i = 0; i < m - budget; ++i) s += len[i];
    return s;
  }

  void recAB(Ctx& c, const Bits& A, const Bits& F) {
    ++nodes_;
    const Graph& G = *c.G;
    int fc = F.count();
    if (fc > k_) return;
    if (A.intersects(c.lowX)) return;
    Bits NA = G.nbr(A);
    Bits Cb = G.compOf(c.b, c.all - A - NA);
    int ac = A.count(), cbc = Cb.count();
    if (ac > cbc) return;
    if (NA == F) {
      c.out->insert(F);
      return;
    }
    if (fc == k_) return;
    int lb = ac + pathBound(G, c.all, A, NA, F, c.b, k_ - fc);
    if (lb > cbc || 2 * lb > G.n - fc) return;

    Bits cand = NA - F;
    Bits forced = cand & G.adj[c.b];
    if (forced.any()) {
      recAB(c, A, F | forced);
      return;
    }
    // The slides allow any v ∈ N(A) \ F.  We take one with the most
    // neighbors in C_b(A): on a sample of 26 PACE 2017 instances this was
    // about 30% faster in total than taking the smallest vertex.
    int v = -1, bestScore = 0;
    cand.forEach([&](int x) {
      int sc = (G.adj[x] & Cb).count();
      if (v < 0 || sc > bestScore) {
        v = x;
        bestScore = sc;
      }
    });
    {
      Bits F2 = F;
      F2.set(v);
      recAB(c, A, F2);
    }
    Bits B = A;
    B.set(v);
    Bits Cb2 = G.compOf(c.b, c.all - G.closedNbr(B));
    Bits S2 = G.nbr(Cb2);
    Bits A2 = G.compOf(c.a, c.all - S2);
    if (A2.intersects(F)) return;
    recAB(c, A2, F);
  }
};

inline std::vector<Bits> listExact(const Graph& G, int k) { return ExactLister(k).list(G); }

// ---------------------------------------------------------------------------
// Heuristic listing (slide 45).
//  Δ0: minimal separators of size <= k contained in a bag of T; we take, for
//      every bag X and every component C of G \ X, the separators close to C
//      (they are contained in N(C) ⊆ X).
//  Δ_{i+1} = Δ_i ∪ { R ∈ Δ_k(G) : R successor in Λ(G) of some S ∈ Δ_i }.
// ---------------------------------------------------------------------------
class HeuristicLister {
 public:
  HeuristicLister(const Graph& G, const TreeDecomp& T, int k) : G_(G), k_(k) {
    for (const auto& bag : T.bags) {
      Bits X = G.empty();
      for (int v : bag) X.set(v);
      for (const Bits& C : G.components(G.all() - X)) closeSeparators(G, C, [&](const Bits& R) { add(R); });
    }
  }
  const std::vector<Bits>& current() const { return list_; }

  // Δ_i -> Δ_{i+1}; returns false if nothing new was added.
  bool expand() {
    size_t end = list_.size();
    for (size_t i = frontier_; i < end; ++i) {
      Bits S = list_[i];
      lambdaSuccessors(G_, S, [&](const Bits& R) { add(R); });
    }
    frontier_ = end;
    return list_.size() > end;
  }

 private:
  const Graph& G_;
  int k_;
  SepSet seen_;
  std::vector<Bits> list_;
  size_t frontier_ = 0;

  void add(const Bits& R) {
    if (R.any() && R.count() <= k_ && seen_.insert(R).second) list_.push_back(R);
  }
};
