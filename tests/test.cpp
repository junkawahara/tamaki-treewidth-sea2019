// Randomized tests:
//  1. listExact(G, k) == { S ∈ allMinimalSeparators(G) : |S| <= k }
//  2. ascend / descend / alternate return valid decompositions whose width equals
//     the treewidth computed by a brute-force O*(2^n) DP (small n), and agree
//     with each other on larger graphs.
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>

#include "../src/solver.hpp"

static Graph randomGraph(std::mt19937& rng, int n, double p) {
  Graph G(n);
  std::uniform_real_distribution<double> U(0, 1);
  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      if (U(rng) < p) G.addEdge(i, j);
  return G;
}

// Random relabeling (the algorithms depend on the vertex order).
static Graph shuffled(const Graph& G, std::mt19937& rng) {
  std::vector<int> perm(G.n);
  for (int i = 0; i < G.n; ++i) perm[i] = i;
  std::shuffle(perm.begin(), perm.end(), rng);
  Graph H(G.n);
  for (int u = 0; u < G.n; ++u) G.adj[u].forEach([&](int v) { H.addEdge(perm[u], perm[v]); });
  return H;
}

// TW(S) = min_{v ∈ S} max(TW(S \ v), |Q(S \ v, v)|)   (Bodlaender et al.)
static int bruteTreewidth(const Graph& G) {
  int n = G.n;
  if (n == 0) return -1;
  std::vector<uint32_t> adj(n, 0);
  for (int v = 0; v < n; ++v) G.adj[v].forEach([&](int u) { adj[v] |= 1u << u; });
  std::vector<int8_t> tw(1u << n, 127);
  tw[0] = -1;
  for (uint32_t S = 1; S < (1u << n); ++S) {
    int best = 127;
    for (int v = 0; v < n; ++v) {
      if (!(S >> v & 1)) continue;
      uint32_t R = S & ~(1u << v);
      // Q(R, v): vertices outside R ∪ {v} reachable from v through R
      uint32_t seen = 1u << v, frontier = 1u << v, q = 0;
      while (frontier) {
        int x = __builtin_ctz(frontier);
        frontier &= frontier - 1;
        uint32_t nb = adj[x] & ~seen;
        seen |= nb;
        q |= nb & ~R;
        frontier |= nb & R;
      }
      int val = std::max<int>(tw[R], __builtin_popcount(q));
      best = std::min(best, val);
    }
    tw[S] = static_cast<int8_t>(best);
  }
  return std::max(0, static_cast<int>(tw[(1u << n) - 1]));
}

// "<mode>" uses the reductions, "<mode>-p0" runs the plain slide algorithm.
static const char* const kModes[] = {"ascend", "descend", "alternate", "ascend-p0", "descend-p0", "alternate-p0"};
static void setMode(Options& o, const std::string& m) {
  o.reduce = m.size() < 3 || m.substr(m.size() - 3) != "-p0";
  o.mode = o.reduce ? m : m.substr(0, m.size() - 3);
}

static int failures = 0;
#define CHECK(cond, msg)                                 \
  do {                                                   \
    if (!(cond)) {                                       \
      ++failures;                                        \
      std::cerr << "FAIL: " << msg << std::endl;         \
    }                                                    \
  } while (0)

static void dump(const Graph& G) {
  std::cerr << "p tw " << G.n << ' ' << G.numEdges() << '\n';
  for (int u = 0; u < G.n; ++u) G.adj[u].forEach([&](int v) {
      if (u < v) std::cerr << u + 1 << ' ' << v + 1 << '\n';
    });
}

static bool sameSets(std::vector<Bits> a, std::vector<Bits> b) {
  auto key = [](const Bits& x) { return x.toVector(); };
  std::set<std::vector<int>> A, B;
  for (auto& x : a) A.insert(key(x));
  for (auto& x : b) B.insert(key(x));
  return A == B && A.size() == a.size() && B.size() == b.size();
}

static void testListing(std::mt19937& rng, int trials, int maxN) {
  long long compared = 0, seps = 0;
  for (int t = 0; t < trials; ++t) {
    int n = 3 + rng() % (maxN - 2);
    double p = 0.15 + (rng() % 60) / 100.0;
    Graph G = randomGraph(rng, n, p);
    std::vector<Bits> all = allMinimalSeparators(G);
    // listExact assumes nothing about connectivity in its recurrence, but the
    // solver only calls it on connected graphs; test both anyway.
    if (G.components(G.all()).size() != 1) continue;
    for (int k = 0; k <= n; ++k) {
      std::vector<Bits> want;
      for (auto& S : all)
        if (S.count() <= k) want.push_back(S);
      std::vector<Bits> got = listExact(G, k);
      ++compared;
      seps += want.size();
      if (!sameSets(want, got)) {
        CHECK(false, "listExact mismatch n=" << n << " k=" << k << " want=" << want.size()
                                             << " got=" << got.size());
        dump(G);
        return;
      }
    }
  }
  std::cerr << "  compared " << compared << " (G,k) pairs, " << seps << " separators" << std::endl;
}

static void testTreewidth(std::mt19937& rng, int trials, int maxN) {
  for (int t = 0; t < trials; ++t) {
    int n = 1 + rng() % maxN;
    double p = 0.1 + (rng() % 70) / 100.0;
    Graph G = randomGraph(rng, n, p);
    int want = bruteTreewidth(G);
    for (const char* mode : kModes) {
      Options o;
      setMode(o, mode);
      TreeDecomp td = Solver(o).solve(G);
      std::string err = validateTD(G, td);
      CHECK(err.empty(), mode << " invalid TD: " << err);
      if (td.width() != want) {
        CHECK(false, mode << " width " << td.width() << " != tw " << want << " (n=" << n << ")");
        dump(G);
        return;
      }
    }
  }
}

static Graph grid(int r, int c) {
  Graph G(r * c);
  for (int i = 0; i < r; ++i)
    for (int j = 0; j < c; ++j) {
      if (i + 1 < r) G.addEdge(i * c + j, (i + 1) * c + j);
      if (j + 1 < c) G.addEdge(i * c + j, i * c + j + 1);
    }
  return G;
}

static void testKnown(std::mt19937& rng) {
  struct Case {
    std::string name;
    Graph G;
    int tw;
  };
  std::vector<Case> cases;
  {
    Graph P(10);
    for (int i = 0; i < 5; ++i) {
      P.addEdge(i, (i + 1) % 5);
      P.addEdge(i, i + 5);
      P.addEdge(5 + i, 5 + (i + 2) % 5);
    }
    cases.push_back({"petersen", P, 4});
  }
  cases.push_back({"grid5x5", grid(5, 5), 5});
  cases.push_back({"grid4x7", grid(4, 7), 4});
  {
    Graph C(12);
    for (int i = 0; i < 12; ++i) C.addEdge(i, (i + 1) % 12);
    cases.push_back({"cycle12", C, 2});
  }
  {
    Graph K(7);
    for (int i = 0; i < 7; ++i)
      for (int j = i + 1; j < 7; ++j) K.addEdge(i, j);
    cases.push_back({"K7", K, 6});
  }
  cases.push_back({"empty5", Graph(5), 0});
  for (auto& c : cases) {
    for (int s = 0; s < 3; ++s) {
      Graph G = s == 0 ? c.G : shuffled(c.G, rng);
      for (const char* mode : kModes) {
        Options o;
        setMode(o, mode);
        TreeDecomp td = Solver(o).solve(G);
        std::string err = validateTD(G, td);
        CHECK(err.empty(), c.name << " " << mode << " invalid: " << err);
        CHECK(td.width() == c.tw, c.name << " " << mode << " width " << td.width() << " != " << c.tw);
      }
    }
  }
}

static void testAgreement(std::mt19937& rng, int trials) {
  for (int t = 0; t < trials; ++t) {
    int n = 20 + rng() % 16;
    double p = 0.08 + (rng() % 20) / 100.0;
    Graph G = randomGraph(rng, n, p);
    int w = -2;
    for (const char* mode : kModes) {
      Options o;
      setMode(o, mode);
      TreeDecomp td = Solver(o).solve(G);
      std::string err = validateTD(G, td);
      CHECK(err.empty(), mode << " invalid TD: " << err);
      if (w == -2) w = td.width();
      CHECK(td.width() == w, "modes disagree on n=" << n << ": " << mode << " " << td.width() << " vs " << w);
    }
  }
}

int main(int argc, char** argv) {
  unsigned seed = argc > 1 ? std::atoi(argv[1]) : 12345;
  std::mt19937 rng(seed);
  std::cerr << "listing..." << std::endl;
  testListing(rng, 400, 15);
  testListing(rng, 100, 24);
  std::cerr << "known graphs..." << std::endl;
  testKnown(rng);
  std::cerr << "treewidth vs brute force..." << std::endl;
  testTreewidth(rng, 400, 14);
  std::cerr << "agreement on larger graphs..." << std::endl;
  testAgreement(rng, 30);
  if (failures) {
    std::cerr << failures << " failure(s)" << std::endl;
    return 1;
  }
  std::cerr << "all tests passed" << std::endl;
  return 0;
}
