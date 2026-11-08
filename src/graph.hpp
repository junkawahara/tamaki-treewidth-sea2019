// Graph, tree decomposition, and PACE 2017 I/O.
#pragma once
#include <algorithm>
#include <charconv>
#include <climits>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

#include "bits.hpp"

struct Graph {
  int n = 0;
  std::vector<Bits> adj;

  Graph() = default;
  explicit Graph(int n_) : n(n_), adj(n_, Bits(n_)) {}

  void addEdge(int u, int v) {
    if (u == v) return;
    adj[u].set(v);
    adj[v].set(u);
  }
  bool adjacent(int u, int v) const { return adj[u].test(v); }
  int degree(int v) const { return adj[v].count(); }

  Bits empty() const { return Bits(n); }
  Bits all() const {
    Bits b(n);
    for (int i = 0; i < n; ++i) b.set(i);
    return b;
  }
  // Open neighborhood N(S).
  Bits nbr(const Bits& S) const {
    Bits r(n);
    S.forEach([&](int v) { r |= adj[v]; });
    r -= S;
    return r;
  }
  // Closed neighborhood N[S].
  Bits closedNbr(const Bits& S) const { return nbr(S) | S; }

  // Connected component containing v inside the vertex set avail (v must be in avail).
  Bits compOf(int v, const Bits& avail) const {
    Bits comp(n);
    comp.set(v);
    Bits frontier = comp, nb(n);
    const size_t W = comp.w.size();
    while (true) {
      std::fill(nb.w.begin(), nb.w.end(), 0);
      frontier.forEach([&](int u) { nb |= adj[u]; });
      bool grew = false;
      for (size_t i = 0; i < W; ++i) {
        uint64_t x = nb.w[i] & avail.w[i] & ~comp.w[i];
        frontier.w[i] = x;
        comp.w[i] |= x;
        grew |= x != 0;
      }
      if (!grew) return comp;
    }
  }
  // Components of G[avail], ordered by their minimum vertex.
  std::vector<Bits> components(Bits avail) const {
    std::vector<Bits> out;
    while (avail.any()) {
      Bits c = compOf(avail.first(), avail);
      avail -= c;
      out.push_back(std::move(c));
    }
    return out;
  }
  // Induced subgraph on vs; vertex vs[i] becomes i.
  Graph induced(const std::vector<int>& vs) const {
    Graph H(static_cast<int>(vs.size()));
    std::vector<int> pos(n, -1);
    for (size_t i = 0; i < vs.size(); ++i) pos[vs[i]] = static_cast<int>(i);
    for (size_t i = 0; i < vs.size(); ++i)
      adj[vs[i]].forEach([&](int u) {
        if (pos[u] >= 0) H.adj[i].set(pos[u]);
      });
    return H;
  }
  long long numEdges() const {
    long long m = 0;
    for (int v = 0; v < n; ++v) m += adj[v].count();
    return m / 2;
  }
  bool isMinimalSeparator(const Bits& S) const {
    int full = 0;
    for (const Bits& c : components(all() - S))
      if (nbr(c) == S && ++full >= 2) return true;
    return false;
  }
};

struct TreeDecomp {
  std::vector<std::vector<int>> bags;
  std::vector<std::pair<int, int>> edges;

  int width() const {
    int w = -1;
    for (const auto& b : bags) w = std::max(w, static_cast<int>(b.size()) - 1);
    return w;
  }
};

// Single bag containing all vertices.
inline TreeDecomp trivialTD(int n) {
  TreeDecomp td;
  std::vector<int> b(n);
  for (int i = 0; i < n; ++i) b[i] = i;
  td.bags.push_back(b);
  return td;
}

// Repeatedly contract tree edges whose one end bag is a subset of the other.
inline TreeDecomp compressTD(const TreeDecomp& in, int n) {
  int nb = static_cast<int>(in.bags.size());
  if (nb <= 1) return in;
  std::vector<Bits> bag(nb, Bits(n));
  for (int i = 0; i < nb; ++i)
    for (int v : in.bags[i]) bag[i].set(v);
  std::vector<std::vector<int>> nbr(nb);
  for (auto [a, b] : in.edges) {
    nbr[a].push_back(b);
    nbr[b].push_back(a);
  }
  std::vector<char> alive(nb, 1);
  bool changed = true;
  while (changed) {
    changed = false;
    for (int i = 0; i < nb; ++i) {
      if (!alive[i]) continue;
      for (int j : nbr[i]) {
        if (!alive[j] || !bag[i].subsetOf(bag[j])) continue;
        // merge i into j
        for (int x : nbr[i]) {
          if (x == j) continue;
          std::replace(nbr[x].begin(), nbr[x].end(), i, j);
          nbr[j].push_back(x);
        }
        nbr[j].erase(std::remove(nbr[j].begin(), nbr[j].end(), i), nbr[j].end());
        nbr[i].clear();
        alive[i] = 0;
        changed = true;
        break;
      }
    }
  }
  std::vector<int> id(nb, -1);
  TreeDecomp out;
  for (int i = 0; i < nb; ++i)
    if (alive[i]) {
      id[i] = static_cast<int>(out.bags.size());
      out.bags.push_back(bag[i].toVector());
    }
  for (int i = 0; i < nb; ++i)
    if (alive[i])
      for (int j : nbr[i])
        if (i < j) out.edges.push_back({id[i], id[j]});
  return out;
}

// Returns an empty string if td is a valid tree decomposition of G, else a message.
// A bag listing a vertex twice is rejected, since width() counts list entries.
inline std::string validateTD(const Graph& G, const TreeDecomp& td) {
  int nb = static_cast<int>(td.bags.size());
  if (G.n > 0 && nb == 0) return "no bags";
  if (nb > 0 && static_cast<int>(td.edges.size()) != nb - 1) return "number of tree edges != bags-1";
  std::vector<Bits> bag(nb, Bits(G.n));
  for (int i = 0; i < nb; ++i)
    for (int v : td.bags[i]) {
      if (v < 0 || v >= G.n) return "vertex out of range";
      if (bag[i].test(v)) return "vertex " + std::to_string(v + 1) + " appears twice in bag " + std::to_string(i + 1);
      bag[i].set(v);
    }
  std::vector<std::vector<int>> adj(nb);
  for (auto [a, b] : td.edges) {
    if (a < 0 || b < 0 || a >= nb || b >= nb || a == b) return "bad tree edge";
    adj[a].push_back(b);
    adj[b].push_back(a);
  }
  // connectivity of the tree
  if (nb > 0) {
    std::vector<char> seen(nb, 0);
    std::vector<int> st{0};
    seen[0] = 1;
    int cnt = 1;
    while (!st.empty()) {
      int x = st.back();
      st.pop_back();
      for (int y : adj[x])
        if (!seen[y]) {
          seen[y] = 1;
          ++cnt;
          st.push_back(y);
        }
    }
    if (cnt != nb) return "tree is not connected";
  }
  for (int u = 0; u < G.n; ++u) {
    // vertex coverage + connectedness of bags containing u
    int start = -1, total = 0;
    for (int i = 0; i < nb; ++i)
      if (bag[i].test(u)) {
        ++total;
        if (start < 0) start = i;
      }
    if (total == 0) return "vertex " + std::to_string(u + 1) + " not covered";
    std::vector<char> seen(nb, 0);
    std::vector<int> st{start};
    seen[start] = 1;
    int cnt = 1;
    while (!st.empty()) {
      int x = st.back();
      st.pop_back();
      for (int y : adj[x])
        if (!seen[y] && bag[y].test(u)) {
          seen[y] = 1;
          ++cnt;
          st.push_back(y);
        }
    }
    if (cnt != total) return "bags containing vertex " + std::to_string(u + 1) + " are not connected";
    bool ok = true;
    G.adj[u].forEach([&](int v) {
      if (!ok || v < u) return;
      bool covered = false;
      for (int i = 0; i < nb && !covered; ++i) covered = bag[i].test(u) && bag[i].test(v);
      if (!covered) ok = false;
    });
    if (!ok) return "an edge at vertex " + std::to_string(u + 1) + " is not covered";
  }
  return "";
}

// ---- PACE 2017 I/O ----
// Lines starting with 'c' are comments; blank lines are skipped.  Every other
// line must match the format exactly, otherwise std::runtime_error is thrown.

namespace pace_io {
inline std::vector<std::string> tokens(const std::string& line) {
  std::istringstream ss(line);
  std::vector<std::string> t;
  std::string s;
  while (ss >> s) t.push_back(s);
  return t;
}
// The whole token must be a decimal integer in [lo, hi].
inline long long parseInt(const std::string& tok, long long lo, long long hi, const std::string& line) {
  long long x = 0;
  const char* end = tok.data() + tok.size();
  auto [p, ec] = std::from_chars(tok.data(), end, x);
  if (ec != std::errc() || p != end || x < lo || x > hi)
    throw std::runtime_error("bad number '" + tok + "' in line: " + line);
  return x;
}
}  // namespace pace_io

// "p tw n m" followed by exactly m edge lines "u v" (1 <= u, v <= n).
inline Graph readGr(std::istream& in) {
  using namespace pace_io;
  std::string line;
  Graph G;
  bool header = false;
  long long m = 0, edges = 0;
  while (std::getline(in, line)) {
    if (!line.empty() && line[0] == 'c') continue;
    std::vector<std::string> t = tokens(line);
    if (t.empty()) continue;
    if (!header) {
      if (t.size() != 4 || t[0] != "p" || t[1] != "tw") throw std::runtime_error("bad p-line: " + line);
      int n = static_cast<int>(parseInt(t[2], 0, INT_MAX, line));
      m = parseInt(t[3], 0, LLONG_MAX, line);
      G = Graph(n);
      header = true;
      continue;
    }
    if (t.size() != 2) throw std::runtime_error("bad edge line: " + line);
    int u = static_cast<int>(parseInt(t[0], 1, G.n, line));
    int v = static_cast<int>(parseInt(t[1], 1, G.n, line));
    if (++edges > m) throw std::runtime_error("more edges than the " + std::to_string(m) + " declared in the p-line");
    G.addEdge(u - 1, v - 1);
  }
  if (!header) throw std::runtime_error("missing p-line");
  if (edges != m)
    throw std::runtime_error("the p-line declares " + std::to_string(m) + " edges but " + std::to_string(edges) +
                             " were read");
  return G;
}

inline TreeDecomp readTd(std::istream& in) {
  std::string line;
  TreeDecomp td;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == 'c') continue;
    std::istringstream ss(line);
    std::string t;
    ss >> t;
    if (t == "s") {
      std::string td_;
      int nb, mx, n;
      ss >> td_ >> nb >> mx >> n;
      td.bags.assign(nb, {});
    } else if (t == "b") {
      int id, v;
      ss >> id;
      while (ss >> v) td.bags.at(id - 1).push_back(v - 1);
    } else {
      int a = std::stoi(t), b;
      ss >> b;
      td.edges.push_back({a - 1, b - 1});
    }
  }
  return td;
}

inline void writeTd(std::ostream& out, const TreeDecomp& td, int n) {
  size_t mx = 0;
  for (const auto& b : td.bags) mx = std::max(mx, b.size());
  out << "s td " << td.bags.size() << ' ' << mx << ' ' << n << '\n';
  for (size_t i = 0; i < td.bags.size(); ++i) {
    out << "b " << i + 1;
    for (int v : td.bags[i]) out << ' ' << v + 1;
    out << '\n';
  }
  for (auto [a, b] : td.edges) out << a + 1 << ' ' << b + 1 << '\n';
}
