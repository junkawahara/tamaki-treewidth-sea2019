// tw: exact treewidth via exact and heuristic lists of minimal separators
// (H. Tamaki, SEA 2019).  Reads a PACE 2017 .gr graph, writes a .td file.
//
// usage: tw [-m ascend|descend|alternate] [-p0] [-v[v]] [input.gr] [> output.td]
#include <fstream>
#include <iostream>
#include <string>

#include "solver.hpp"

static void usage() {
  std::cerr << "usage: tw [-m ascend|descend|alternate] [-p0] [-v] [-vv] [input.gr]\n"
               "  -p0: disable the safe reductions (preprocessing)\n"
               "  reads PACE 2017 .gr from the file (or stdin), writes .td to stdout\n";
}

int main(int argc, char** argv) {
  std::ios::sync_with_stdio(false);
  Options opt;
  std::string path;
  bool havePath = false;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "-m" && i + 1 < argc) {
      opt.mode = argv[++i];
      if (opt.mode != "ascend" && opt.mode != "descend" && opt.mode != "alternate") {
        usage();
        return 2;
      }
    } else if (a == "-p0") {
      opt.reduce = false;
    } else if (a == "-v") {
      opt.verbose = 1;
    } else if (a == "-vv") {
      opt.verbose = 2;
    } else if (a == "-h" || a == "--help") {
      usage();
      return 0;
    } else if ((!a.empty() && a[0] == '-') || havePath) {  // unknown option or a second input file
      usage();
      return 2;
    } else {
      path = a;
      havePath = true;
    }
  }
  Graph G;
  try {
    if (!havePath) {
      G = readGr(std::cin);
    } else {
      std::ifstream in(path);
      if (!in) {
        std::cerr << "cannot open " << path << "\n";
        return 1;
      }
      G = readGr(in);
    }
  } catch (const std::exception& e) {
    std::cerr << "input error: " << e.what() << "\n";
    return 1;
  }
  Solver solver(opt);
  TreeDecomp td = solver.solve(G);
  if (opt.verbose) std::cerr << "c treewidth " << td.width() << std::endl;
  writeTd(std::cout, td, G.n);
  // A failed write (e.g. a full disk) may only show up when the buffer is flushed.
  std::cout.flush();
  if (!std::cout) {
    std::cerr << "write error: could not write the tree decomposition\n";
    return 1;
  }
  return 0;
}
