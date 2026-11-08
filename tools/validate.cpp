// td-validate: checks that a .td file is a tree decomposition of a .gr graph.
// usage: td-validate graph.gr decomposition.td   -> prints "OK width=W" or an error
#include <fstream>
#include <iostream>

#include "../src/graph.hpp"

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: td-validate graph.gr decomposition.td\n";
    return 2;
  }
  std::ifstream g(argv[1]), t(argv[2]);
  if (!g || !t) {
    std::cerr << "cannot open input\n";
    return 2;
  }
  try {
    Graph G = readGr(g);
    TreeDecomp td = readTd(t);
    std::string err = validateTD(G, td);
    if (!err.empty()) {
      std::cout << "INVALID: " << err << "\n";
      return 1;
    }
    std::cout << "OK width=" << td.width() << "\n";
  } catch (const std::exception& e) {
    std::cout << "INVALID: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
