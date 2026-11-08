CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra -march=native
HDRS := $(wildcard src/*.hpp)

all: tw td-validate

tw: src/main.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ src/main.cpp

td-validate: tools/validate.cpp src/graph.hpp src/bits.hpp
	$(CXX) $(CXXFLAGS) -o $@ tools/validate.cpp

test-bin: tests/test.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -g -o $@ tests/test.cpp

test: test-bin
	./test-bin

clean:
	rm -f tw td-validate test-bin

.PHONY: all test clean
