# tw — Exact treewidth via lists of minimal separators (a reimplementation of Tamaki, SEA 2019)

A C++20 implementation of the method of H. Tamaki, "Computing treewidth via exact and heuristic lists of minimal separators"
(SEA 2019 / talk slides at Shonan Meeting 144, March 7, 2019).
The paper itself was not available, so the details were inferred from the slides (see "Guesses and deviations" below).

All of the code in this repository was implemented by Claude Opus 5.5.

## Build and run

```sh
make                 # build tw and td-validate
make test            # randomized tests (checked against brute-force enumeration and a brute-force treewidth DP)

./tw [-m ascend|descend|alternate] [-p0] [-v|-vv] [input.gr] > output.td
./td-validate input.gr output.td      # validate the tree decomposition ("OK width=W")
```

* Input: PACE 2017 `.gr` (read from standard input if no file is given); output: PACE 2017 `.td` (written to standard output)
* `-m`: search method. The default is `ascend` (a required feature); `descend` and `alternate` are also implemented
* `-p0`: disable preprocessing (the safe reduction rules) and solve with the algorithm of the slides alone
* `-v` / `-vv`: print progress to standard error as `c ...` lines

## Files

| File | Contents |
|---|---|
| `src/bits.hpp` | Variable-length bitset for vertex sets |
| `src/graph.hpp` | Graphs, connected components, tree decompositions, PACE I/O, validation and compression of tree decompositions |
| `src/msdp.hpp` | Part 1: `msDP(G, k, Δ)` |
| `src/minseps.hpp` | Part 2: `listExact(G, k)`; Part 3: `HeuristicLister` (listHeuristic); successors in Λ(G); full enumeration by Berry et al. (for testing) |
| `src/heuristics.hpp` | Initial tree decomposition by min-fill (upper bound), MMD+ lower bound |
| `src/reduce.hpp` | Preprocessing: removal of simplicial / almost simplicial vertices and restoration of the tree decomposition |
| `src/solver.hpp` | Ascend / Descend / Alternate, processing per connected component |
| `tests/test.cpp` | Randomized tests |
| `tools/validate.cpp` | `.td` validator |

## Algorithm overview

### Part 1: msDP(G, k, Δ) (slides 13–29)

We fix a total order 0 < 1 < … on the vertices and compare vertex sets by their minimum vertex. A connected set C is **inbound**
if G∖N(C) has a full component D with min D < min C.

1. L = { C : C is inbound and N(C) ∈ Δ } (for each S in Δ, the full components of G∖S except the one with the smallest minimum vertex)
2. Process L in increasing order of |C|, and mark C as feasible if `isFeasible(C)` is true
3. Finally, decide by `isFeasible(V(G))` (N(V) = ∅)

```
isFeasible(C):
  if |N[C]| ≤ k+1: return true                         # leaf bag N[C]
  for marked D (min C ∈ D, D ⊊ C, |N(C) ∪ N(D)| ≤ k+1):
     if allFeasible(N(C) ∪ N(D), C): return true
  return allFeasible(N(C) ∪ {min C}, C)
allFeasible(S, C):
  if |S| > k+1: return false
  for each component D of G∖S inside C:
     if N(D) = S, require isFeasible(D) recursively (Case 2)
     otherwise, require that D is marked (Case 1, inbound)
  return true
```

Potential maximal cliques are not enumerated; the bags S appear implicitly in the recursion. We record the S used on success
(or the leaf bag N[C]), and on YES we reconstruct the tree decomposition from these records (the subtrees for the components
of G∖S inside C are attached as children of the bag S).

### Part 2: listExact(G, k) (slides 33–44)

For a connected A (a ∈ A, b ∉ N[A], and N(A) is a minimal a-b separator with full components A and C_b(A)) and F ⊆ N(A), let

  Δ_ab(A, F) = { minimal a-b separators S with |S| ≤ k, F ⊆ S, A ⊆ C_a, C_b ∩ N[A] = ∅, |C_a| ≤ |C_b| }.

We enumerate it recursively by

  Δ_ab(A, F) = Δ_ab(A, F ∪ {v}) ∪ Δ_ab(A′, F)   (v ∈ N(A)∖F, and N(A′) is the minimal (A∪{v})-b separator close to A ∪ {v}).

The pruning rules are:

* If |F| > k, |A| > |C_b(A)|, or (|F| = k and N(A) ≠ F), the set is empty; if N(A) = F, output {F}
* The vertices of N(A)∖F adjacent to b must be in S (they are all added to F at once, without branching)
* In the second branch, the set is empty if A′ ∩ F ≠ ∅
* The slides allow any v ∈ N(A)∖F for branching; we choose a vertex with the most neighbors in C_b(A)
  (we compared 6 selection rules on 26 PACE 2017 instances, and this one was about 30% faster in total than choosing the smallest-numbered vertex)
* **Vertex-disjoint paths**: greedily extend vertex-disjoint paths from the vertices of N(A)∖F through V∖N[A]∖{b}. Since S∖F can cut
  at most k − |F| of these paths, all vertices of the shortest (#paths − (k−|F|)) paths lie on the a side. If the resulting lower bound
  lb = |A| + (the number of those vertices) exceeds |C_b(A)|, or 2·lb > n − |F|, the set is empty (the "15 vertices" example on slide 41)

**Divide and conquer** (slides 42–44): for a min-fill vertex v, let X = N(v). Every minimal separator S either
(i) crosses X (two non-adjacent vertices a, b of X lie in different full components of G∖S), or
(ii) is local to some component C of G∖N[v] (S = N(C) or S ⊆ N[C]).
For (i), we enumerate Δ_ab for each ordered pair (a, b) of non-adjacent vertices of X (using a = min(C_a ∩ X) for pruning).
For (ii), we recurse on the local graph G[N[C]] in which N(C) is made into a clique. Finally, we re-check that every output set is a minimal separator of G.

### Part 3: listHeuristic(G, T, k) (slide 45)

* Δ0: for each bag X of T and each component C of G∖X, the minimal separators close to C (⊆ N(C) ⊆ X) of size at most k
* Δ_{i+1} = Δ_i ∪ { successors in Λ(G) of the elements of Δ_i, of size at most k } (a successor is, for a full component C and v ∈ S, the separator close to C ∪ {v})

### Overall search (slide 4)

* **Ascend**: starting from the lower bound (MMD+), increase k by one and stop when `msDP(G, k, listExact(G, k))` returns YES
* **Descend**: starting from the min-fill tree decomposition (of width w), repeatedly try to improve it by running msDP with k = w−1 on Δ0, Δ1, ….
  When the list stops growing, check optimality with `msDP(G, w−1, listExact(G, w−1))` (if YES, keep improving)
* **Alternate**: alternate one Ascend step (verifying the lower bound) and one Descend step (extending the list by one level and running msDP)
  so that their cumulative running times stay equal. Stop when the lower and upper bounds meet

The graph is solved one connected component at a time. Any width up to the largest width obtained for the earlier components is acceptable,
so that width is used as a lower bound for the remaining components.

## Guesses and deviations (since the paper itself was not available)

All of the changes preserve correctness (whenever the answer is YES, a valid tree decomposition of width at most k is obtained).

1. **Line 2 of isFeasible**: the slides have "N(C) ∈ Δ and |N[C]| ≤ k+1", but we dropped the condition N(C) ∈ Δ (a single bag N[C] always gives a valid tree decomposition).
2. **Non-full components in allFeasible**: a component D is accepted even if it is not marked, as long as |N[D]| ≤ k+1 (for the same reason). The `assert D is inbound` of the slides is not checked; we only check whether D is marked.
3. **Memoization**: successes are memoized permanently, and failures are memoized until the next new set is marked.
4. **Base conditions of Δ_ab**: among the conditions on slide 35, those that we judged could miss separators if used as they are were replaced.
   * `min(C_b) ≠ b` is not used (C_b shrinks during the recursion, so testing it with the minimum of the current C_b(A) could lose solutions). Instead, we use only `a = min(C_a ∩ X)` on the a side.
   * The size condition (|V|−k)/2 was changed to (n−|F|)/2 (so that it also holds for solutions with |S| < k). The condition |A| + (|N(A)|−k) is subsumed by the lower bound from vertex-disjoint paths (it corresponds to paths of length 1).
   * The pairs (a, b) are enumerated as **ordered pairs** of non-adjacent vertices of X = N(v); some of the symmetry is removed by |C_a| ≤ |C_b|, and duplicates are removed with a hash set.
5. **Local graph**: for a component C of G∖N[v], the local graph is G[N[C]] with N(C) made into a clique (N(C) is a minimal separator with full components {v} and C, so the standard realization lemma applies).
6. **Definition of Δ0**: we interpret "minimal separators local to a bag" as "minimal separators contained in a bag X and close to a component outside X".
7. **Preprocessing (an addition not in the slides)**: safe reductions (Bodlaender–Koster) that remove simplicial vertices and almost simplicial vertices of degree ≤ the lower bound are enabled by default. They can be disabled with `-p0`.
8. **Choice of the branching vertex v**: determined experimentally, as described above.
9. **Resource balancing in Alternate**: the slides only say "some resource balancing", so we keep the cumulative running times equal.

## Tests

`make test` checks the following with different random seeds (a seed can be given as an argument: `./test-bin 42`).

* `listExact(G, k)` coincides exactly with the minimal separators of size at most k among all minimal separators enumerated by the algorithm of Berry et al. (n ≤ 24)
* For each of the three search methods, both with and without preprocessing (`-p0`), the output is a valid tree decomposition whose width equals the treewidth computed by a brute-force treewidth DP (n ≤ 14)
* The known treewidths of the Petersen graph, grids, cycles, complete graphs and edgeless graphs are obtained (including versions with randomly relabeled vertices)
* The three methods give the same results on random graphs with n = 20–35

## Experimental results (for reference)

Results on all 200 instances of the PACE 2017 exact track with a time limit of 30 seconds per instance
(2 runs in parallel on a 4-core machine; outputs were checked with `td-validate` and compared with the published optimal widths).

| Method | Solved within 30 s | Of which within 10 s | Wrong width | Invalid decomposition |
|---|---|---|---|---|
| ascend | 162 / 200 | 136 | 0 | 0 |
| alternate | 171 / 200 | 150 | 0 | 0 |

The first version (branching vertex = the smallest-numbered one, without the `compOf` optimization) solved 144 and 155 instances, respectively.

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for the full text.

The license covers the code and documentation in this repository. The underlying algorithm is due to the work of H. Tamaki cited above;
this repository is not an official implementation by Tamaki.
