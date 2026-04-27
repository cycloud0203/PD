# Fixed-Outline Floorplanning with B\*-Tree and Simulated Annealing

## 1. Problem Description

Given a set of rectangular hard blocks, a netlist, and a fixed chip outline, find a non-overlapping placement of all blocks within the outline that minimizes the cost function:

```
cost = alpha * area + (1 - alpha) * wirelength
```

where `alpha` is a user-specified weight.

## 2. Overall Algorithm

Our floorplanner follows a three-stage pipeline:

1. **Multi-start sampling** -- quickly evaluate many random initial B\*-tree configurations to pick the best starting points.
2. **Multi-threaded Simulated Annealing (SA)** -- launch parallel SA workers from the top-K seeds, each optimizing the cost function under the outline constraint.
3. **Post-refinement** -- apply steepest-descent local search (swap + rotation enumeration) and a global shift refinement to further reduce wirelength.

## 3. Data Structures

### 3.1 B\*-Tree Representation

The B\*-tree is a binary tree that compactly encodes a non-overlapping, compacted floorplan. Each internal/leaf node represents one block. The tree structure determines the placement:

- **Left child** of node `v`: the block placed directly to the *right* of `v` (horizontally adjacent).
- **Right child** of node `v`: the block placed directly *above* `v` (vertically adjacent at the same x-coordinate as `v`).

Each tree node stores:

```cpp
struct BTreeNode {
    int parent;
    int left;
    int right;
    int blockId;
    bool rotated;
};
```

The packing procedure traverses the tree in a DFS order and maintains a **horizontal contour** (a staircase-shaped skyline) to determine the y-coordinate of each block in O(n) time per packing.

### 3.2 Perturbation Moves

Five types of SA perturbation moves are supported, each with a configurable weight:

| Move | Description |
|---|---|
| **Rotate** | Flip one block's width/height. |
| **Swap** | Exchange the block IDs (and rotation flags) of two tree nodes. |
| **Delete-Insert** | Remove a node from the tree and re-insert it at a random (or net-aware) position. |
| **Rebuild** | Completely re-randomize the entire tree. |
| **Swap-Children** | Swap the left and right children of one node. |

The move mix is tuned per circuit size:

```cpp
static PerturbConfig forCircuit(int n, int numNets) {
    double density = (numNets > 0 && n > 0) ? (double)numNets / n : 0;
    if (n <= 15 && density > 10.0)
        return {8, 2, 8, 0, 4, 90};   // heavy net-aware delete-insert
    if (n <= 12) return {6, 5, 5, 1, 5, 70};
    if (n <= 15) return {7, 6, 6, 1, 0, 0};
    if (n <= 50) return {6, 4, 8, 2, 3, 50};
    return {5, 2, 11, 2, 2, 50};
}
```

### 3.3 Net-Aware Perturbation

For the delete-insert move, blocks can be re-inserted *near* a net-adjacent neighbor with a configurable probability (`netAwarePercent`). This exploits locality: connected blocks placed nearby tend to produce lower wirelength.

The adjacency list `blockAdj` is pre-built and sorted by connectivity weight (most-connected neighbors first):

```cpp
// After building blockAdj, sort by descending edge count
for (auto& adj : _blockAdj) {
    std::sort(adj.begin(), adj.end());
    std::vector<std::pair<int,int>> weighted;
    // ... count duplicate edges as weight ...
    std::sort(weighted.begin(), weighted.end());
    adj.clear();
    for (auto& [w, id] : weighted) adj.push_back(id);
}
```

### 3.4 Incremental Wirelength Cache (`WLCache`)

Full HPWL recomputation costs O(total pins) per move. Our `WLCache` reduces this by tracking which blocks actually moved and only recomputing the nets that are affected:

```cpp
class WLCache {
    // Per-block: reference center coordinates
    std::vector<int> _refCX2, _refCY2;
    // Per-net: cached HPWL contribution
    std::vector<long long> _netWL;
    // Dirty flags for incremental update
    std::vector<char> _netDirty;
    // ...
};
```

On each call to `compute()`:
1. Detect which blocks moved (center changed).
2. Mark their incident nets as dirty.
3. Recompute only the dirty nets and update the total.
4. On `accept()`, commit the new values; on `undoPerturb()`, the old values are still valid.

This reduces wirelength evaluation from O(total pins) to O(affected pins), a large speedup for circuits with many blocks.

### 3.5 Thread-Safe Global Best (`GlobalBest`)

All SA threads share a single `GlobalBest` structure protected by a mutex:

```cpp
struct GlobalBest {
    std::mutex mtx;
    double bestCost = 1e18;
    bool feasible = false;
    std::vector<BTreeNode> nodes;
    int root = -1;
    int width = 0, height = 0;
    std::vector<int> bestX, bestY, bestX2, bestY2;

    bool tryUpdate(double cost, bool feas, ...) {
        std::lock_guard<std::mutex> lock(mtx);
        bool doUpdate = false;
        if (feas && !feasible) doUpdate = true;
        else if (feas == feasible && cost < bestCost) doUpdate = true;
        if (doUpdate) { /* copy new solution */ }
        return doUpdate;
    }
};
```

Feasible solutions always dominate infeasible ones, regardless of cost.

## 4. Multi-Threading Strategy

The floorplanner launches multiple SA threads in parallel:

```cpp
unsigned hwThreads = std::thread::hardware_concurrency();
// ...
int numThreads = std::min((int)hwThreads, maxThreads);

std::vector<std::thread> workers;
for (int t = 0; t < numThreads; ++t) {
    workers.emplace_back([...]() {
        SAEngine engine(_sd, _gb, baseSeed + (unsigned)t, t, threadAlpha);
        engine.run(deadline);  // or runFromState(...)
    });
}
for (auto& w : workers) w.join();
```

Key design choices:

- **Top-K seeding**: the best K initial configurations from multi-start sampling are assigned to the first K threads, so the best starting points get the most SA budget.
- **Multi-alpha exploration**: remaining threads run with different `alpha` values (0.1, 0.3, 0.5, 0.7, 0.9) to explore diverse trade-offs. This helps escape local optima that a single alpha would get stuck in.
- **Each thread has its own `BTree`, `WLCache`, and `std::mt19937` RNG** -- no contention except the lightweight `GlobalBest::tryUpdate()` mutex.

## 5. Simulated Annealing Details

Each SA thread runs a two-phase schedule:

### Phase 1: Full SA (area + wirelength + outline penalty)

- **Initial temperature** computed adaptively: sample 500 random moves and set T so that 95% of uphill moves are accepted initially.
- **Adaptive cooling**: the cooling rate depends on the acceptance ratio:
  - `> 80%` acceptance: fast cooling (0.85x)
  - `15%--80%` acceptance: base cooling (auto-tuned to fill the time budget)
  - `< 15%` acceptance: slow cooling
- **Outline penalty**: solutions exceeding the outline are penalized proportionally to overflow, not rejected outright, so SA can traverse through infeasible space.

### Phase 2: Wirelength Refinement SA

After Phase 1, the thread reloads its local best feasible solution and runs a lower-temperature SA that only accepts feasible moves, focusing on minimizing the cost within the outline.

## 6. Post-Refinement

### 6.1 Steepest-Descent Local Search

After all SA threads finish, the main thread loads the global best solution and applies a deterministic steepest-descent search:

- For small circuits (n <= 12): exhaustively enumerate all 2^n rotation masks combined with all O(n^2) pairwise block-ID swaps, greedily applying the best-improving move until convergence.
- For larger circuits: greedy single-swap and single-rotation hill climbing.

### 6.2 Global Shift Refinement

As a final step, the entire placement is shifted within the outline to minimize wirelength. The optimal (dx, dy) offset is found by ternary search along each axis:

```cpp
void Floorplanner::applyGlobalShiftRefinement() {
    // Find feasible shift ranges [loX, hiX] and [loY, hiY]
    int bestDx = optimizeAxisShift(..., loX, hiX, true);
    int bestDy = optimizeAxisShift(..., loY, hiY, false);
    // Apply if wirelength improves
    if (newWL < oldWL) {
        for (int i = 0; i < n; ++i) {
            _gb.bestX[i] += bestDx;  _gb.bestY[i] += bestDy;
            _gb.bestX2[i] += bestDx; _gb.bestY2[i] += bestDy;
        }
    }
}
```

This is especially helpful for circuits with many fixed terminals (I/O pads).

## 7. Experimental Results

All benchmarks run with `alpha = 0.5`. Results from the evaluator:

| Benchmark | Area | Wirelength | Runtime (s) | Normalized Cost | Score |
|---|---|---|---|---|---|
| ami33 | 1,252,440 | 73,166 | 132.88 | 0.0870 | 7.58 |
| ami49 | 38,100,832 | 810,133 | 180.84 | 0.1832 | 8.02 |
| apte | 47,313,280 | 698,659 | 42.16 | 0.0271 | 8.49 |
| hp | 9,852,528 | 201,906 | 37.14 | 0.0197 | 8.48 |
| xerox | 20,424,474 | 495,778 | 42.15 | 0.0324 | 8.50 |

Score = 0.8 * Quality Score + 0.2 * Runtime Score.

### Floorplan Visualizations

#### ami33

![ami33 floorplan](output/viz/ami33.png)

#### ami49

![ami49 floorplan](output/viz/ami49.png)

#### apte

![apte floorplan](output/viz/apte.png)

#### hp

![hp floorplan](output/viz/hp.png)

#### xerox

![xerox floorplan](output/viz/xerox.png)

## 8. Summary of Key Techniques

1. **B\*-tree** provides a compact, efficient, and easy-to-perturb representation for non-overlapping placements with O(n) packing.
2. **Simulated Annealing** with adaptive cooling, outline-overflow penalty, and multi-phase scheduling (packing then wirelength refinement) effectively explores the solution space.
3. **Multi-threading** with diverse initial seeds and alpha values exploits modern multi-core CPUs, achieving solution diversity and better global optima.
4. **Incremental wirelength cache** (`WLCache`) avoids redundant HPWL recomputation, providing a significant speedup per SA step.
5. **Net-aware perturbation** biases move selection toward placing connected blocks nearby, improving wirelength.
6. **Post-refinement** via steepest-descent local search and global shift ternary search squeezes out the last few percent of cost reduction.

## 9. Source File Overview

| File | Description |
|---|---|
| `src/main.cpp` | Entry point; parses arguments, invokes floorplanner, writes output. |
| `src/floorplanner.h/cpp` | Top-level orchestrator: parsing, normalization, multi-thread launch, post-refinement. |
| `src/btree.h/cpp` | B\*-tree data structure: packing, perturbation, undo, initialization strategies. |
| `src/sa_engine.h/cpp` | Simulated Annealing engine: temperature schedule, SA loop, WL refinement phase. |
| `src/wl_cache.h` | Incremental HPWL cache for fast wirelength evaluation. |
| `src/global_best.h` | Thread-safe global best solution tracker. |
| `src/shared_data.h` | Read-only shared data (block dimensions, nets, adjacency) for all threads. |
| `src/module.h/cpp` | Block, Terminal, and Net class definitions. |
