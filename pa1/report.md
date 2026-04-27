# PA1 Report: Data Structures and Findings

## 1. Program overview

This program implements a hypergraph partitioner based on the Fiduccia-Mattheyses (FM) algorithm.  
Input is parsed as a set of nets and cells, and the solver searches for a legal bi-partition that minimizes cutsize under the given balance factor.

The implementation extends a standard FM flow with:
- multiple initialization strategies,
- iterative local search (ILS) perturbation and refinement,
- random restarts under a time budget.

## 2. Data structures used

### 2.1 Cell structure

Each cell is represented by:
- `gain`: current FM gain.
- `part`: current partition (`0` or `1`).
- `locked`: whether the cell has already been moved in current FM pass.
- `bPrev`, `bNext`: links for an intrusive doubly-linked bucket list.

Why this design:
- storing bucket links inside each cell avoids extra node allocation,
- insert/remove in bucket list is O(1),
- good cache locality during gain updates.

### 2.2 Net structure

Each net stores:
- `dist[2]`: number of incident cells in partition 0 and partition 1.

Why this design:
- FM gain updates only need partition distribution counts,
- cutsize check is direct (`dist[0] > 0 && dist[1] > 0`),
- memory footprint is small and update cost is low.

### 2.3 Adjacency lists

Main containers:
- `_cellNets`: list of nets connected to each cell.
- `_netCells`: list of cells connected to each net.
- `_cellName`: cell id to name mapping.

Why this design:
- bidirectional traversal is needed by FM:
  - from moved cell to affected nets,
  - from affected net to neighboring cells.
- vectors provide compact storage and fast iteration.

### 2.4 Bucket lists

Per partition:
- `_bHead[2]`: array of bucket heads indexed by shifted gain.
- `_maxGIdx[2]`: current highest non-empty gain bucket.
- gain index conversion:
  - `gainIdx(gain) = gain + _maxPinNum`
  - bucket size `= 2 * _maxPinNum + 1`

Why this design:
- max-gain candidate extraction is near O(1) amortized,
- supports rapid dynamic gain updates after each move.

### 2.5 Move history and best solution

- `_moveStack`: stores move order in one FM pass, used for rollback.
- `_bestPart`: stores globally best partition found.
- `_bestCutSize`: best cutsize value.

Why this design:
- FM pass needs rollback after best prefix gain point,
- global best tracking allows multi-start and ILS frameworks.

## 3. Algorithm flow

1. Parse input and build hypergraph.
2. Remove single-cell nets (cannot contribute to cutsize).
3. Compute balance bounds from `bFactor`.
4. Run multi-phase search within time budget:
   - deterministic initial partition + FM,
   - BFS-based initial partitions + FM,
   - ILS perturbation around best solution + FM,
   - random restarts + FM.
5. Output global best partition.

Inside each FM pass:
- initialize net distribution and all cell gains,
- repeatedly choose legal highest-gain movable cell,
- update neighboring gains incrementally,
- keep best prefix by accumulated gain,
- rollback moves after best prefix.

## 4. Findings from this assignment

### 4.1 Effective implementation choices

1. **Intrusive bucket list was critical**  
   O(1) remove/insert and direct max-bucket tracking significantly reduced pass overhead compared with searching all free cells.

2. **Incremental gain update is essential**  
   Updating only cells on touched nets made each move proportional to local connectivity, avoiding expensive full recomputation.

3. **Multi-start strategy improved stability**  
   Deterministic start alone can be trapped in local minima. BFS starts, perturbation, and random restarts improved best-case and average quality.

4. **Time-budgeted search gives practical control**  
   The solver scales search effort with problem size and avoids unbounded runtime.

### 4.2 Practical trade-offs observed

- More aggressive perturbation helps escape local minima, but too large perturbation destroys good structure and slows convergence.
- Early stopping in FM pass reduces wasted moves when no improvement is likely, but if too aggressive it can miss delayed improvements.
- Sorting adjacency lists slightly improves locality and iteration behavior, especially in large instances.

### 4.3 Lessons learned

- Data-structure design has direct impact on algorithm performance in FM.
- A strong local optimizer (FM pass) still needs a global exploration strategy for robust final quality.
- Keeping implementation details simple (vectors + intrusive links + counters) can provide both speed and maintainability.

## 5. Conclusion

This implementation combines an efficient FM core with lightweight metaheuristics (BFS starts, ILS perturbation, random restarts).  
The chosen data structures support fast incremental updates and legal move selection, enabling better cutsize results under practical runtime constraints.
