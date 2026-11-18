# Hamiltonian Path – integrated version (IT388)

Team: Jeff Penn, Alysa Solomon, Sam Mate

## Build

```
make            # hampath (OpenMP), run_bench, gen_synthetic
make mpi        # hampath_mpi with mpicxx (needs an MPI install, e.g. on Expanse)
make test       # smoke benchmark
make clean
```


## Running the solver

```
./hampath data/tsphcp/SH_64.hcp                          # Hamiltonian path (default for graphs)
./hampath data/tsphcp/SH_64.hcp --mode cycle --threads 8
./hampath data/graphdb/sam_test4.hcp --algo sam-recursive
./hampath data/graphdb/20_NH.g6 --mode cycle --batch --expect non-hamiltonian
./hampath data/graphdb/18_H.g6 --graph 12 --mode cycle --print-matrix --create-code
./hampath --visualize-demo
./hampath data/tsplib/ulysses22.tsp --expect 7013        # TSP branch and bound
./hampath data/tsplib/ulysses16.tsp --mode shp           # shortest Hamiltonian path
./hampath data/tsphcp/SH_64.hcp --mode cycle --algo hcp-as-tsp
mpirun -np 5 ./hampath_mpi data/tsphcp/GP0m_1006.hcp --mode cycle --threads 16
```

Modes:

- `path`: Hamiltonian path (the project's main problem).
- `cycle`: Hamiltonian cycle (what the published benchmark answers refer to).
- `tsp`: optimal tour on a TSPLIB file.
- `shp`: shortest Hamiltonian path, the modified TSP from the plan.

Algorithms for `path`/`cycle`:

- `pruned` (default): the parallel solver.
- `sam-iterative` / `sam-recursive`: Sam's baselines.
- `hcp-as-tsp`: runs the plan's transformations and solves with branch and bound.

Every run ends with one machine-readable line:

```
RESULT,instance,mode,algorithm,ranks,threads,vertices,edges,status,cost,expected,verdict,seconds,nodes,lower_bound,verified
```

`verified` comes from an independent checker (`src/verify.cpp`): every vertex appears exactly once, every edge exists, and the cost is recomputed. `verdict` compares the answer with `--expect`: `MATCH`, `MISMATCH`, `UNSOLVED` (time limit), or `GAP=x%` for an unproven TSP tour.

## parallel solver:

**Pruned DFS (`src/hamsearch.cpp`).** This is Sam's DFS made iterative, so 8,000-vertex graphs don't overflow the stack, with the following pruning added:

- Next vertex chosen by fewest unvisited neighbours (Warnsdorff ordering).
- A branch is dropped when an unvisited vertex can no longer get the 2 edges it needs (1 at a path end), when two vertices would both have to be the path's end, or when the unvisited vertices split into disconnected pieces.
- Cycle search starts at the lowest-degree vertex and only explores one direction of each cycle.

**OpenMP.** The search tree is expanded breadth-first until there are `frontier × threads` partial paths. The threads then take them with `#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)`. The graph is read-only and each thread keeps its own visited state. The stop flag, the best cost, the solution and the node counter are plain shared variables, changed only inside `#pragma omp critical`. Each thread syncs with them every 1,024 search nodes, so it notices another thread's answer within that many steps.

The OpenMP code uses only the IT 388 Test #2 syntax list:

- `parallel for`, `critical`, `num_threads`, `schedule(dynamic, chunk)` and a scalar `reduction(+ : …)` in batch mode
- `omp_get_wtime()` for all timing and `omp_get_num_procs()` for the default thread count
- no atomics, locks or other OpenMP features

**Branch and bound (`src/tsp.cpp`).**

- Starting upper bound: nearest-neighbour + 2-opt.
- Lower bound: Held-Karp 1-tree with subgradient multipliers at the root, then a minimum spanning tree over the unvisited vertices at every node.
- The best cost is shared through `#pragma omp critical`. Each thread keeps a copy that it refreshes every 1,024 nodes and whenever it finds a better tour.
- `shp` adds a dummy vertex joined to everything with weight 0 (transformation 3).
- `hcp-as-tsp` adds a universal vertex for path mode (transformation 1) and uses 0/1 weights (transformation 2): the graph is Hamiltonian exactly when the optimum is 0.

**MPI (`src/mpi_driver.cpp`).** This is the manager-worker queue from the Intellectual Challenge section. It uses only functions from the IT 388 formula sheet: `MPI_Init`, `MPI_Finalize`, `MPI_Comm_size`, `MPI_Comm_rank`, `MPI_Barrier`, `MPI_Wtime`, `MPI_Send`, `MPI_Recv` and `MPI_Reduce`. It does not use `MPI_Probe`, `MPI_Iprobe`, `MPI_Get_count` or `MPI_Init_thread`.

- Rank 0 builds a queue of partial paths. Each worker sends a request with `MPI_Send`, and rank 0 answers with `MPI_Recv(…, MPI_ANY_SOURCE, …)` and the next partial path, or a stop message.
- A request also carries the worker's result from its last task: a Hamiltonian path, or a better TSP tour. Rank 0 sends the current best cost back with each new task, so bounds spread between ranks every time a worker asks for work.
- Without `MPI_Iprobe`, a busy worker can't be interrupted. Instead, each task has a budget of 524,288 search nodes. When the budget runs out, the worker sends the unexplored part of its subtree back to rank 0, which puts it at the front of the queue. No work is repeated, and every worker reports to rank 0 at least that often. So the run stops shortly after any worker finds a path. Without the budget, a worker stuck in a huge subtree could keep the whole job running long after the answer was found.
- Variable-length data goes in two messages: a fixed-size header with the length, then the data. This avoids needing `MPI_Get_count`.
- Only the main thread of each rank makes MPI calls, always outside OpenMP parallel regions.
- In TSP mode the budget isn't needed. Branch-and-bound tasks always finish, and a better tour only tightens the bound.
- Node counts are summed with `MPI_Reduce`, and times are measured with `MPI_Wtime` between `MPI_Barrier` calls.

## Benchmarking against public sets

`bench/instances.csv` lists the public instances and their published answers:

| suite | instances | published answer | source |
|---|---|---|---|
| `fhcp` | 45 graphs, 64–8,544 vertices | all Hamiltonian | Baniasadi, Ejov, Haythorpe, Rossomakhine, *A new benchmark set for TSP and HCP*, arXiv:1806.09285 (Flinders) |
| `graphdb` | 39,635 cubic graphs on 18 vertices; 14,498 on 20 | 18_H all Hamiltonian, 20_NH all non-Hamiltonian | Flinders Hamiltonian Cycle Project graph database (GENREG, Meringer 1999) |
| `tsplib` | 30 classic instances | published optimal tour lengths | TSPLIB95 (Reinelt) |
| `smoke` | small mix of the above | | |
| `synthetic` | from `./gen_synthetic` | planted cycle (Hamiltonian by construction); random TSP must give the same optimum in every configuration | seeded |

Only burma14, ulysses16 and ulysses22 are included in `data/tsplib`. To add the others, download the `.tsp` files from the TSPLIB95 site (Heidelberg) into `data/tsplib/` with the names used in `bench/instances.csv`. The runner skips any that are missing.

There are two runners that take the same options. `bench/run_bench.py` (Python 3.6+) is the one used for Expanse.

- S(p) = T(1) / T(p), with p = MPI ranks × OpenMP threads
- E(p) = S(p) / p
- the serial fraction f, solved from Amdahl's law with the measured speedup: f = (1/S − 1/p) / (1 − 1/p), plus the matching Amdahl limit 1 / (f + (1 − f)/p)

It also draws S(p) and E(p) plots when matplotlib is available.

`./run_bench` is the C++ version of the same runner, without the formulas and plots.

```
python3 bench/run_bench.py --suite smoke --threads 1,2,4 --plots results/plots
python3 bench/run_bench.py --suite graphdb,fhcp,tsplib --threads 1,2,4,8,16 --repeat 3 --time-limit 120 --out results/public.csv
python3 bench/run_bench.py --suite tsplib --ranks 3,5,9 --threads 8 --launcher "mpirun -np {ranks} --map-by slot:PE={threads} --bind-to core" --baseline results/public.csv --out results/mpi.csv
python3 bench/run_bench.py --suite tsplib --fetch-tsplib
```




The C++ runner:

```
./run_bench --suite smoke --threads 1,2,4
./run_bench --suite graphdb,fhcp,tsplib --threads 1,2,4,8,16 --repeat 3 --time-limit 120 --out results/public.csv
./run_bench --suite graphdb --algos sam-iterative,sam-recursive,pruned,hcp-as-tsp --threads 1 --out results/algorithms.csv
./gen_synthetic && ./run_bench --manifest bench/synthetic.csv --suite synthetic --threads 1,2,4,8
./run_bench --suite fhcp,tsplib --ranks 3,5,9 --threads 8 --launcher "mpirun -np {ranks}"
```

Each run writes a CSV and a Markdown summary. The CSV has the status, the verdict against the published answer, median time, search nodes, speedup against 1 core for the same algorithm, and efficiency. The summary has a correctness table per suite and a check that all configurations agree on every proven optimum.


## Files

```
src/common.*        status codes, omp_get_wtime timer, shared state guarded by omp critical
src/graph.*         Graph (Sam) + WeightedGraph
src/io.*            .hcp/.tsp/.txt/.asc/.g6/.tour readers
src/sam_dfs.*       Sam's iterative and recursive DFS baselines
src/hamsearch.*     pruned DFS + OpenMP work distribution
src/tsp.*           branch and bound, Held-Karp bound, plan transformations
src/mpi_driver.*    MPI manager-worker (compiled with -DUSE_MPI)
src/verify.*        independent answer checker
src/visualize.*     Alysa's adjacency-matrix class, print() and createCode()
src/main.cpp        command line
bench/run_bench.py  benchmark runner with S(p), E(p) and Amdahl serial fraction (used by sbatch.sb)
bench/run_bench.cpp the same runner in C++
bench/gen_synthetic.cpp  seeded synthetic instance generator
bench/instances.csv public instances with published answers
data/               benchmark data from Sam's TestData plus three TSPLIB files
sbatch.sb           Expanse Slurm job: build, then run every benchmark suite
```
