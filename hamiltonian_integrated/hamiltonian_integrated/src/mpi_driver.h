#ifndef HAM_MPI_DRIVER_H
#define HAM_MPI_DRIVER_H

#ifdef USE_MPI

#include "common.h"
#include "graph.h"

SearchResult mpiSolveHamiltonian(const Graph& graph, bool cycle, int threads, double timeLimit, int frontierFactor);
SearchResult mpiSolveTsp(const WeightedGraph& weights, int threads, double timeLimit, int frontierFactor,
                         long long knownLowerBound);

#endif

#endif
