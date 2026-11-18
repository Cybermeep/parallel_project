#ifndef HAM_VERIFY_H
#define HAM_VERIFY_H

#include <string>
#include <vector>

#include "graph.h"

bool verifyHamiltonian(const Graph& graph, const std::vector<int>& path, bool cycle, std::string& reason);
bool verifyTour(const WeightedGraph& weights, const std::vector<int>& tour, bool closed, long long claimedCost,
                std::string& reason);
long long pathCost(const WeightedGraph& weights, const std::vector<int>& path, bool closed);

#endif
