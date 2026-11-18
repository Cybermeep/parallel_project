#ifndef HAM_TSP_H
#define HAM_TSP_H

#include <vector>

#include "common.h"
#include "graph.h"

class TspBranchAndBound {
public:
    explicit TspBranchAndBound(const WeightedGraph& weights);

    long long heuristicTour(std::vector<int>& tour, int maxStarts) const;
    double heldKarpRoot(long long upperBound, int iterations);
    double rootBound() const { return rootLowerBound; }
    long long tourCost(const std::vector<int>& tour) const;

    std::vector<std::vector<int>> expand(std::vector<std::vector<int>> frontier, size_t target, int maxLevels,
                                         SharedControl& control) const;
    void solvePrefix(const std::vector<int>& prefix, SharedControl& control) const;
    double prefixBound(const std::vector<int>& prefix) const;
    void setProofTarget(long long target) { proofTarget = target; }
    long long getProofTarget() const { return proofTarget; }
    const std::vector<double>& multipliers() const { return pi; }
    void setMultipliers(const std::vector<double>& values, double lowerBound) { pi = values; rootLowerBound = lowerBound; }

private:
    struct State {
        std::vector<char> visited;
        std::vector<int> path;
        std::vector<int> unvisited;
        std::vector<double> key;
        std::vector<char> inTree;
        long long cost = 0;
        double piRemaining = 0.0;
        long long local = 0;
        long long bestSeen = 0;
        bool stopped = false;
    };

    const WeightedGraph& wg;
    int n;
    std::vector<double> pi;
    double rootLowerBound;
    long long proofTarget;

    double reduced(int i, int j) const { return (double)wg.at(i, j) + pi[i] + pi[j]; }
    void initState(State& s, const std::vector<int>& prefix) const;
    double bound(State& s) const;
    void recurse(State& s, SharedControl& control) const;
    void twoOpt(std::vector<int>& tour) const;
};

SearchResult solveTsp(const WeightedGraph& weights, int threads, double timeLimit, int frontierFactor,
                      long long knownLowerBound);

WeightedGraph addZeroDummy(const WeightedGraph& weights);
WeightedGraph hamiltonicityToTsp(const Graph& graph);
Graph pathToCycleGraph(const Graph& graph);
std::vector<int> removeVertexFromCycle(const std::vector<int>& cycle, int vertex);

#endif
