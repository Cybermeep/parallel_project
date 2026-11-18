#ifndef HAM_HAMSEARCH_H
#define HAM_HAMSEARCH_H

#include <vector>

#include "common.h"
#include "graph.h"

class HamSearch {
public:
    HamSearch(const Graph& graph, bool cycle);

    bool trivialAnswer(SearchResult& result) const;
    std::vector<std::vector<int>> roots() const;
    std::vector<std::vector<int>> expand(std::vector<std::vector<int>> frontier, size_t target, int maxLevels,
                                         SharedControl& control) const;
    void solvePrefix(const std::vector<int>& prefix, SharedControl& control) const;
    bool solvePrefixBudget(const std::vector<int>& prefix, SharedControl& control, long long budget,
                           std::vector<long long>& leftover) const;

private:
    struct State {
        std::vector<char> visited;
        std::vector<int> avail;
        std::vector<int> path;
        std::vector<int> stamp;
        std::vector<int> queue;
        int stampCounter = 0;
        int closersAbove = 0;
    };

    const Graph& g;
    bool cycle;
    bool pruning;
    int start;
    int connectivityEvery;

    void initState(State& s, const std::vector<int>& prefix) const;
    void visit(State& s, int u) const;
    void unvisit(State& s) const;
    bool feasible(State& s, int previousTail) const;
    bool connectedRemainder(State& s) const;
    bool complete(const State& s) const;
    void candidates(const State& s, std::vector<int>& out) const;
};

SearchResult solveHamiltonian(const Graph& graph, bool cycle, int threads, double timeLimit, int frontierFactor);

#endif
