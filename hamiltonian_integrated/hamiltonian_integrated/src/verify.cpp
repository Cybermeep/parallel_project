#include "verify.h"

static bool isPermutation(int n, const std::vector<int>& path, std::string& reason) {
    if ((int)path.size() != n) {
        reason = "length " + std::to_string(path.size()) + " != " + std::to_string(n);
        return false;
    }
    std::vector<char> seen(n, 0);
    for (int v : path) {
        if (v < 0 || v >= n) { reason = "vertex out of range"; return false; }
        if (seen[v]) { reason = "vertex " + std::to_string(v + 1) + " repeated"; return false; }
        seen[v] = 1;
    }
    return true;
}

bool verifyHamiltonian(const Graph& graph, const std::vector<int>& path, bool cycle, std::string& reason) {
    if (!isPermutation(graph.V, path, reason)) return false;
    for (size_t i = 0; i + 1 < path.size(); i++)
        if (!graph.hasEdge(path[i], path[i + 1])) {
            reason = "missing edge " + std::to_string(path[i] + 1) + "-" + std::to_string(path[i + 1] + 1);
            return false;
        }
    if (cycle && graph.V > 1 && !graph.hasEdge(path.back(), path.front())) {
        reason = "missing closing edge";
        return false;
    }
    reason = "ok";
    return true;
}

long long pathCost(const WeightedGraph& weights, const std::vector<int>& path, bool closed) {
    long long total = 0;
    for (size_t i = 0; i + 1 < path.size(); i++) total += weights.at(path[i], path[i + 1]);
    if (closed && path.size() > 1) total += weights.at(path.back(), path.front());
    return total;
}

bool verifyTour(const WeightedGraph& weights, const std::vector<int>& tour, bool closed, long long claimedCost,
                std::string& reason) {
    if (!isPermutation(weights.n, tour, reason)) return false;
    long long actual = pathCost(weights, tour, closed);
    if (actual != claimedCost) {
        reason = "cost mismatch " + std::to_string(actual) + " vs " + std::to_string(claimedCost);
        return false;
    }
    reason = "ok";
    return true;
}
