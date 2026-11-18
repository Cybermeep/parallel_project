#include "tsp.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>


namespace {
const double INF = std::numeric_limits<double>::infinity();
long long ceilBound(double value) { return (long long)std::ceil(value - 1e-6); }
}

TspBranchAndBound::TspBranchAndBound(const WeightedGraph& weights)
    : wg(weights), n(weights.n), pi(weights.n, 0.0), rootLowerBound(0.0),
      proofTarget(std::numeric_limits<long long>::min()) {}

long long TspBranchAndBound::tourCost(const std::vector<int>& tour) const {
    long long total = 0;
    for (size_t i = 0; i < tour.size(); i++) total += wg.at(tour[i], tour[(i + 1) % tour.size()]);
    return total;
}

void TspBranchAndBound::twoOpt(std::vector<int>& t) const {
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 0; i < n - 2; i++)
            for (int j = i + 2; j < n; j++) {
                if (i == 0 && j == n - 1) continue;
                int a = t[i], b = t[i + 1], c = t[j], d = t[(j + 1) % n];
                long long delta = wg.at(a, c) + wg.at(b, d) - wg.at(a, b) - wg.at(c, d);
                if (delta < 0) {
                    std::reverse(t.begin() + i + 1, t.begin() + j + 1);
                    improved = true;
                }
            }
    }
}

long long TspBranchAndBound::heuristicTour(std::vector<int>& best, int maxStarts) const {
    long long bestCost = std::numeric_limits<long long>::max();
    int starts = std::min(n, std::max(1, maxStarts));
    for (int s = 0; s < starts; s++) {
        std::vector<int> tour = {s};
        std::vector<char> used(n, 0);
        used[s] = 1;
        for (int k = 1; k < n; k++) {
            int last = tour.back(), choice = -1;
            for (int v = 0; v < n; v++)
                if (!used[v] && (choice < 0 || wg.at(last, v) < wg.at(last, choice))) choice = v;
            used[choice] = 1;
            tour.push_back(choice);
        }
        if (n >= 4) twoOpt(tour);
        long long c = tourCost(tour);
        if (c < bestCost) {
            bestCost = c;
            best = tour;
        }
    }
    auto zero = std::find(best.begin(), best.end(), 0);
    std::rotate(best.begin(), zero, best.end());
    return bestCost;
}

double TspBranchAndBound::heldKarpRoot(long long upperBound, int iterations) {
    if (n < 3) return rootLowerBound = 0.0;
    std::vector<double> current(n, 0.0), bestPi(n, 0.0), key(n);
    std::vector<int> parent(n), degree(n);
    std::vector<char> inTree(n);
    double bestBound = -INF, lambda = 2.0;
    int stall = 0;
    for (int it = 0; it < iterations; it++) {
        auto red = [&](int i, int j) { return (double)wg.at(i, j) + current[i] + current[j]; };
        std::fill(degree.begin(), degree.end(), 0);
        std::fill(inTree.begin(), inTree.end(), 0);
        std::fill(key.begin(), key.end(), INF);
        key[1] = 0.0;
        parent[1] = -1;
        double total = 0.0;
        for (int k = 1; k < n; k++) {
            int u = -1;
            for (int v = 1; v < n; v++)
                if (!inTree[v] && (u < 0 || key[v] < key[u])) u = v;
            inTree[u] = 1;
            total += key[u];
            if (parent[u] >= 0) { degree[u]++; degree[parent[u]]++; }
            for (int v = 1; v < n; v++)
                if (!inTree[v] && red(u, v) < key[v]) { key[v] = red(u, v); parent[v] = u; }
        }
        int e1 = -1, e2 = -1;
        for (int v = 1; v < n; v++) {
            if (e1 < 0 || red(0, v) < red(0, e1)) { e2 = e1; e1 = v; }
            else if (e2 < 0 || red(0, v) < red(0, e2)) e2 = v;
        }
        total += red(0, e1) + red(0, e2);
        degree[0] = 2;
        degree[e1]++;
        degree[e2]++;
        double bound = total - 2.0 * std::accumulate(current.begin(), current.end(), 0.0);
        if (bound > bestBound + 1e-9) { bestBound = bound; bestPi = current; stall = 0; }
        else if (++stall >= std::max(5, n / 4)) { lambda *= 0.5; stall = 0; }
        if (ceilBound(bestBound) >= upperBound) break;
        double norm = 0.0;
        for (int v = 0; v < n; v++) norm += (double)(degree[v] - 2) * (degree[v] - 2);
        if (norm == 0.0 || lambda < 1e-5) break;
        double step = lambda * ((double)upperBound - bound) / norm;
        if (step <= 0) step = lambda / norm;
        for (int v = 0; v < n; v++) current[v] += step * (degree[v] - 2);
    }
    pi = bestPi;
    rootLowerBound = std::max(0.0, bestBound);
    return rootLowerBound;
}

void TspBranchAndBound::initState(State& s, const std::vector<int>& prefix) const {
    s.visited.assign(n, 0);
    s.key.assign(n, INF);
    s.inTree.assign(n, 0);
    s.path.clear();
    s.cost = 0;
    s.piRemaining = std::accumulate(pi.begin(), pi.end(), 0.0);
    for (int v : prefix) {
        if (!s.path.empty()) s.cost += wg.at(s.path.back(), v);
        s.path.push_back(v);
        s.visited[v] = 1;
        s.piRemaining -= pi[v];
    }
}

double TspBranchAndBound::bound(State& s) const {
    int start = s.path.front(), tail = s.path.back();
    s.unvisited.clear();
    for (int v = 0; v < n; v++)
        if (!s.visited[v]) s.unvisited.push_back(v);
    if (s.unvisited.empty()) return (double)(s.cost + wg.at(tail, start));
    double tree = 0.0;
    int first = s.unvisited[0];
    for (int v : s.unvisited) { s.key[v] = INF; s.inTree[v] = 0; }
    s.key[first] = 0.0;
    for (size_t k = 0; k < s.unvisited.size(); k++) {
        int u = -1;
        for (int v : s.unvisited)
            if (!s.inTree[v] && (u < 0 || s.key[v] < s.key[u])) u = v;
        s.inTree[u] = 1;
        tree += s.key[u];
        for (int v : s.unvisited)
            if (!s.inTree[v]) {
                double r = reduced(u, v);
                if (r < s.key[v]) s.key[v] = r;
            }
    }
    double tailMin = INF, startMin = INF, startSecond = INF;
    for (int v : s.unvisited) {
        tailMin = std::min(tailMin, reduced(tail, v));
        double r = reduced(start, v);
        if (r < startMin) { startSecond = startMin; startMin = r; }
        else if (r < startSecond) startSecond = r;
    }
    double connect = start == tail ? startMin + (s.unvisited.size() > 1 ? startSecond : startMin) : tailMin + startMin;
    return (double)s.cost + tree + connect - pi[tail] - pi[start] - 2.0 * s.piRemaining;
}

double TspBranchAndBound::prefixBound(const std::vector<int>& prefix) const {
    State s;
    initState(s, prefix);
    return bound(s);
}

void TspBranchAndBound::recurse(State& s, SharedControl& control) const {
    if (s.stopped) return;
    s.local++;
    if ((s.local & 1023) == 0 && control.checkpoint(1024, s.bestSeen)) {
        s.stopped = true;
        return;
    }
    int start = s.path.front(), tail = s.path.back();
    if ((int)s.path.size() == n) {
        long long total = s.cost + wg.at(tail, start);
        if (total < s.bestSeen && control.offerTour(total, s.path, proofTarget, s.bestSeen)) s.stopped = true;
        return;
    }
    std::vector<int> children;
    for (int v = 0; v < n; v++)
        if (!s.visited[v]) children.push_back(v);
    std::sort(children.begin(), children.end(), [&](int a, int b) { return wg.at(tail, a) < wg.at(tail, b); });
    for (int u : children) {
        long long edge = wg.at(tail, u);
        if (s.cost + edge >= s.bestSeen) continue;
        s.cost += edge;
        s.path.push_back(u);
        s.visited[u] = 1;
        s.piRemaining -= pi[u];
        double lb = bound(s);
        if (ceilBound(lb) < s.bestSeen) recurse(s, control);
        s.piRemaining += pi[u];
        s.visited[u] = 0;
        s.path.pop_back();
        s.cost -= edge;
        if (s.stopped) return;
    }
}

std::vector<std::vector<int>> TspBranchAndBound::expand(std::vector<std::vector<int>> frontier, size_t target,
                                                        int maxLevels, SharedControl& control) const {
    State s;
    int levels = 0;
    while (!frontier.empty() && frontier.size() < target && levels < maxLevels && !control.stop) {
        std::vector<std::pair<double, std::vector<int>>> next;
        for (const auto& prefix : frontier) {
            if ((int)prefix.size() >= n - 1) {
                next.push_back({prefixBound(prefix), prefix});
                continue;
            }
            initState(s, prefix);
            for (int u = 0; u < n; u++) {
                if (s.visited[u]) continue;
                std::vector<int> child = prefix;
                child.push_back(u);
                control.nodes++;
                double lb = prefixBound(child);
                if (ceilBound(lb) < control.best) next.push_back({lb, child});
            }
        }
        std::stable_sort(next.begin(), next.end(),
                         [](const auto& a, const auto& b) { return a.first < b.first; });
        frontier.clear();
        for (auto& item : next) frontier.push_back(std::move(item.second));
        levels++;
    }
    return frontier;
}

void TspBranchAndBound::solvePrefix(const std::vector<int>& prefix, SharedControl& control) const {
    State s;
    initState(s, prefix);
    s.stopped = control.checkpoint(0, s.bestSeen);
    if (s.stopped || ceilBound(bound(s)) >= s.bestSeen) return;
    recurse(s, control);
    control.checkpoint(s.local & 1023, s.bestSeen);
}

SearchResult solveTsp(const WeightedGraph& weights, int threads, double timeLimit, int frontierFactor,
                      long long knownLowerBound) {
    SearchResult result;
    result.algorithm = "branch-and-bound";
    result.threads = threads;
    double t0 = wallTime();
    int n = weights.n;
    if (n <= 2) {
        result.status = SearchStatus::Optimal;
        for (int v = 0; v < n; v++) result.path.push_back(v);
        result.cost = n == 2 ? weights.at(0, 1) + weights.at(1, 0) : 0;
        result.seconds = wallTime() - t0;
        return result;
    }
    TspBranchAndBound bb(weights);
    std::vector<int> tour;
    long long upper = bb.heuristicTour(tour, 50);
    SharedControl control;
    control.deadline.limit = timeLimit;
    control.deadline.start = t0;
    control.best = upper;
    control.solution = tour;
    double rootLb = bb.heldKarpRoot(upper, std::min(3000, 100 + 20 * n));
    result.lowerBound = std::max(rootLb, (double)knownLowerBound);
    bb.setProofTarget(std::max(ceilBound(rootLb), knownLowerBound));
    if (upper > bb.getProofTarget()) {
        std::vector<std::vector<int>> frontier = {{0}};
        if (threads > 1) frontier = bb.expand(frontier, (size_t)frontierFactor * threads, 4, control);
        long long count = (long long)frontier.size();
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (long long i = 0; i < count; i++) {
            if (control.isStopped()) continue;
            bb.solvePrefix(frontier[i], control);
        }
    }
    result.nodes = control.nodes;
    result.cost = control.best;
    result.path = control.solution;
    result.status = control.timedOut ? SearchStatus::Feasible : SearchStatus::Optimal;
    result.seconds = wallTime() - t0;
    return result;
}

WeightedGraph addZeroDummy(const WeightedGraph& weights) {
    WeightedGraph out(weights.n + 1);
    out.name = weights.name + "+dummy";
    out.weightType = weights.weightType;
    for (int i = 0; i < weights.n; i++)
        for (int j = 0; j < weights.n; j++) out.set(i, j, weights.at(i, j));
    return out;
}

WeightedGraph hamiltonicityToTsp(const Graph& graph) {
    WeightedGraph out(graph.V);
    out.name = graph.name + "-as-tsp";
    out.weightType = "HCP_0_1";
    for (int i = 0; i < graph.V; i++)
        for (int j = 0; j < graph.V; j++) out.set(i, j, i == j ? 0 : (graph.hasEdge(i, j) ? 0 : 1));
    return out;
}

Graph pathToCycleGraph(const Graph& graph) {
    Graph out(graph.V + 1, graph.directed);
    out.name = graph.name + "+universal";
    for (int u = 0; u < graph.V; u++)
        for (int v : graph.adjList[u])
            if (graph.directed || u < v) out.addEdge(u, v);
    for (int u = 0; u < graph.V; u++) {
        out.addEdge(graph.V, u);
        if (graph.directed) out.addEdge(u, graph.V);
    }
    out.finalize();
    return out;
}

std::vector<int> removeVertexFromCycle(const std::vector<int>& cycle, int vertex) {
    std::vector<int> path;
    auto it = std::find(cycle.begin(), cycle.end(), vertex);
    if (it == cycle.end()) return cycle;
    size_t k = (size_t)(it - cycle.begin());
    for (size_t i = 1; i < cycle.size(); i++) path.push_back(cycle[(k + i) % cycle.size()]);
    return path;
}
