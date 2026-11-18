#include "hamsearch.h"

#include <algorithm>


HamSearch::HamSearch(const Graph& graph, bool cycleMode)
    : g(graph), cycle(cycleMode), pruning(!graph.directed), start(0), connectivityEvery(graph.V <= 256 ? 1 : 16) {
    for (int v = 1; v < g.V; v++)
        if (g.degree(v) < g.degree(start)) start = v;
}

bool HamSearch::trivialAnswer(SearchResult& result) const {
    int n = g.V;
    if (n == 0) { result.status = SearchStatus::None; return true; }
    if (n == 1) { result.status = SearchStatus::Found; result.path = {0}; return true; }
    if (g.directed) return false;
    if (!g.isConnected()) { result.status = SearchStatus::None; return true; }
    int leaves = 0;
    for (int v = 0; v < n; v++) {
        if (g.degree(v) == 0) { result.status = SearchStatus::None; return true; }
        if (g.degree(v) == 1) leaves++;
    }
    if (cycle && (n < 3 || leaves > 0)) { result.status = SearchStatus::None; return true; }
    if (!cycle && leaves > 2) { result.status = SearchStatus::None; return true; }
    return false;
}

std::vector<std::vector<int>> HamSearch::roots() const {
    std::vector<std::vector<int>> result;
    if (g.V == 0) return result;
    if (cycle) {
        result.push_back({start});
        return result;
    }
    if (!g.directed)
        for (int v = 0; v < g.V; v++)
            if (g.degree(v) == 1) {
                result.push_back({v});
                return result;
            }
    std::vector<int> order(g.V);
    for (int v = 0; v < g.V; v++) order[v] = v;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return g.degree(a) < g.degree(b); });
    for (int v : order) result.push_back({v});
    return result;
}

void HamSearch::visit(State& s, int u) const {
    if (cycle && !g.directed) {
        if (s.path.size() == 1) {
            s.closersAbove = 0;
            for (int w : g.adjList[start])
                if (w > u && !s.visited[w]) s.closersAbove++;
        } else if (s.path.size() >= 2 && u > s.path[1] && g.hasEdge(u, start)) {
            s.closersAbove--;
        }
    }
    s.visited[u] = 1;
    s.path.push_back(u);
    if (pruning)
        for (int w : g.adjList[u]) s.avail[w]--;
}

void HamSearch::unvisit(State& s) const {
    int u = s.path.back();
    s.path.pop_back();
    s.visited[u] = 0;
    if (pruning)
        for (int w : g.adjList[u]) s.avail[w]++;
    if (cycle && !g.directed && s.path.size() >= 2 && u > s.path[1] && g.hasEdge(u, start)) s.closersAbove++;
}

void HamSearch::initState(State& s, const std::vector<int>& prefix) const {
    s.visited.assign(g.V, 0);
    s.avail.resize(g.V);
    for (int v = 0; v < g.V; v++) s.avail[v] = g.degree(v);
    s.stamp.assign(g.V, 0);
    s.stampCounter = 0;
    s.queue.clear();
    s.path.clear();
    s.path.reserve(g.V);
    s.closersAbove = 0;
    for (int v : prefix) visit(s, v);
}

bool HamSearch::complete(const State& s) const {
    if ((int)s.path.size() != g.V) return false;
    if (!cycle || g.V == 1) return true;
    int tail = s.path.back();
    if (!g.hasEdge(tail, s.path.front())) return false;
    if (!g.directed && g.V >= 3 && tail < s.path[1]) return false;
    return true;
}

bool HamSearch::connectedRemainder(State& s) const {
    int remaining = g.V - (int)s.path.size();
    int seed = -1;
    for (int v = 0; v < g.V; v++)
        if (!s.visited[v]) { seed = v; break; }
    if (seed < 0) return true;
    s.stampCounter++;
    s.queue.clear();
    s.queue.push_back(seed);
    s.stamp[seed] = s.stampCounter;
    size_t head = 0;
    while (head < s.queue.size()) {
        int u = s.queue[head++];
        for (int w : g.adjList[u])
            if (!s.visited[w] && s.stamp[w] != s.stampCounter) {
                s.stamp[w] = s.stampCounter;
                s.queue.push_back(w);
            }
    }
    return (int)s.queue.size() == remaining;
}

bool HamSearch::feasible(State& s, int previousTail) const {
    int remaining = g.V - (int)s.path.size();
    if (remaining == 0 || !pruning) return true;
    int tail = s.path.back();
    if (cycle && s.path.size() >= 2 && s.closersAbove == 0) return false;
    s.stampCounter++;
    int ones = 0;
    auto check = [&](int w) {
        if (s.visited[w] || s.stamp[w] == s.stampCounter) return true;
        s.stamp[w] = s.stampCounter;
        int potential = s.avail[w] + (g.hasEdge(w, tail) ? 1 : 0) + (cycle && g.hasEdge(w, start) ? 1 : 0);
        if (cycle) return potential >= 2;
        if (potential == 0) return false;
        if (potential == 1) {
            if (s.avail[w] == 0 && remaining > 1) return false;
            if (++ones >= 2) return false;
        }
        return true;
    };
    for (int w : g.adjList[tail])
        if (!check(w)) return false;
    if (previousTail >= 0)
        for (int w : g.adjList[previousTail])
            if (!check(w)) return false;
    if (remaining >= 2 && (s.path.size() % connectivityEvery) == 0 && !connectedRemainder(s)) return false;
    return true;
}

void HamSearch::candidates(const State& s, std::vector<int>& out) const {
    out.clear();
    int tail = s.path.back();
    for (int w : g.adjList[tail])
        if (!s.visited[w]) out.push_back(w);
    if (pruning)
        std::stable_sort(out.begin(), out.end(), [&](int a, int b) { return s.avail[a] < s.avail[b]; });
    else
        std::stable_sort(out.begin(), out.end(), [&](int a, int b) { return g.degree(a) < g.degree(b); });
}

std::vector<std::vector<int>> HamSearch::expand(std::vector<std::vector<int>> frontier, size_t target, int maxLevels,
                                                 SharedControl& control) const {
    State s;
    std::vector<int> cand;
    int levels = 0;
    while (!frontier.empty() && frontier.size() < target && levels < maxLevels && !control.stop) {
        std::vector<std::vector<int>> next;
        for (const auto& prefix : frontier) {
            initState(s, prefix);
            if ((int)s.path.size() == g.V) {
                if (complete(s)) { control.offerSolution(s.path); return {}; }
                continue;
            }
            candidates(s, cand);
            int previousTail = s.path.back();
            for (int u : cand) {
                visit(s, u);
                control.nodes++;
                if ((int)s.path.size() == g.V) {
                    if (complete(s)) { control.offerSolution(s.path); return {}; }
                } else if (feasible(s, previousTail)) {
                    next.push_back(s.path);
                }
                unvisit(s);
            }
        }
        frontier.swap(next);
        levels++;
    }
    return frontier;
}

void HamSearch::solvePrefix(const std::vector<int>& prefix, SharedControl& control) const {
    std::vector<long long> unused;
    solvePrefixBudget(prefix, control, 0, unused);
}

bool HamSearch::solvePrefixBudget(const std::vector<int>& prefix, SharedControl& control, long long budget,
                                  std::vector<long long>& leftover) const {
    State s;
    initState(s, prefix);
    if ((int)s.path.size() == g.V) {
        if (complete(s)) control.offerSolution(s.path);
        return true;
    }
    std::vector<std::vector<int>> cand(g.V + 1);
    std::vector<size_t> index(g.V + 1, 0);
    int base = (int)s.path.size();
    int depth = base;
    candidates(s, cand[depth]);
    long long local = 0;
    long long bestSeen = 0;
    bool stopped = control.isStopped();
    while (!stopped) {
        if (index[depth] >= cand[depth].size()) {
            if (depth == base) break;
            unvisit(s);
            depth--;
            continue;
        }
        if (budget > 0 && local >= budget) {
            long long tasks = 0;
            for (int d = base; d <= depth; d++) tasks += (long long)(cand[d].size() - index[d]);
            leftover.push_back((long long)s.path.size());
            for (int v : s.path) leftover.push_back(v);
            leftover.push_back(tasks);
            for (int d = base; d <= depth; d++)
                for (size_t j = index[d]; j < cand[d].size(); j++) {
                    leftover.push_back(d);
                    leftover.push_back(cand[d][j]);
                }
            control.checkpoint(local & 1023, bestSeen);
            return false;
        }
        int u = cand[depth][index[depth]++];
        int previousTail = s.path.back();
        visit(s, u);
        local++;
        if ((local & 1023) == 0 && control.checkpoint(1024, bestSeen)) break;
        if ((int)s.path.size() == g.V) {
            if (complete(s)) {
                control.offerSolution(s.path);
                break;
            }
            unvisit(s);
            continue;
        }
        if (!feasible(s, previousTail)) {
            unvisit(s);
            continue;
        }
        depth++;
        candidates(s, cand[depth]);
        index[depth] = 0;
    }
    control.checkpoint(local & 1023, bestSeen);
    return true;
}

SearchResult solveHamiltonian(const Graph& graph, bool cycle, int threads, double timeLimit, int frontierFactor) {
    SearchResult result;
    result.algorithm = "pruned-dfs";
    result.threads = threads;
    double t0 = wallTime();
    HamSearch search(graph, cycle);
    if (search.trivialAnswer(result)) {
        result.seconds = wallTime() - t0;
        return result;
    }
    SharedControl control;
    control.deadline.limit = timeLimit;
    control.deadline.start = t0;
    std::vector<std::vector<int>> frontier = search.roots();
    if (threads > 1) frontier = search.expand(frontier, (size_t)frontierFactor * threads, 64, control);
    if (!control.found) {
        long long count = (long long)frontier.size();
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (long long i = 0; i < count; i++) {
            if (control.isStopped()) continue;
            search.solvePrefix(frontier[i], control);
        }
    }
    result.nodes = control.nodes;
    if (control.found) {
        result.status = SearchStatus::Found;
        result.path = control.solution;
    } else if (control.timedOut) {
        result.status = SearchStatus::Timeout;
    } else {
        result.status = SearchStatus::None;
    }
    result.seconds = wallTime() - t0;
    return result;
}
