#ifdef USE_MPI

#include "mpi_driver.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <vector>


#include "hamsearch.h"
#include "tsp.h"

namespace {

enum Tag { TAG_REQUEST = 1, TAG_DATA = 2, TAG_WORK = 3, TAG_STOP = 4 };
enum Report { REPORT_NOTHING = 0, REPORT_PATH = 1, REPORT_TOUR = 2 };

const int HEADER = 5;
const long long TASK_BUDGET = 1LL << 19;

struct Task {
    int pathId;
    int length;
    int vertex;
};

struct ManagerOutcome {
    bool found = false;
    bool timedOut = false;
    long long best = 0;
    std::vector<int> solution;
};

void addLeftovers(const std::vector<long long>& extra, std::vector<std::vector<int>>& pool, std::deque<Task>& queue) {
    size_t k = 0;
    while (k < extra.size()) {
        long long length = extra[k++];
        std::vector<int> path(extra.begin() + k, extra.begin() + k + length);
        k += length;
        long long count = extra[k++];
        int id = (int)pool.size();
        pool.push_back(path);
        std::vector<Task> tasks;
        for (long long t = 0; t < count; t++, k += 2) tasks.push_back({id, (int)extra[k], (int)extra[k + 1]});
        std::stable_sort(tasks.begin(), tasks.end(), [](const Task& a, const Task& b) { return a.length > b.length; });
        for (auto it = tasks.rbegin(); it != tasks.rend(); ++it) queue.push_front(*it);
    }
}

ManagerOutcome runManager(const std::vector<std::vector<int>>& frontier, SharedControl& control, bool optimization,
                          long long proofTarget, int size, int vertices) {
    ManagerOutcome outcome;
    outcome.found = control.found;
    outcome.best = control.best;
    outcome.solution = control.solution;
    bool stop = outcome.found || (optimization && outcome.best <= proofTarget);
    std::vector<std::vector<int>> pool;
    std::deque<Task> queue;
    for (const auto& prefix : frontier) {
        queue.push_back({(int)pool.size(), (int)prefix.size(), -1});
        pool.push_back(prefix);
    }
    int active = size - 1;
    int capacity = HEADER + vertices + 1;
    std::vector<long long> header(capacity);
    while (active > 0) {
        MPI_Status status;
        MPI_Recv(header.data(), capacity, MPI_LONG_LONG, MPI_ANY_SOURCE, TAG_REQUEST, MPI_COMM_WORLD, &status);
        int worker = status.MPI_SOURCE;
        long long kind = header[1], cost = header[2], length = header[3], extraLength = header[4];
        if (header[0]) outcome.timedOut = true;
        if (extraLength > 0) {
            std::vector<long long> extra(extraLength);
            MPI_Recv(extra.data(), (int)extraLength, MPI_LONG_LONG, worker, TAG_DATA, MPI_COMM_WORLD, &status);
            if (!stop) addLeftovers(extra, pool, queue);
        }
        if (kind == REPORT_PATH && !outcome.found) {
            outcome.found = true;
            outcome.solution.assign(header.begin() + HEADER, header.begin() + HEADER + length);
            stop = true;
        } else if (kind == REPORT_TOUR && cost < outcome.best) {
            outcome.best = cost;
            outcome.solution.assign(header.begin() + HEADER, header.begin() + HEADER + length);
            if (cost <= proofTarget) stop = true;
        }
        if (control.deadline.expired()) {
            outcome.timedOut = true;
            stop = true;
        }
        if (stop || queue.empty()) {
            long long done = 0;
            MPI_Send(&done, 1, MPI_LONG_LONG, worker, TAG_STOP, MPI_COMM_WORLD);
            active--;
        } else {
            Task task = queue.front();
            queue.pop_front();
            const std::vector<int>& base = pool[task.pathId];
            std::vector<long long> work = {outcome.best, (long long)task.length + (task.vertex >= 0 ? 1 : 0)};
            for (int i = 0; i < task.length; i++) work.push_back(base[i]);
            if (task.vertex >= 0) work.push_back(task.vertex);
            MPI_Send(work.data(), (int)work.size(), MPI_LONG_LONG, worker, TAG_WORK, MPI_COMM_WORLD);
        }
    }
    return outcome;
}

void sendReport(int timedOut, long long kind, long long cost, const std::vector<int>& path,
                const std::vector<long long>& extra) {
    std::vector<long long> header = {timedOut, kind, cost, (long long)path.size(), (long long)extra.size()};
    for (int v : path) header.push_back(v);
    MPI_Send(header.data(), (int)header.size(), MPI_LONG_LONG, 0, TAG_REQUEST, MPI_COMM_WORLD);
    if (!extra.empty())
        MPI_Send(extra.data(), (int)extra.size(), MPI_LONG_LONG, 0, TAG_DATA, MPI_COMM_WORLD);
}

void runHamWorker(const HamSearch& search, SharedControl& control, int threads, int frontierFactor, int vertices) {
    long long kind = REPORT_NOTHING;
    std::vector<int> path;
    std::vector<long long> extra;
    std::vector<long long> buffer(vertices + 3);
    while (true) {
        sendReport(control.timedOut ? 1 : 0, kind, 0, path, extra);
        MPI_Status status;
        MPI_Recv(buffer.data(), (int)buffer.size(), MPI_LONG_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) break;
        kind = REPORT_NOTHING;
        path.clear();
        extra.clear();
        std::vector<int> prefix(buffer.begin() + 2, buffer.begin() + 2 + buffer[1]);
        std::vector<std::vector<int>> local = {prefix};
        if (threads > 1) local = search.expand(local, (size_t)frontierFactor * threads, 8, control);
        long long count = (long long)local.size();
        std::vector<std::vector<long long>> leftovers(count);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (long long i = 0; i < count; i++) {
            if (control.isStopped()) continue;
            search.solvePrefixBudget(local[i], control, TASK_BUDGET, leftovers[i]);
        }
        if (control.found) {
            kind = REPORT_PATH;
            path = control.solution;
        } else if (!control.stop) {
            for (auto& part : leftovers) extra.insert(extra.end(), part.begin(), part.end());
        }
    }
}

void runTspWorker(const TspBranchAndBound& bb, SharedControl& control, int threads, int frontierFactor,
                  int vertices) {
    long long reported = control.best;
    long long kind = REPORT_NOTHING, cost = 0;
    std::vector<int> path;
    std::vector<long long> none;
    std::vector<long long> buffer(vertices + 3);
    while (true) {
        sendReport(control.timedOut ? 1 : 0, kind, cost, path, none);
        MPI_Status status;
        MPI_Recv(buffer.data(), (int)buffer.size(), MPI_LONG_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == TAG_STOP) break;
        kind = REPORT_NOTHING;
        path.clear();
        control.lowerBest(buffer[0]);
        std::vector<int> prefix(buffer.begin() + 2, buffer.begin() + 2 + buffer[1]);
        std::vector<std::vector<int>> local = {prefix};
        if (threads > 1) local = bb.expand(local, (size_t)frontierFactor * threads, 8, control);
        long long count = (long long)local.size();
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
        for (long long i = 0; i < count; i++) {
            if (control.isStopped()) continue;
            bb.solvePrefix(local[i], control);
        }
        if (control.best < reported) {
            reported = control.best;
            kind = REPORT_TOUR;
            cost = reported;
            path = control.solution;
        }
    }
}

long long reduceNodes(long long local) {
    long long total = 0;
    MPI_Reduce(&local, &total, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    return total;
}

}

SearchResult mpiSolveHamiltonian(const Graph& graph, bool cycle, int threads, double timeLimit, int frontierFactor) {
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    SearchResult result;
    result.algorithm = "pruned-dfs-mpi";
    result.threads = threads;
    result.ranks = size;
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    HamSearch search(graph, cycle);
    if (search.trivialAnswer(result)) {
        result.seconds = MPI_Wtime() - t0;
        return result;
    }
    if (size == 1) {
        result = solveHamiltonian(graph, cycle, threads, timeLimit, frontierFactor);
        result.algorithm = "pruned-dfs-mpi";
        result.ranks = 1;
        return result;
    }
    SharedControl control;
    control.deadline.limit = timeLimit;
    control.deadline.start = wallTime();
    if (rank == 0) {
        std::vector<std::vector<int>> frontier =
            search.expand(search.roots(), (size_t)frontierFactor * (size - 1), 64, control);
        ManagerOutcome outcome = runManager(frontier, control, false, 0, size, graph.V);
        result.status = outcome.found ? SearchStatus::Found
                                      : (outcome.timedOut ? SearchStatus::Timeout : SearchStatus::None);
        result.path = outcome.solution;
    } else {
        runHamWorker(search, control, threads, frontierFactor, graph.V);
    }
    result.nodes = reduceNodes(control.nodes);
    MPI_Barrier(MPI_COMM_WORLD);
    result.seconds = MPI_Wtime() - t0;
    return result;
}

SearchResult mpiSolveTsp(const WeightedGraph& weights, int threads, double timeLimit, int frontierFactor,
                         long long knownLowerBound) {
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size == 1 || weights.n <= 3) {
        SearchResult r = solveTsp(weights, threads, timeLimit, frontierFactor, knownLowerBound);
        r.algorithm = "branch-and-bound-mpi";
        r.ranks = size;
        return r;
    }
    SearchResult result;
    result.algorithm = "branch-and-bound-mpi";
    result.threads = threads;
    result.ranks = size;
    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    TspBranchAndBound bb(weights);
    std::vector<int> tour;
    long long upper = bb.heuristicTour(tour, 50);
    double rootLb = bb.heldKarpRoot(upper, std::min(3000, 100 + 20 * weights.n));
    long long target = std::max((long long)std::ceil(rootLb - 1e-6), knownLowerBound);
    bb.setProofTarget(target);
    SharedControl control;
    control.deadline.limit = timeLimit;
    control.deadline.start = wallTime();
    control.best = upper;
    control.solution = tour;
    result.lowerBound = std::max(rootLb, (double)knownLowerBound);
    if (rank == 0) {
        std::vector<std::vector<int>> frontier;
        if (upper > target) frontier = bb.expand({{0}}, (size_t)frontierFactor * (size - 1) * 4, 4, control);
        ManagerOutcome outcome = runManager(frontier, control, true, target, size, weights.n);
        result.cost = outcome.best;
        result.path = outcome.solution;
        result.status = outcome.timedOut ? SearchStatus::Feasible : SearchStatus::Optimal;
    } else {
        runTspWorker(bb, control, threads, frontierFactor, weights.n);
    }
    result.nodes = reduceNodes(control.nodes);
    MPI_Barrier(MPI_COMM_WORLD);
    result.seconds = MPI_Wtime() - t0;
    return result;
}

#endif
