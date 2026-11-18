#include "common.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <chrono>
#endif

std::string statusName(SearchStatus s) {
    switch (s) {
        case SearchStatus::Found: return "FOUND";
        case SearchStatus::None: return "NONE";
        case SearchStatus::Timeout: return "TIMEOUT";
        case SearchStatus::Optimal: return "OPTIMAL";
        case SearchStatus::Feasible: return "FEASIBLE";
        default: return "ERROR";
    }
}

std::string modeName(SearchMode m) {
    switch (m) {
        case SearchMode::Path: return "path";
        case SearchMode::Cycle: return "cycle";
        case SearchMode::Tsp: return "tsp";
        default: return "shp";
    }
}

double wallTime() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now().time_since_epoch()).count();
#endif
}

bool SharedControl::checkpoint(long long newNodes, long long& bestSeen) {
    bool expired = deadline.expired();
    bool stopped;
#pragma omp critical
    {
        nodes += newNodes;
        if (expired && !stop) {
            timedOut = true;
            stop = true;
        }
        bestSeen = best;
        stopped = stop;
    }
    return stopped;
}

bool SharedControl::isStopped() {
    bool stopped;
#pragma omp critical
    stopped = stop;
    return stopped;
}

long long SharedControl::currentBest() {
    long long value;
#pragma omp critical
    value = best;
    return value;
}

bool SharedControl::offerSolution(const std::vector<int>& path) {
    bool accepted = false;
#pragma omp critical
    {
        if (!found) {
            solution = path;
            found = true;
            stop = true;
            accepted = true;
        }
    }
    return accepted;
}

bool SharedControl::offerTour(long long cost, const std::vector<int>& tour, long long proofTarget,
                              long long& bestSeen) {
    bool stopped;
#pragma omp critical
    {
        if (cost < best) {
            best = cost;
            solution = tour;
            if (best <= proofTarget) stop = true;
        }
        bestSeen = best;
        stopped = stop;
    }
    return stopped;
}

void SharedControl::lowerBest(long long cost) {
#pragma omp critical
    {
        if (cost < best) best = cost;
    }
}
