#ifndef HAM_COMMON_H
#define HAM_COMMON_H

#include <limits>
#include <string>
#include <vector>

enum class SearchMode { Path, Cycle, Tsp, ShortestPath };

enum class SearchStatus { Found, None, Timeout, Optimal, Feasible, Error };

std::string statusName(SearchStatus s);
std::string modeName(SearchMode m);

double wallTime();

struct Deadline {
    double limit = 0.0;
    double start = 0.0;
    bool active() const { return limit > 0.0; }
    bool expired() const { return active() && wallTime() - start > limit; }
};

struct SharedControl {
    bool stop = false;
    bool found = false;
    bool timedOut = false;
    long long nodes = 0;
    long long best = std::numeric_limits<long long>::max();
    std::vector<int> solution;
    Deadline deadline;

    bool checkpoint(long long newNodes, long long& bestSeen);
    bool isStopped();
    long long currentBest();
    bool offerSolution(const std::vector<int>& path);
    bool offerTour(long long cost, const std::vector<int>& tour, long long proofTarget, long long& bestSeen);
    void lowerBest(long long cost);
};

struct SearchResult {
    SearchStatus status = SearchStatus::None;
    std::vector<int> path;
    long long cost = -1;
    long long nodes = 0;
    double seconds = 0.0;
    double lowerBound = 0.0;
    std::string algorithm;
    int threads = 1;
    int ranks = 1;
};

#endif
