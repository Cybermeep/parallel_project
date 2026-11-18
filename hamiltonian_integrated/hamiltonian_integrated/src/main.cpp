#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "common.h"
#include "graph.h"
#include "hamsearch.h"
#include "io.h"
#include "mpi_driver.h"
#include "sam_dfs.h"
#include "tsp.h"
#include "verify.h"
#include "visualize.h"

namespace {

struct Options {
    std::string file;
    std::string mode;
    std::string algo = "auto";
    std::string expect;
    std::string tourFile;
    std::string csvFile;
    std::string label;
    int threads = 1;
    double timeLimit = 0.0;
    int frontier = 16;
    int start = 1;
    bool allStarts = false;
    bool directed = false;
    int vertices = 0;
    int graphIndex = 1;
    bool batch = false;
    int first = 1;
    int count = 0;
    bool printMatrix = false;
    bool createCode = false;
    bool visualizeDemo = false;
    bool printPath = false;
    bool quiet = false;
};

int rankId = 0;
int rankCount = 1;

void usage() {
    std::cout
        << "Usage: hampath <graph-file> [options]\n"
        << "  --mode path|cycle|tsp|shp       problem to solve (default: path for graphs, tsp for weighted files)\n"
        << "  --algo auto|pruned|sam-iterative|sam-recursive|hcp-as-tsp|bnb\n"
        << "  --threads N                     OpenMP threads per rank\n"
        << "  --time-limit S                  wall-clock limit in seconds (0 = none)\n"
        << "  --frontier F                    partial paths generated per thread for load balancing\n"
        << "  --start V                       start vertex for sam-* solvers (1-indexed)\n"
        << "  --all-starts                    sam-* solvers try every start vertex in path mode\n"
        << "  --directed                      treat edge lists as directed\n"
        << "  --vertices N                    force vertex count for plain edge-list files\n"
        << "  --graph K                       pick graph K from a multi-graph file (.txt/.asc/.g6)\n"
        << "  --batch [--first K --count C]   solve every graph in a multi-graph file\n"
        << "  --expect hamiltonian|non-hamiltonian|<optimal-cost>\n"
        << "  --tour FILE                     known optimal TSPLIB .tour file for comparison\n"
        << "  --print-matrix                  print the adjacency matrix (Alysa's print())\n"
        << "  --create-code                   print the Graphviz shell command (Alysa's createCode())\n"
        << "  --visualize-demo                run Alysa's 4-vertex visualizeGraph demo (no file needed)\n"
        << "  --print-path                    always print the full solution\n"
        << "  --csv FILE                      append the RESULT row to a CSV file\n"
        << "  --label NAME                    instance label used in RESULT rows\n"
        << "  --quiet                         only print the RESULT row\n";
}

bool parseArgs(int argc, char** argv, Options& o) {
#ifdef _OPENMP
    o.threads = omp_get_num_procs();
#endif
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + a);
            return argv[++i];
        };
        if (a == "-h" || a == "--help") return false;
        else if (a == "--mode") o.mode = value();
        else if (a == "--algo") o.algo = value();
        else if (a == "--threads") o.threads = std::max(1, std::stoi(value()));
        else if (a == "--time-limit") o.timeLimit = std::stod(value());
        else if (a == "--frontier") o.frontier = std::max(1, std::stoi(value()));
        else if (a == "--start") o.start = std::stoi(value());
        else if (a == "--all-starts") o.allStarts = true;
        else if (a == "--directed") o.directed = true;
        else if (a == "--vertices") o.vertices = std::stoi(value());
        else if (a == "--graph") o.graphIndex = std::stoi(value());
        else if (a == "--batch") o.batch = true;
        else if (a == "--first") o.first = std::stoi(value());
        else if (a == "--count") o.count = std::stoi(value());
        else if (a == "--expect") o.expect = value();
        else if (a == "--tour") o.tourFile = value();
        else if (a == "--csv") o.csvFile = value();
        else if (a == "--label") o.label = value();
        else if (a == "--print-matrix") o.printMatrix = true;
        else if (a == "--create-code") o.createCode = true;
        else if (a == "--visualize-demo") o.visualizeDemo = true;
        else if (a == "--print-path") o.printPath = true;
        else if (a == "--quiet") o.quiet = true;
        else if (!a.empty() && a[0] == '-') throw std::runtime_error("unknown option " + a);
        else o.file = a;
    }
    return !o.file.empty() || o.visualizeDemo;
}

std::string formatPath(const std::vector<int>& path, bool cycle) {
    std::ostringstream ss;
    for (size_t i = 0; i < path.size(); i++) ss << (i ? " -> " : "") << path[i] + 1;
    if (cycle && !path.empty()) ss << " -> " << path.front() + 1;
    return ss.str();
}

std::string verdictFor(const std::string& expect, SearchMode mode, const SearchResult& r) {
    if (expect.empty()) return "UNKNOWN";
    if (expect == "hamiltonian") {
        if (r.status == SearchStatus::Found) return "MATCH";
        if (r.status == SearchStatus::None) return "MISMATCH";
        return "UNSOLVED";
    }
    if (expect == "non-hamiltonian") {
        if (mode == SearchMode::Path) return r.status == SearchStatus::None ? "MATCH" : "NA";
        if (r.status == SearchStatus::None) return "MATCH";
        if (r.status == SearchStatus::Found) return "MISMATCH";
        return "UNSOLVED";
    }
    long long known = std::stoll(expect);
    if (r.cost < 0) return "UNSOLVED";
    if (r.cost < known) return "MISMATCH";
    if (r.status == SearchStatus::Optimal) return r.cost == known ? "MATCH" : "MISMATCH";
    std::ostringstream ss;
    ss << "GAP=" << std::fixed << std::setprecision(2) << 100.0 * (double)(r.cost - known) / std::max(1LL, known) << "%";
    return ss.str();
}

void emitResult(const Options& o, const std::string& instance, SearchMode mode, const SearchResult& r, int n,
                long long m, const std::string& status, const std::string& verdict, const std::string& verified) {
    std::ostringstream row;
    row << instance << "," << modeName(mode) << "," << r.algorithm << "," << r.ranks << "," << r.threads << "," << n
        << "," << m << "," << status << "," << r.cost << "," << (o.expect.empty() ? "-" : o.expect) << "," << verdict
        << "," << std::fixed << std::setprecision(6) << r.seconds << "," << r.nodes << "," << std::setprecision(2)
        << r.lowerBound << "," << verified;
    std::cout << "RESULT," << row.str() << std::endl;
    if (!o.csvFile.empty()) {
        bool fresh = !std::ifstream(o.csvFile).good();
        std::ofstream csv(o.csvFile, std::ios::app);
        if (fresh)
            csv << "instance,mode,algorithm,ranks,threads,vertices,edges,status,cost,expected,verdict,seconds,nodes,"
                   "lower_bound,verified\n";
        csv << row.str() << "\n";
    }
}

SearchMode parseMode(const Options& o, bool weighted) {
    std::string m = o.mode.empty() ? (weighted ? "tsp" : "path") : o.mode;
    if (m == "path") return SearchMode::Path;
    if (m == "cycle") return SearchMode::Cycle;
    if (m == "tsp") return SearchMode::Tsp;
    if (m == "shp") return SearchMode::ShortestPath;
    throw std::runtime_error("unknown mode " + m);
}

SearchResult solveHamiltonicityViaTsp(const Graph& g, bool cycle, const Options& o, bool distributed) {
    Graph augmented = cycle ? Graph() : pathToCycleGraph(g);
    const Graph& target = cycle ? g : augmented;
    WeightedGraph wg = hamiltonicityToTsp(target);
    SearchResult r;
#ifdef USE_MPI
    r = distributed ? mpiSolveTsp(wg, o.threads, o.timeLimit, o.frontier, 0) : solveTsp(wg, o.threads, o.timeLimit, o.frontier, 0);
#else
    (void)distributed;
    r = solveTsp(wg, o.threads, o.timeLimit, o.frontier, 0);
#endif
    r.algorithm = "hcp-as-tsp";
    if (r.cost == 0) {
        r.status = SearchStatus::Found;
        if (!cycle) r.path = removeVertexFromCycle(r.path, g.V);
    } else {
        r.status = r.status == SearchStatus::Optimal ? SearchStatus::None : SearchStatus::Timeout;
        r.path.clear();
    }
    return r;
}

SearchResult solveGraphInstance(const Graph& g, SearchMode mode, const Options& o, bool distributed) {
    bool cycle = mode == SearchMode::Cycle;
    std::string algo = o.algo == "auto" ? "pruned" : o.algo;
    Deadline deadline;
    deadline.limit = o.timeLimit;
    deadline.start = wallTime();
    SearchResult r;
    if (algo == "pruned") {
#ifdef USE_MPI
        if (distributed) return mpiSolveHamiltonian(g, cycle, o.threads, o.timeLimit, o.frontier);
#endif
        return solveHamiltonian(g, cycle, o.threads, o.timeLimit, o.frontier);
    }
    if (algo == "hcp-as-tsp") return solveHamiltonicityViaTsp(g, cycle, o, distributed);
    if (algo == "sam-iterative" || algo == "sam-recursive") {
        if (distributed && rankId != 0) {
            r.algorithm = algo;
            return r;
        }
        int startVertex = std::min(std::max(o.start, 1), std::max(1, g.V)) - 1;
        r = algo == "sam-iterative" ? dfsIterative(g, cycle, startVertex, o.allStarts, deadline)
                                    : dfsRecursive(g, cycle, startVertex, o.allStarts, deadline);
        return r;
    }
    throw std::runtime_error("algorithm " + algo + " does not apply to mode " + modeName(mode));
}

int runBatch(const Options& o, std::vector<Instance>& instances) {
    SearchMode mode = parseMode(o, false);
    int total = (int)instances.size();
    int firstIndex = std::max(1, o.first) - 1;
    int lastIndex = o.count > 0 ? std::min(total, firstIndex + o.count) : total;
    Options single = o;
    single.threads = 1;
    long long foundCount = 0, noneCount = 0, timeoutCount = 0, invalidCount = 0, nodeCount = 0;
    std::vector<int> mismatches;
    double t0 = wallTime();
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    t0 = MPI_Wtime();
#endif
#pragma omp parallel for schedule(dynamic, 16) num_threads(o.threads) \
    reduction(+ : foundCount, noneCount, timeoutCount, invalidCount, nodeCount)
    for (int i = firstIndex; i < lastIndex; i++) {
        if ((i - firstIndex) % rankCount != rankId) continue;
        const Graph& g = instances[i].graph;
        SearchResult r = solveGraphInstance(g, mode, single, false);
        std::string reason;
        if (r.status == SearchStatus::Found && !verifyHamiltonian(g, r.path, mode == SearchMode::Cycle, reason))
            invalidCount++;
        else if (r.status == SearchStatus::Found)
            foundCount++;
        else if (r.status == SearchStatus::None)
            noneCount++;
        else
            timeoutCount++;
        nodeCount += r.nodes;
        std::string v = verdictFor(o.expect, mode, r);
        if (v == "MISMATCH") {
#pragma omp critical
            mismatches.push_back(i + 1);
        }
    }
    long long totals[5] = {foundCount, noneCount, timeoutCount, invalidCount, nodeCount};
#ifdef USE_MPI
    long long counts[5] = {foundCount, noneCount, timeoutCount, invalidCount, nodeCount};
    MPI_Reduce(counts, totals, 5, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    int localMismatch = (int)mismatches.size(), totalMismatch = 0;
    MPI_Reduce(&localMismatch, &totalMismatch, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    int totalMismatch = (int)mismatches.size();
#endif
#ifdef USE_MPI
    double seconds = MPI_Wtime() - t0;
#else
    double seconds = wallTime() - t0;
#endif
    if (rankId != 0) return 0;
    long long solved = lastIndex - firstIndex;
    std::string verdict = "UNKNOWN";
    if (o.expect == "hamiltonian") verdict = totals[0] == solved ? "MATCH" : (totalMismatch ? "MISMATCH" : "UNSOLVED");
    else if (o.expect == "non-hamiltonian" && mode == SearchMode::Cycle)
        verdict = totals[1] == solved ? "MATCH" : (totalMismatch ? "MISMATCH" : "UNSOLVED");
    if (totals[3] > 0) verdict = "INVALID";
    if (!o.quiet) {
        std::cout << "Batch " << o.file << " graphs " << firstIndex + 1 << ".." << lastIndex << " (" << solved
                  << "), mode " << modeName(mode) << std::endl;
        std::cout << "  found=" << totals[0] << " none=" << totals[1] << " timeout=" << totals[2]
                  << " invalid=" << totals[3] << " mismatches=" << totalMismatch << std::endl;
        if (!mismatches.empty()) {
            std::sort(mismatches.begin(), mismatches.end());
            std::cout << "  first mismatching graphs on rank 0:";
            for (size_t k = 0; k < std::min<size_t>(10, mismatches.size()); k++) std::cout << " " << mismatches[k];
            std::cout << std::endl;
        }
        std::cout << "  time " << seconds << " s" << std::endl;
    }
    SearchResult summary;
    summary.algorithm = "batch-" + (o.algo == "auto" ? std::string("pruned") : o.algo);
    summary.threads = o.threads;
    summary.ranks = rankCount;
    summary.seconds = seconds;
    summary.nodes = totals[4];
    std::ostringstream status;
    status << "BATCH(found=" << totals[0] << ";none=" << totals[1] << ";timeout=" << totals[2] << ")";
    std::string label = o.label.empty() ? baseName(o.file) + "[" + std::to_string(firstIndex + 1) + "-" + std::to_string(lastIndex) + "]" : o.label;
    int n = instances.empty() ? 0 : instances[firstIndex].graph.V;
    emitResult(o, label, mode, summary, n, (long long)solved, status.str(), verdict, totals[3] ? "no" : "yes");
    return verdict == "MISMATCH" || verdict == "INVALID" ? 2 : 0;
}

int runSingle(const Options& o, Instance& inst) {
    SearchMode mode = parseMode(o, inst.weighted);
    Graph& g = inst.graph;
    bool distributed = rankCount > 1;
    if (rankId == 0 && !o.quiet) {
        std::cout << "Instance " << inst.name << ": " << g.V << " vertices, " << g.E << " edges"
                  << (inst.weighted ? " (weighted " + inst.weights.weightType + ")" : "") << ", mode " << modeName(mode)
                  << ", " << rankCount << " rank(s) x " << o.threads << " thread(s)" << std::endl;
        if (o.printMatrix) printAdjacencyMatrix(g);
        if (o.createCode) createCode(g);
    }
    SearchResult r;
    std::string verified = "-";
    std::string reason;
    if (mode == SearchMode::Path || mode == SearchMode::Cycle) {
        //run the algorithm to get the path(s)
        r = solveGraphInstance(g, mode, o, distributed);
    } else {
        if (!inst.weighted) {
            if (rankId == 0)
                std::cerr << "Modes tsp/shp need a weighted TSPLIB file. For Hamiltonicity through TSP use "
                             "--mode cycle --algo hcp-as-tsp."
                          << std::endl;
            return 1;
        }
        WeightedGraph problem = mode == SearchMode::Tsp ? inst.weights : addZeroDummy(inst.weights);
#ifdef USE_MPI
        r = distributed ? mpiSolveTsp(problem, o.threads, o.timeLimit, o.frontier, 0)
                        : solveTsp(problem, o.threads, o.timeLimit, o.frontier, 0);
#else
        r = solveTsp(problem, o.threads, o.timeLimit, o.frontier, 0);
#endif
        if (mode == SearchMode::ShortestPath && !r.path.empty()) r.path = removeVertexFromCycle(r.path, inst.weights.n);
    }
    r.ranks = rankCount;
    if (o.threads > 0) r.threads = o.threads;
    if (rankId != 0) return 0;

    if (mode == SearchMode::Path || mode == SearchMode::Cycle) {
        if (r.status == SearchStatus::Found) {
            verified = verifyHamiltonian(g, r.path, mode == SearchMode::Cycle, reason) ? "yes" : "no:" + reason;
            if (inst.weighted) r.cost = pathCost(inst.weights, r.path, mode == SearchMode::Cycle);
        }
    } else if (!r.path.empty()) {
        verified = verifyTour(inst.weights, r.path, mode == SearchMode::Tsp, r.cost, reason) ? "yes" : "no:" + reason;
    }

    Options effective = o;
    if (effective.expect.empty() && !o.tourFile.empty() && inst.weighted) {
        std::vector<int> known = readTourFile(o.tourFile);
        std::string tourReason;
        if (verifyTour(inst.weights, known, true, pathCost(inst.weights, known, true), tourReason))
            effective.expect = std::to_string(pathCost(inst.weights, known, true));
    }
    std::string verdict = verdictFor(effective.expect, mode, r);
    if (verified.rfind("no", 0) == 0) verdict = "INVALID";

    if (!o.quiet) {
        std::cout << "Algorithm " << r.algorithm << " finished in " << std::fixed << std::setprecision(6) << r.seconds
                  << " s (" << r.nodes << " search nodes)" << std::endl;
        std::cout << "Status: " << statusName(r.status);
        if (mode == SearchMode::Tsp || mode == SearchMode::ShortestPath || (inst.weighted && r.cost >= 0))
            std::cout << ", cost " << r.cost;
        if (r.lowerBound > 0) std::cout << ", root lower bound " << std::setprecision(2) << r.lowerBound;
        std::cout << std::endl;
        //print path(s)
        if (!r.path.empty() && (o.printPath || r.path.size() <= 60))
            std::cout << formatPath(r.path, mode == SearchMode::Cycle || mode == SearchMode::Tsp) << std::endl;
        if (!effective.expect.empty()) std::cout << "Expected " << effective.expect << ": " << verdict << std::endl;
        if (verified != "-") std::cout << "Independent verification: " << verified << std::endl;
    }

    std::string status = statusName(r.status);
    std::string label = o.label.empty() ? inst.name : o.label;
    emitResult(effective, label, mode, r, g.V, g.E, status, verdict, verified);
    return verdict == "MISMATCH" || verdict == "INVALID" ? 2 : 0;
}

}

int main(int argc, char* argv[]) {
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankId);
    MPI_Comm_size(MPI_COMM_WORLD, &rankCount);
#endif
    int code = 0;
    try {
        Options o;
        if (!parseArgs(argc, argv, o)) {
            if (rankId == 0) usage();
            code = 1;
        } else if (o.visualizeDemo && o.file.empty()) {
            if (rankId == 0) code = visualizeDemo();
        } else {
            double readStart = wallTime();
            //process file and convert to graph
            std::vector<Instance> instances = loadInstances(o.file, o.directed, o.vertices);
            if (instances.empty()) throw std::runtime_error("no graphs found in " + o.file);
            if (rankId == 0 && !o.quiet)
                std::cout << "Read " << instances.size() << " instance(s) from " << o.file << " in " << std::fixed
                          << std::setprecision(4) << wallTime() - readStart << " s" << std::endl;
            if (o.batch) {
                code = runBatch(o, instances);
            } else {
                int index = std::min(std::max(o.graphIndex, 1), (int)instances.size()) - 1;
                code = runSingle(o, instances[index]);
            }
        }
    } catch (const std::exception& e) {
        if (rankId == 0) std::cerr << "Error: " << e.what() << std::endl;
        code = 1;
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return code;
}
