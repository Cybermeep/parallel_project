#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include <chrono>

namespace {

struct Options {
    std::string bin = "./hampath";
    std::string mpiBin = "./hampath_mpi";
    std::string launcher = "mpirun -np {ranks}";
    std::vector<std::string> manifests;
    std::vector<std::string> suites = {"smoke"};
    std::string filter;
    std::vector<std::string> algos = {"pruned"};
    std::vector<int> threads = {1, 2, 4};
    std::vector<int> ranks = {1};
    int repeat = 1;
    double timeLimit = 60.0;
    int frontier = 16;
    std::string out = "results/bench.csv";
    std::string summary;
    bool dryRun = false;
};

struct Instance {
    std::string suite, name, file, mode, expect, args, source;
};

struct Row {
    Instance inst;
    std::string algo;
    int ranks = 1, threads = 1;
    std::map<std::string, std::string> result;
    double seconds = 0.0;
    double wall = 0.0;
    double speedup = -1.0;
    double efficiency = -1.0;
    std::string command;
};

const std::vector<std::string> RESULT_FIELDS = {"instance", "mode", "algorithm", "ranks", "threads", "vertices",
                                                "edges", "status", "cost", "expected", "verdict", "seconds",
                                                "nodes", "lower_bound", "verified"};

std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> parts;
    std::string cur;
    std::istringstream ss(s);
    while (std::getline(ss, cur, sep)) parts.push_back(cur);
    if (!s.empty() && s.back() == sep) parts.push_back("");
    return parts;
}

std::vector<int> splitInts(const std::string& s) {
    std::vector<int> v;
    for (auto& p : split(s, ','))
        if (!p.empty()) v.push_back(std::stoi(p));
    return v;
}

std::vector<std::string> splitWords(const std::string& s) {
    std::vector<std::string> v;
    for (auto& p : split(s, ','))
        if (!p.empty()) v.push_back(p);
    return v;
}

std::string quote(const std::string& s) {
    std::string q = "'";
    for (char c : s) {
        if (c == '\'') q += "'\\''";
        else q += c;
    }
    return q + "'";
}

std::string trimLine(std::string s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back();
    return s;
}

void usage() {
    std::cout << "Usage: run_bench [options]\n"
              << "  --bin PATH              OpenMP binary (default ./hampath)\n"
              << "  --mpi-bin PATH          MPI binary used when ranks > 1 (default ./hampath_mpi)\n"
              << "  --launcher CMD          MPI launcher, {ranks} is replaced (default \"mpirun -np {ranks}\")\n"
              << "  --manifest FILE         instance list, repeatable (default bench/instances.csv)\n"
              << "  --suite a,b|all         suites to run (default smoke)\n"
              << "  --filter TEXT           only instances whose name contains TEXT\n"
              << "  --algos a,b             pruned,sam-iterative,sam-recursive,hcp-as-tsp (default pruned)\n"
              << "  --threads 1,2,4         OpenMP thread counts\n"
              << "  --ranks 1               MPI rank counts\n"
              << "  --repeat N              repetitions per configuration, median time is reported\n"
              << "  --time-limit S          per-run limit in seconds (default 60)\n"
              << "  --frontier F            partial paths per thread\n"
              << "  --out FILE              CSV output (default results/bench.csv)\n"
              << "  --summary FILE          Markdown summary (default next to the CSV)\n"
              << "  --dry-run               print commands only\n";
}

bool parseArgs(int argc, char** argv, Options& o) {
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value for " + a);
            return argv[++i];
        };
        if (a == "-h" || a == "--help") return false;
        else if (a == "--bin") o.bin = value();
        else if (a == "--mpi-bin") o.mpiBin = value();
        else if (a == "--launcher") o.launcher = value();
        else if (a == "--manifest") o.manifests.push_back(value());
        else if (a == "--suite") o.suites = splitWords(value());
        else if (a == "--filter") o.filter = value();
        else if (a == "--algos") o.algos = splitWords(value());
        else if (a == "--threads") o.threads = splitInts(value());
        else if (a == "--ranks") o.ranks = splitInts(value());
        else if (a == "--repeat") o.repeat = std::max(1, std::stoi(value()));
        else if (a == "--time-limit") o.timeLimit = std::stod(value());
        else if (a == "--frontier") o.frontier = std::stoi(value());
        else if (a == "--out") o.out = value();
        else if (a == "--summary") o.summary = value();
        else if (a == "--dry-run") o.dryRun = true;
        else throw std::runtime_error("unknown option " + a);
    }
    if (o.manifests.empty()) o.manifests.push_back("bench/instances.csv");
    return true;
}

std::vector<Instance> loadManifests(const Options& o) {
    std::vector<Instance> list;
    bool all = std::find(o.suites.begin(), o.suites.end(), "all") != o.suites.end();
    for (const auto& path : o.manifests) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot open manifest " + path);
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)) {
            line = trimLine(line);
            if (line.empty()) continue;
            auto f = split(line, ',');
            if (f.size() < 6) continue;
            Instance inst{f[0], f[1], f[2], f[3], f[4], f[5], ""};
            for (size_t k = 6; k < f.size(); k++) inst.source += (k > 6 ? "," : "") + f[k];
            if (!all && std::find(o.suites.begin(), o.suites.end(), inst.suite) == o.suites.end()) continue;
            if (!o.filter.empty() && inst.name.find(o.filter) == std::string::npos) continue;
            list.push_back(inst);
        }
    }
    return list;
}

bool fileExists(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0;
}

std::string buildCommand(const Options& o, const Instance& inst, const std::string& algo, int ranks, int threads) {
    std::ostringstream cmd;
    cmd << "OMP_NUM_THREADS=" << threads << " ";
    if (ranks > 1) {
        std::string l = o.launcher;
        size_t pos = l.find("{ranks}");
        if (pos != std::string::npos) l.replace(pos, 7, std::to_string(ranks));
        cmd << l << " " << quote(o.mpiBin);
    } else {
        cmd << quote(o.bin);
    }
    cmd << " " << quote(inst.file) << " --mode " << inst.mode << " --threads " << threads << " --time-limit "
        << o.timeLimit << " --frontier " << o.frontier << " --label " << quote(inst.name) << " --quiet";
    bool batch = inst.args.find("--batch") != std::string::npos;
    if (inst.mode == "path" || inst.mode == "cycle") {
        cmd << " --algo " << algo;
        if (algo.rfind("sam", 0) == 0 && inst.mode == "path" && !batch) cmd << " --all-starts";
    }
    if (!inst.expect.empty()) cmd << " --expect " << inst.expect;
    if (!inst.args.empty()) cmd << " " << inst.args;
    cmd << " 2>&1";
    return cmd.str();
}

bool runCommand(const std::string& cmd, std::map<std::string, std::string>& result, double& wall, std::string& error) {
    auto t0 = std::chrono::steady_clock::now();
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        error = "popen failed";
        return false;
    }
    std::string last, all;
    char buffer[4096];
    while (fgets(buffer, sizeof(buffer), pipe)) {
        std::string line = trimLine(buffer);
        all += line + "\n";
        if (line.rfind("RESULT,", 0) == 0) last = line;
    }
    pclose(pipe);
    wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (last.empty()) {
        error = all.size() > 400 ? all.substr(all.size() - 400) : all;
        return false;
    }
    auto values = split(last.substr(7), ',');
    for (size_t k = 0; k < RESULT_FIELDS.size() && k < values.size(); k++) result[RESULT_FIELDS[k]] = values[k];
    return true;
}

bool solvedStatus(const std::string& status) {
    if (status.rfind("BATCH(", 0) == 0) return status.find("timeout=0") != std::string::npos;
    return status == "FOUND" || status == "NONE" || status == "OPTIMAL";
}

void makeParent(const std::string& path) {
    size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) return;
    std::string dir = path.substr(0, slash);
    std::string cur;
    for (auto& part : split(dir, '/')) {
        cur += part + "/";
        if (!part.empty()) mkdir(cur.c_str(), 0755);
    }
}

std::string fmt(double v, int digits) {
    if (v < 0) return "-";
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(digits) << v;
    return ss.str();
}

void writeCsv(const std::string& path, const std::vector<Row>& rows) {
    makeParent(path);
    std::ofstream out(path);
    out << "suite,instance,mode,algo,algorithm,ranks,threads,cores,vertices,edges,status,cost,expected,verdict,"
           "seconds,wall_seconds,nodes,lower_bound,verified,speedup,efficiency,source\n";
    for (const auto& r : rows) {
        auto g = [&](const std::string& k) {
            auto it = r.result.find(k);
            return it == r.result.end() ? std::string() : it->second;
        };
        out << r.inst.suite << "," << r.inst.name << "," << r.inst.mode << "," << r.algo << "," << g("algorithm")
            << "," << r.ranks << "," << r.threads << "," << r.ranks * r.threads << "," << g("vertices") << ","
            << g("edges") << "," << g("status") << "," << g("cost") << "," << g("expected") << "," << g("verdict")
            << "," << fmt(r.seconds, 6) << "," << fmt(r.wall, 4) << "," << g("nodes") << "," << g("lower_bound")
            << "," << g("verified") << "," << fmt(r.speedup, 4) << "," << fmt(r.efficiency, 4) << ","
            << r.inst.source << "\n";
    }
}

void writeSummary(const std::string& path, const std::vector<Row>& rows) {
    makeParent(path);
    std::ofstream out(path);
    out << "# Benchmark summary\n\n";
    out << "Every answer is compared with the published result for its public instance "
           "(Hamiltonicity from the Flinders sets, optimal tour length from TSPLIB).\n\n";
    out << "## Correctness against published answers\n\n";
    out << "| suite | runs | MATCH | MISMATCH | UNSOLVED / GAP / UNKNOWN | INVALID |\n|---|---|---|---|---|---|\n";
    std::vector<std::string> order;
    std::map<std::string, std::map<std::string, int>> counts;
    for (const auto& r : rows) {
        if (!counts.count(r.inst.suite)) order.push_back(r.inst.suite);
        std::string v = r.result.count("verdict") ? r.result.at("verdict") : "ERROR";
        std::string key = v == "MATCH" || v == "MISMATCH" || v == "INVALID" ? v : "OTHER";
        counts[r.inst.suite][key]++;
        counts[r.inst.suite]["runs"]++;
    }
    for (const auto& s : order) {
        auto& c = counts[s];
        out << "| " << s << " | " << c["runs"] << " | " << c["MATCH"] << " | " << c["MISMATCH"] << " | " << c["OTHER"]
            << " | " << c["INVALID"] << " |\n";
    }
    for (const auto& s : order) {
        out << "\n## Suite `" << s << "`\n\n";
        out << "| instance | n | published | algo | ranks x threads | status | cost | verdict | median s | speedup | efficiency" ;
        for (const auto& r : rows) {
            if (r.inst.suite != s) continue;
            auto g = [&](const std::string& k) {
                auto it = r.result.find(k);
                return it == r.result.end() ? std::string("-") : it->second;
            };
            std::string cost = g("cost") == "-1" ? "-" : g("cost");
            out << "| " << r.inst.name << " | " << g("vertices") << " | " << (r.inst.expect.empty() ? "-" : r.inst.expect)
                << " | " << r.algo << " | " << r.ranks << " x " << r.threads << " | " << g("status") << " | " << cost
                << " | " << g("verdict") << " | " << fmt(r.seconds, 4) << " | " << fmt(r.speedup, 2) << " | "
                << fmt(r.efficiency, 2) << " |\n";
        }
    }
    std::map<std::string, std::set<std::string>> optimal;
    for (const auto& r : rows)
        if (r.result.count("status") && r.result.at("status") == "OPTIMAL") optimal[r.inst.name].insert(r.result.at("cost"));
    std::vector<std::string> disagree;
    for (auto& kv : optimal)
        if (kv.second.size() > 1) disagree.push_back(kv.first);
    out << "\n## Consistency\n\n";
    if (disagree.empty()) out << "All configurations agree on every proven optimum.\n";
    else {
        out << "Configurations disagree on:";
        for (auto& d : disagree) out << " " << d;
        out << "\n";
    }
}

}

int main(int argc, char** argv) {
    Options o;
    try {
        if (!parseArgs(argc, argv, o)) {
            usage();
            return 1;
        }
        std::vector<Instance> instances = loadManifests(o);
        if (instances.empty()) {
            std::cout << "no instances selected" << std::endl;
            return 1;
        }
        std::vector<Row> rows;
        std::vector<std::string> missing;
        for (const auto& inst : instances) {
            if (!fileExists(inst.file)) {
                missing.push_back(inst.name);
                continue;
            }
            bool graphMode = inst.mode == "path" || inst.mode == "cycle";
            bool batch = inst.args.find("--batch") != std::string::npos;
            std::vector<std::string> algos = graphMode ? o.algos : std::vector<std::string>{"bnb"};
            size_t firstRow = rows.size();
            for (const auto& algo : algos)
                for (int r : o.ranks)
                    for (int t : o.threads) {
                        if (algo.rfind("sam", 0) == 0 && (r > 1 || t > 1) && !batch) continue;
                        std::string cmd = buildCommand(o, inst, algo, r, t);
                        if (o.dryRun) {
                            std::cout << cmd << std::endl;
                            continue;
                        }
                        std::vector<std::pair<double, std::map<std::string, std::string>>> runs;
                        double wall = 0.0;
                        for (int k = 0; k < o.repeat; k++) {
                            std::map<std::string, std::string> res;
                            std::string error;
                            double w = 0.0;
                            if (!runCommand(cmd, res, w, error)) {
                                res["status"] = "ERROR";
                                res["verdict"] = "ERROR";
                                res["seconds"] = std::to_string(w);
                                std::cout << "  error: " << error << std::endl;
                            }
                            wall = w;
                            runs.push_back({std::stod(res.count("seconds") ? res["seconds"] : "0"), res});
                        }
                        std::sort(runs.begin(), runs.end(),
                                  [](const auto& a, const auto& b) { return a.first < b.first; });
                        Row row;
                        row.inst = inst;
                        row.algo = algo;
                        row.ranks = r;
                        row.threads = t;
                        row.result = runs[runs.size() / 2].second;
                        row.seconds = runs[runs.size() / 2].first;
                        row.wall = wall;
                        row.command = cmd;
                        rows.push_back(row);
                        std::cout << std::left << std::setw(10) << inst.suite << std::setw(30) << inst.name.substr(0, 29)
                                  << std::setw(15) << algo << std::setw(7)
                                  << (std::to_string(r) + "x" + std::to_string(t)) << std::setw(38)
                                  << row.result["status"].substr(0, 37) << std::setw(11) << row.result["verdict"]
                                  << std::right << std::setw(11) << fmt(row.seconds, 4) << "s" << std::endl;
                    }
            for (size_t i = firstRow; i < rows.size(); i++) {
                Row& row = rows[i];
                const Row* base = nullptr;
                for (size_t j = firstRow; j < rows.size(); j++)
                    if (rows[j].algo == row.algo && rows[j].ranks == 1 && rows[j].threads == 1) base = &rows[j];
                if (base && solvedStatus(base->result.at("status")) && solvedStatus(row.result["status"]) &&
                    row.seconds > 0) {
                    row.speedup = base->seconds / row.seconds;
                    row.efficiency = row.speedup / (row.ranks * row.threads);
                }
            }
        }
        if (o.dryRun) return 0;
        if (!missing.empty()) {
            std::cout << "skipped " << missing.size() << " instance(s) with no data file:";
            for (auto& m : missing) std::cout << " " << m;
            std::cout << std::endl;
        }
        writeCsv(o.out, rows);
        std::string summary = o.summary;
        if (summary.empty()) {
            size_t dot = o.out.find_last_of('.');
            summary = (dot == std::string::npos ? o.out : o.out.substr(0, dot)) + "_summary.md";
        }
        writeSummary(summary, rows);
        std::cout << "wrote " << o.out << " and " << summary << std::endl;
        int bad = 0;
        for (auto& r : rows)
            if (r.result["verdict"] == "MISMATCH" || r.result["verdict"] == "INVALID") bad++;
        if (bad) {
            std::cout << bad << " run(s) disagree with the published answer" << std::endl;
            return 2;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
