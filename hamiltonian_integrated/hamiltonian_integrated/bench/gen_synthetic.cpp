#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

namespace {

std::vector<int> parseSizes(const std::string& s) {
    std::vector<int> v;
    std::stringstream ss(s);
    std::string part;
    while (std::getline(ss, part, ','))
        if (!part.empty()) v.push_back(std::stoi(part));
    return v;
}

void writeTsp(const std::string& path, const std::string& name, int n, std::mt19937_64& rng) {
    std::uniform_real_distribution<double> coord(0.0, 1000.0);
    std::ofstream f(path);
    f << "NAME: " << name << "\nTYPE: TSP\nCOMMENT: uniform random points\nDIMENSION: " << n
      << "\nEDGE_WEIGHT_TYPE: EUC_2D\nNODE_COORD_SECTION\n";
    f << std::fixed << std::setprecision(2);
    for (int i = 1; i <= n; i++) {
        double x = coord(rng), y = coord(rng);
        f << i << " " << x << " " << y << "\n";
    }
    f << "EOF\n";
}

void writePlanted(const std::string& path, const std::string& name, int n, double extraDegree, std::mt19937_64& rng) {
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    std::set<std::pair<int, int>> edges;
    for (int i = 0; i < n; i++) {
        int a = order[i], b = order[(i + 1) % n];
        edges.insert({std::min(a, b), std::max(a, b)});
    }
    std::uniform_int_distribution<int> pick(0, n - 1);
    size_t target = (size_t)n + (size_t)(n * extraDegree / 2);
    while (edges.size() < target) {
        int a = pick(rng), b = pick(rng);
        if (a != b) edges.insert({std::min(a, b), std::max(a, b)});
    }
    std::ofstream f(path);
    f << "NAME: " << name << "\nTYPE: HCP\nCOMMENT: random graph with a planted Hamiltonian cycle\nDIMENSION: " << n
      << "\nEDGE_DATA_FORMAT: EDGE_LIST\nEDGE_DATA_SECTION\n";
    for (auto& e : edges) f << e.first + 1 << " " << e.second + 1 << "\n";
    f << "-1\nEOF\n";
}

}

int main(int argc, char** argv) {
    unsigned long long seed = 388;
    std::string tspSizes = "40,50,60,70", hcpSizes = "150,200,300";
    double extraDegree = 2.0;
    std::string outDir = "data/synthetic", manifest = "bench/synthetic.csv";
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        std::string v = i + 1 < argc ? argv[i + 1] : "";
        if (a == "--seed") { seed = std::stoull(v); i++; }
        else if (a == "--tsp-sizes") { tspSizes = v; i++; }
        else if (a == "--hcp-sizes") { hcpSizes = v; i++; }
        else if (a == "--extra-degree") { extraDegree = std::stod(v); i++; }
        else if (a == "--out") { outDir = v; i++; }
        else if (a == "--manifest") { manifest = v; i++; }
        else {
            std::cout << "Usage: gen_synthetic [--seed S] [--tsp-sizes 40,50] [--hcp-sizes 150,200] "
                         "[--extra-degree D] [--out DIR] [--manifest FILE]\n";
            return 1;
        }
    }
    mkdir(outDir.c_str(), 0755);
    std::vector<std::string> rows = {"suite,instance,file,mode,expect,args,source"};
    for (int n : parseSizes(tspSizes)) {
        std::mt19937_64 rng(seed * 1000 + n);
        std::string name = "rand" + std::to_string(n) + "_s" + std::to_string(seed);
        std::string path = outDir + "/" + name + ".tsp";
        writeTsp(path, name, n, rng);
        rows.push_back("synthetic," + name + "," + path +
                       ",tsp,,,synthetic uniform EUC_2D (no published optimum; configurations must agree)");
    }
    for (int n : parseSizes(hcpSizes)) {
        std::mt19937_64 rng(seed * 7919 + n);
        std::ostringstream nm;
        nm << "planted" << n << "_d" << extraDegree << "_s" << seed;
        std::string path = outDir + "/" + nm.str() + ".hcp";
        writePlanted(path, nm.str(), n, extraDegree, rng);
        rows.push_back("synthetic," + nm.str() + "," + path +
                       ",cycle,hamiltonian,,synthetic planted Hamiltonian cycle (Hamiltonian by construction)");
    }
    std::ofstream m(manifest);
    for (auto& r : rows) m << r << "\n";
    std::cout << "wrote " << rows.size() - 1 << " instances to " << outDir << " and manifest " << manifest << std::endl;
    return 0;
}
