#include "io.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

static std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) a++;
    while (b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

static std::string upper(std::string s) {
    for (char& c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

std::string fileExtension(const std::string& path) {
    size_t dot = path.find_last_of('.');
    size_t slash = path.find_last_of("/\\");
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";
    std::string ext = path.substr(dot + 1);
    for (char& c : ext) c = (char)std::tolower((unsigned char)c);
    return ext;
}

std::string baseName(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string file = slash == std::string::npos ? path : path.substr(slash + 1);
    size_t dot = file.find_last_of('.');
    return dot == std::string::npos ? file : file.substr(0, dot);
}

static bool splitKeyValue(const std::string& line, std::string& key, std::string& value) {
    size_t colon = line.find(':');
    if (colon == std::string::npos) {
        key = upper(trim(line));
        value = "";
        return false;
    }
    key = upper(trim(line.substr(0, colon)));
    value = trim(line.substr(colon + 1));
    return true;
}

static long long nint(double v) { return (long long)(v + 0.5); }

static double geoRadians(double coordinate) {
    const double PI = 3.141592;
    double deg = std::trunc(coordinate);
    double min = coordinate - deg;
    return PI * (deg + 5.0 * min / 3.0) / 180.0;
}

static long long tsplibDistance(const std::string& type, double xi, double yi, double xj, double yj) {
    double dx = xi - xj, dy = yi - yj;
    if (type == "EUC_2D") return nint(std::sqrt(dx * dx + dy * dy));
    if (type == "CEIL_2D") return (long long)std::ceil(std::sqrt(dx * dx + dy * dy));
    if (type == "MAN_2D") return nint(std::fabs(dx) + std::fabs(dy));
    if (type == "MAX_2D") return std::max(nint(std::fabs(dx)), nint(std::fabs(dy)));
    if (type == "ATT") {
        double r = std::sqrt((dx * dx + dy * dy) / 10.0);
        long long t = nint(r);
        return t < r ? t + 1 : t;
    }
    if (type == "GEO") {
        const double RRR = 6378.388;
        double lati = geoRadians(xi), loni = geoRadians(yi);
        double latj = geoRadians(xj), lonj = geoRadians(yj);
        double q1 = std::cos(loni - lonj);
        double q2 = std::cos(lati - latj);
        double q3 = std::cos(lati + latj);
        return (long long)(RRR * std::acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    }
    throw std::runtime_error("unsupported EDGE_WEIGHT_TYPE " + type);
}

static void fillExplicit(WeightedGraph& wg, const std::string& format, const std::vector<long long>& values) {
    int n = wg.n;
    size_t k = 0;
    auto next = [&]() {
        if (k >= values.size()) throw std::runtime_error("EDGE_WEIGHT_SECTION too short");
        return values[k++];
    };
    if (format == "FULL_MATRIX") {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) wg.set(i, j, next());
        return;
    }
    for (int i = 0; i < n; i++) wg.set(i, i, 0);
    if (format == "UPPER_ROW") {
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "LOWER_ROW") {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "UPPER_DIAG_ROW") {
        for (int i = 0; i < n; i++)
            for (int j = i; j < n; j++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "LOWER_DIAG_ROW") {
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "UPPER_COL") {
        for (int j = 0; j < n; j++)
            for (int i = 0; i < j; i++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "LOWER_COL") {
        for (int j = 0; j < n; j++)
            for (int i = j + 1; i < n; i++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "UPPER_DIAG_COL") {
        for (int j = 0; j < n; j++)
            for (int i = 0; i <= j; i++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else if (format == "LOWER_DIAG_COL") {
        for (int j = 0; j < n; j++)
            for (int i = j; i < n; i++) { long long v = next(); wg.set(i, j, v); wg.set(j, i, v); }
    } else {
        throw std::runtime_error("unsupported EDGE_WEIGHT_FORMAT " + format);
    }
}

/*
 * Used to convert a .hcp file from the Flinders TSP and HCP Benchmark set into a graph.
 */
Instance readTsplibFile(const std::string& path, bool directed) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Instance inst;
    inst.name = baseName(path);
    int dimension = -1;
    std::string weightType = "EUC_2D", weightFormat = "FULL_MATRIX", edgeFormat = "EDGE_LIST";
    std::vector<double> xs, ys;
    std::vector<long long> explicitValues;
    std::vector<std::pair<int, int>> edges;
    std::string line;
    while (std::getline(in, line)) {
        std::string key, value;
        splitKeyValue(line, key, value);
        if (key.empty()) continue;
        if (key == "EOF") break;
        if (key == "NAME") continue;
        else if (key == "TYPE") inst.type = upper(value);
        else if (key.rfind("DIMEN", 0) == 0) dimension = std::stoi(value); //get size
        else if (key == "EDGE_WEIGHT_TYPE") weightType = upper(value);
        else if (key == "EDGE_WEIGHT_FORMAT") weightFormat = upper(value);
        else if (key == "EDGE_DATA_FORMAT") edgeFormat = upper(value);
        else if (key == "NODE_COORD_SECTION") {
            if (dimension < 0) throw std::runtime_error("DIMENSION missing before NODE_COORD_SECTION");
            xs.assign(dimension, 0.0);
            ys.assign(dimension, 0.0);
            for (int i = 0; i < dimension; i++) {
                int id;
                double x, y;
                in >> id >> x >> y;
                if (id >= 1 && id <= dimension) { xs[id - 1] = x; ys[id - 1] = y; }
            }
        } else if (key == "EDGE_WEIGHT_SECTION") {
            if (dimension < 0) throw std::runtime_error("DIMENSION missing before EDGE_WEIGHT_SECTION");
            size_t needed = 0;
            if (weightFormat == "FULL_MATRIX") needed = (size_t)dimension * dimension;
            else if (weightFormat.find("DIAG") != std::string::npos) needed = (size_t)dimension * (dimension + 1) / 2;
            else needed = (size_t)dimension * (dimension - 1) / 2;
            long long v;
            while (explicitValues.size() < needed && in >> v) explicitValues.push_back(v);
        } else if (key.rfind("EDGE_DATA_S", 0) == 0) {
            if (edgeFormat == "ADJ_LIST") {
                long long token;
                int current = -1;
                while (in >> token) {
                    if (token == -1) {
                        if (current == -1) break;
                        current = -1;
                    } else if (current == -1) {
                        current = (int)token;
                    } else {
                        edges.push_back({current - 1, (int)token - 1});
                    }
                }
            } else {
                std::string tok;
                while (in >> tok) {
                    if (tok == "-1" || upper(tok) == "EOF") break;
                    std::string tok2;
                    if (!(in >> tok2)) break;
                    edges.push_back({std::stoi(tok) - 1, std::stoi(tok2) - 1});
                }
            }
        } else if (key == "DISPLAY_DATA_SECTION" || key == "TOUR_SECTION" || key == "FIXED_EDGES_SECTION") {
            std::string tok;
            while (in >> tok && tok != "-1" && upper(tok) != "EOF") {}
        }
    }
    if (dimension <= 0) throw std::runtime_error("DIMENSION missing in " + path);
    if (inst.type.empty()) inst.type = edges.empty() ? "TSP" : "HCP";

    if (inst.type == "HCP" || !edges.empty()) {
        inst.weighted = false;
        inst.graph = Graph(dimension, directed);
        inst.graph.name = inst.name;
        for (auto& e : edges) inst.graph.addEdge(e.first, e.second);
        inst.graph.finalize();
        return inst;
    }

    inst.weighted = true;
    WeightedGraph wg(dimension);
    wg.name = inst.name;
    wg.weightType = weightType;
    if (weightType == "EXPLICIT") {
        fillExplicit(wg, weightFormat, explicitValues);
    } else {
        if (xs.empty()) throw std::runtime_error("NODE_COORD_SECTION missing in " + path);
        for (int i = 0; i < dimension; i++)
            for (int j = 0; j < dimension; j++)
                wg.set(i, j, i == j ? 0 : tsplibDistance(weightType, xs[i], ys[i], xs[j], ys[j]));
        wg.x = xs;
        wg.y = ys;
    }
    inst.weights = wg;
    inst.graph = wg.toGraph((long long)4e18);
    return inst;
}

/*
 * Used to convert a .txt file from the Flinders Graph Database into a graph. //TODO: this just gets the first graph right now, the file has multiple
 */
std::vector<Graph> readEdgeListCollection(const std::string& path, bool directed, int forcedVertices) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<Graph> graphs;
    std::vector<std::pair<int, int>> edges;
    int maxVertex = 0;
    auto flush = [&]() {
        if (edges.empty()) return;
        int n = forcedVertices > 0 ? forcedVertices : maxVertex;
        Graph g(n, directed);
        g.name = baseName(path) + "#" + std::to_string(graphs.size() + 1);
        for (auto& e : edges) g.addEdge(e.first - 1, e.second - 1);
        g.finalize();
        graphs.push_back(std::move(g));
        edges.clear();
        maxVertex = 0;
    };
    std::string line;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty()) { flush(); continue; } //TODO: reset graph, add to a list
        std::istringstream ss(t);
        int u, v;
        if (ss >> u >> v) {
            edges.push_back({u, v});
            maxVertex = std::max(maxVertex, std::max(u, v));
        }
    }
    flush();
    return graphs;
}

std::vector<Graph> readGenregAscCollection(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<Graph> graphs;
    std::vector<std::vector<int>> rows;
    std::string title;
    auto flush = [&]() {
        if (rows.empty()) return;
        int n = (int)rows.size();
        Graph g(n, false);
        g.name = baseName(path) + "#" + (title.empty() ? std::to_string(graphs.size() + 1) : title);
        for (int u = 0; u < n; u++)
            for (int v : rows[u]) g.addEdge(u, v - 1);
        g.finalize();
        graphs.push_back(std::move(g));
        rows.clear();
    };
    std::string line;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty()) continue;
        if (upper(t).rfind("GRAPH", 0) == 0) {
            flush();
            std::istringstream ss(t.substr(5));
            std::string id;
            ss >> id;
            if (!id.empty() && id.back() == ':') id.pop_back();
            title = id;
            continue;
        }
        size_t colon = t.find(':');
        if (colon == std::string::npos) continue;
        std::string head = trim(t.substr(0, colon));
        if (head.empty() || !std::all_of(head.begin(), head.end(), [](unsigned char c) { return std::isdigit(c); }))
            continue;
        int u = std::stoi(head);
        std::istringstream ss(t.substr(colon + 1));
        if ((int)rows.size() < u) rows.resize(u);
        int v;
        while (ss >> v) rows[u - 1].push_back(v);
    }
    flush();
    return graphs;
}

std::vector<Graph> readGraph6Collection(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<Graph> graphs;
    std::string line;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.rfind(">>graph6<<", 0) == 0) t = t.substr(10);
        if (t.empty()) continue;
        size_t pos = 0;
        long long n;
        if ((unsigned char)t[0] != 126) {
            n = (unsigned char)t[0] - 63;
            pos = 1;
        } else if ((unsigned char)t[1] != 126) {
            n = 0;
            for (int i = 1; i <= 3; i++) n = (n << 6) | ((unsigned char)t[i] - 63);
            pos = 4;
        } else {
            n = 0;
            for (int i = 2; i <= 7; i++) n = (n << 6) | ((unsigned char)t[i] - 63);
            pos = 8;
        }
        Graph g((int)n, false);
        g.name = baseName(path) + "#" + std::to_string(graphs.size() + 1);
        long long bit = 0;
        for (int j = 1; j < n; j++)
            for (int i = 0; i < j; i++, bit++) {
                size_t byteIndex = pos + (size_t)(bit / 6);
                if (byteIndex >= t.size()) continue;
                int value = (unsigned char)t[byteIndex] - 63;
                if ((value >> (5 - bit % 6)) & 1) g.addEdge(i, j);
            }
        g.finalize();
        graphs.push_back(std::move(g));
    }
    return graphs;
}

std::vector<int> readTourFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<int> tour;
    std::string line;
    bool inSection = false;
    while (std::getline(in, line)) {
        std::string t = upper(trim(line));
        if (!inSection) {
            if (t.rfind("TOUR_SECTION", 0) == 0) inSection = true;
            continue;
        }
        std::istringstream ss(t);
        long long v;
        while (ss >> v) {
            if (v == -1) return tour;
            tour.push_back((int)v - 1);
        }
    }
    return tour;
}

std::vector<Instance> loadInstances(const std::string& path, bool directed, int forcedVertices) {
    std::string ext = fileExtension(path);
    std::vector<Instance> result;
    auto wrap = [&](std::vector<Graph>&& graphs) {
        for (auto& g : graphs) {
            Instance inst;
            inst.name = g.name;
            inst.type = "HCP";
            inst.graph = std::move(g);
            result.push_back(std::move(inst));
        }
    };
    if (ext == "txt" || ext == "el" || ext == "edges") {
        std::ifstream probe(path);
        std::string first;
        while (std::getline(probe, first) && trim(first).empty()) {}
        std::string key, value;
        splitKeyValue(first, key, value);
        if (key == "NAME") result.push_back(readTsplibFile(path, directed));
        else wrap(readEdgeListCollection(path, directed, forcedVertices));
    } else if (ext == "asc") {
        wrap(readGenregAscCollection(path));
    } else if (ext == "g6") {
        wrap(readGraph6Collection(path));
    } else {
        result.push_back(readTsplibFile(path, directed));
    }
    return result;
}
