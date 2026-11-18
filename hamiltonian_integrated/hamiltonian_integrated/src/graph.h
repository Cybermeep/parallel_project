#ifndef HAM_GRAPH_H
#define HAM_GRAPH_H

#include <cstdint>
#include <string>
#include <vector>

class Graph {
//this code is from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack
public:
    int V;
    long long E;
    bool directed;
    std::string name;
    std::vector<std::vector<int>> adjList;

    explicit Graph(int numVertices = 0, bool isDirected = false);

    void addEdge(int source, int destination);
    void finalize();
    bool hasEdge(int u, int v) const;
    int degree(int v) const { return (int)adjList[v].size(); }
    bool isConnected() const;

private:
    int words;
    std::vector<uint64_t> bits;
};

class WeightedGraph {
public:
    int n;
    std::string name;
    std::string weightType;
    std::vector<long long> w;
    std::vector<double> x;
    std::vector<double> y;

    explicit WeightedGraph(int numVertices = 0);
    long long at(int i, int j) const { return w[(size_t)i * n + j]; }
    void set(int i, int j, long long value) { w[(size_t)i * n + j] = value; }
    Graph toGraph(long long missingWeight) const;
};

#endif
