#include "graph.h"

#include <algorithm>

//this code is from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack

Graph::Graph(int numVertices, bool isDirected)
    : V(numVertices), E(0), directed(isDirected), adjList(numVertices),
      words((numVertices + 63) / 64), bits((size_t)numVertices * ((numVertices + 63) / 64), 0) {}

void Graph::addEdge(int source, int destination) {
    if (source == destination || source < 0 || destination < 0 || source >= V || destination >= V) return;
    adjList[source].push_back(destination);
    bits[(size_t)source * words + destination / 64] |= (uint64_t)1 << (destination % 64);
    if (!directed) //for undirected graphs, we add the edge going back the other way
    {
        adjList[destination].push_back(source);
        bits[(size_t)destination * words + source / 64] |= (uint64_t)1 << (source % 64);
    }
}

void Graph::finalize() {
    E = 0;
    for (auto& list : adjList) {
        std::sort(list.begin(), list.end());
        list.erase(std::unique(list.begin(), list.end()), list.end());
        E += (long long)list.size();
    }
    if (!directed) E /= 2;
}

bool Graph::hasEdge(int u, int v) const {
    return (bits[(size_t)u * words + v / 64] >> (v % 64)) & 1ULL;
}

bool Graph::isConnected() const {
    if (V <= 1) return true;
    std::vector<std::vector<int>> undirectedAdj;
    const std::vector<std::vector<int>>* adj = &adjList;
    if (directed) {
        undirectedAdj.assign(V, {});
        for (int u = 0; u < V; u++)
            for (int v : adjList[u]) {
                undirectedAdj[u].push_back(v);
                undirectedAdj[v].push_back(u);
            }
        adj = &undirectedAdj;
    }
    std::vector<char> seen(V, 0);
    std::vector<int> stack = {0};
    seen[0] = 1;
    int count = 1;
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (int v : (*adj)[u])
            if (!seen[v]) {
                seen[v] = 1;
                count++;
                stack.push_back(v);
            }
    }
    return count == V;
}

WeightedGraph::WeightedGraph(int numVertices) : n(numVertices), w((size_t)numVertices * numVertices, 0) {}

Graph WeightedGraph::toGraph(long long missingWeight) const {
    Graph g(n, false);
    g.name = name;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (at(i, j) < missingWeight) g.addEdge(i, j);
    g.finalize();
    return g;
}
