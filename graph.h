#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef GRAPH
#define GRAPH

class Graph {
//this code is from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack
public:
    int V;
    int E;
    std::unordered_map<int, std::vector<int>> adjList;

    Graph(int numVertices);
    
    void addVertex(int vertex);

    void addEdge(int source, int destination, bool directed);
};

#endif