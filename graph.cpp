#include "graph.h"

//this code is from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack

Graph::Graph(int numVertices) {
    for (int v = 0; v < numVertices; v++)
        adjList[v] = std::vector<int>();
};

void Graph::addVertex(int vertex) {
    adjList[vertex] = std::vector<int>();
};

void Graph::addEdge(int source, int destination) {
    adjList[source].push_back(destination);
    // For undirected graphs, you would also add:
    // adjList[destination].push_back(source);
};