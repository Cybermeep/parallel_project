#include "graph.h"

//this code is from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack

Graph::Graph(int numVertices) {
    V = numVertices;
    E = 0;
    for (int v = 0; v < numVertices; v++)
        adjList[v] = std::vector<int>();
};

void Graph::addVertex(int vertex) {
    adjList[vertex] = std::vector<int>();
};

void Graph::addEdge(int source, int destination, bool directed) {
    this->E++;
    adjList[source].push_back(destination);
    if (!directed) //for undirected graphs, we add the edge going back the other way
    {
        adjList[destination].push_back(source);
        E++;
    }
    
};