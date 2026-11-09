#ifndef SOLVER
#define SOLVER

#include <vector>

using namespace std;

vector<vector<int>> dfsIterative(Graph& graph, int startNode, bool findAll);

vector<vector<int>> dfsRecursive(Graph& graph, int startNode, bool findAll);

void dfsRecursiveMain(Graph& graph, int node, std::unordered_set<int>& visited, vector<int> path, bool findAll);

#endif