#include <iostream>
#include <stack>
#include <vector>
#include "graph.h"
#include "solver.h"

void dfs(Graph& graph, int startNode) {
//this code is from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack
    std::unordered_set<int> visited; // Keep track of visited nodes
    std::stack<int> stack; // Use a stack

    stack.push(startNode); // Add the starting node

    std::cout << "DFS Traversal (Iterative): ";
    while (!stack.empty()) {
        int node = stack.top(); // Get the top node
        stack.pop(); // Remove the top node

        if (visited.find(node) == visited.end()) {
            visited.insert(node);
            std::cout << node << " "; // Process the node

            // Add neighbors to the stack in reverse order to maintain DFS order
            std::vector<int> neighbors = graph.adjList[node];
            for (int i = neighbors.size() - 1; i >= 0; i--) {
                char neighbor = neighbors[i];
                if (visited.find(neighbor) == visited.end()) {
                    stack.push(neighbor);
                }
            }
        }
    }
    std::cout << std::endl;
}