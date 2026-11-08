#include <iostream>
#include <stack>
#include <vector>
#include "graph.h"
#include "solver.h"

using namespace std;

struct node {
    int v;
    int parent;
};

/*
 * Uses iterative DFS with a stack to find Hamiltonian paths.
 */
vector<vector<int>> dfs(Graph& graph, int startNode, bool findAll) {
//modified, original code from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack
    unordered_set<int> visited; // Keep track of visited nodes
    stack<node> stack; // Use a stack

    vector<int> path;
    vector<vector<int>> allPaths;

    node start = { startNode, -1 };
    stack.push(start); // Add the starting node
    int prev = -1;

    while (!stack.empty()) {
        node current = stack.top(); // Get the top node
        stack.pop(); // Remove the top node

        if (current.parent != prev) //we have gone up the tree, so we must undo some of our choices
        {
            int choice = path.back();
            int count = 0;
            while (choice != current.parent)
            {
                visited.erase(choice); //mark as unvisited
                path.pop_back(); //remove from path
                choice = path.back();
                count++;
            }
        }

        prev = current.v;

        if (visited.find(current.v) == visited.end()) { //we have not visited this node already
            visited.insert(current.v);
            //cout << current.v << " "; // Process the node
            path.push_back(current.v); //add it to the path

            // Add neighbors to the stack in reverse order to maintain DFS order
            vector<int> neighbors = graph.adjList[current.v];
            int neighborsAdded = 0;
            for (int i = neighbors.size() - 1; i >= 0; i--) {
                int neighbor = neighbors[i];
                if (visited.find(neighbor) == visited.end()) { //if we have not visited the node already, add it
                    node n = { neighbor, current.v };
                    stack.push(n);
                    neighborsAdded++;
                }
            }
            if (neighborsAdded == 0) //leaf node in the dfs tree
            {
                if (path.size() == graph.V) //we have found a hamiltonian path
                {
                    //save path
                    allPaths.push_back(path);
                    if (!findAll) //if we're just looking for the first hamiltonian path, not all, we can stop here
                        break;
                }
            }
        }
    }

    return allPaths;

    cout << endl;
}