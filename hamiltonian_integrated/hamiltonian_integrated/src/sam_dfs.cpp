#include "sam_dfs.h"

#include <stack>

namespace {

struct StackEntry {
    int v;
    int depth;
};

bool closes(const Graph& g, bool cycle, const std::vector<int>& path) {
    if (!cycle) return true;
    if (g.V <= 1) return true;
    return g.hasEdge(path.back(), path.front());
}

/*
 * Uses iterative DFS with a stack to find Hamiltonian paths.
 */
int iterativeFrom(const Graph& graph, bool cycle, int start, const Deadline& deadline,
                  long long& nodes, std::vector<int>& out) {
//modified, original code from: https://youcademy.org/graph-dfs-implementation-cpp#iterative-dfs-implementation-using-a-stack
    std::vector<char> visited(graph.V, 0); // Keep track of visited nodes
    std::vector<int> path;
    std::stack<StackEntry> stack; // Use a stack
    stack.push({start, 0}); // Add the starting node
    while (!stack.empty()) {
        StackEntry current = stack.top(); // Get the top node
        stack.pop(); // Remove the top node
        while ((int)path.size() > current.depth) //we have gone up the tree, so we must undo some of our choices
        {
            visited[path.back()] = 0; //mark as unvisited
            path.pop_back(); //remove from path
        }
        if (visited[current.v]) continue;
        visited[current.v] = 1; //we have not visited this node already
        path.push_back(current.v); //add it to the path
        nodes++;
        if ((nodes & 4095) == 0 && deadline.expired()) return -1;
        if ((int)path.size() == graph.V) //we have found a hamiltonian path
        {
            if (closes(graph, cycle, path)) {
                out = path; //save path
                return 1; //if we're just looking for the first hamiltonian path, not all, we can stop here
            }
            continue;
        }
        // Add neighbors to the stack in reverse order to maintain DFS order
        const std::vector<int>& neighbors = graph.adjList[current.v];
        for (int i = (int)neighbors.size() - 1; i >= 0; i--)
            if (!visited[neighbors[i]]) stack.push({neighbors[i], current.depth + 1}); //if we have not visited the node already, add it
    }
    return 0;
}

struct RecursiveState {
    const Graph& graph;
    bool cycle;
    const Deadline& deadline;
    std::vector<char> visited;
    std::vector<int> path;
    long long nodes = 0;
    bool timedOut = false;
};

/*
 * Recursive method: Uses recursive DFS to find Hamiltonian paths.
 */
bool recursiveMain(RecursiveState& s, int node) {
    //modified, original code from: https://youcademy.org/graph-dfs-implementation-cpp
    s.visited[node] = 1; //add to path
    s.path.push_back(node); //add to visited
    s.nodes++;
    if ((s.nodes & 4095) == 0 && s.deadline.expired()) s.timedOut = true;
    if (!s.timedOut) {
        if ((int)s.path.size() == s.graph.V) // base case: we have reached a leaf node
        {
            if (closes(s.graph, s.cycle, s.path)) return true;
        } else {
            for (int neighbor : s.graph.adjList[node]) //process neighbors
                if (!s.visited[neighbor] && recursiveMain(s, neighbor)) return true; //we have not visited the node already //TODO: do we need this if condition?
        }
    }
    s.visited[node] = 0;
    s.path.pop_back();
    return false;
}

}

SearchResult dfsIterative(const Graph& graph, bool cycle, int startNode, bool allStarts, const Deadline& deadline) {
    SearchResult result;
    result.algorithm = "sam-iterative";
    double t0 = wallTime();
    int first = allStarts && !cycle ? 0 : startNode;
    int last = allStarts && !cycle ? graph.V - 1 : startNode;
    for (int s = first; s <= last && graph.V > 0; s++) {
        int r = iterativeFrom(graph, cycle, s, deadline, result.nodes, result.path);
        if (r == 1) { result.status = SearchStatus::Found; break; }
        if (r == -1) { result.status = SearchStatus::Timeout; break; }
    }
    result.seconds = wallTime() - t0;
    return result;
}

/*
 * Wrapper method: Uses recursive DFS to find Hamiltonian paths.
 */
SearchResult dfsRecursive(const Graph& graph, bool cycle, int startNode, bool allStarts, const Deadline& deadline) {
    SearchResult result;
    result.algorithm = "sam-recursive";
    double t0 = wallTime();
    RecursiveState s{graph, cycle, deadline, std::vector<char>(graph.V, 0), {}};
    int first = allStarts && !cycle ? 0 : startNode;
    int last = allStarts && !cycle ? graph.V - 1 : startNode;
    for (int st = first; st <= last && graph.V > 0; st++) {
        if (recursiveMain(s, st)) {
            result.status = SearchStatus::Found;
            result.path = s.path;
            break;
        }
        if (s.timedOut) { result.status = SearchStatus::Timeout; break; }
    }
    result.nodes = s.nodes;
    result.seconds = wallTime() - t0;
    return result;
}
