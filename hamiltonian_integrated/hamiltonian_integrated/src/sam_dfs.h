#ifndef HAM_SAM_DFS_H
#define HAM_SAM_DFS_H

#include "common.h"
#include "graph.h"

SearchResult dfsIterative(const Graph& graph, bool cycle, int startNode, bool allStarts, const Deadline& deadline);
SearchResult dfsRecursive(const Graph& graph, bool cycle, int startNode, bool allStarts, const Deadline& deadline);

#endif
