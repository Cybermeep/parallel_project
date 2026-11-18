#ifndef HAM_IO_H
#define HAM_IO_H

#include <string>
#include <vector>

#include "graph.h"

struct Instance {
    std::string name;
    std::string type;
    bool weighted = false;
    Graph graph;
    WeightedGraph weights;
};

std::string fileExtension(const std::string& path);
std::string baseName(const std::string& path);

Instance readTsplibFile(const std::string& path, bool directed);
std::vector<Graph> readEdgeListCollection(const std::string& path, bool directed, int forcedVertices);
std::vector<Graph> readGenregAscCollection(const std::string& path);
std::vector<Graph> readGraph6Collection(const std::string& path);
std::vector<int> readTourFile(const std::string& path);

std::vector<Instance> loadInstances(const std::string& path, bool directed, int forcedVertices);

#endif
