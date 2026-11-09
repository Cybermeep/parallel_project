#include "setup.h"
#include <fstream>
#include <iostream>
#include <string>
#include <regex>
#include <sstream>

using namespace std;

/*
 * Used to convert a .hcp file from the Flinders TSP and HCP Benchmark set into a graph.
 */
Graph convertHCPFileToGraph(char* filename, bool directed)
{
    fstream fileStream;
    fileStream.open(filename);

    int numV;
    string line;
    //skip first 3 lines
    getline(fileStream, line);
    getline(fileStream, line);
    getline(fileStream, line);
    //get size
    getline(fileStream, line, ':');
    fileStream >> numV;
    //skip next 2 lines
    getline(fileStream, line);
    getline(fileStream, line);
    getline(fileStream, line);

    Graph g = Graph(numV);

    int u, v;
    while (fileStream >> u >> v)
    {
        g.addEdge(u, v, directed);
    }

    fileStream.close();

    return g;
}

/*
 * Used to convert a .txt file from the Flinders Graph Database into a graph. //TODO: this just gets the first graph right now, the file has multiple
 */
Graph convertTxtFileToGraph(char* filename, bool directed)
{
    fstream fileStream;
    fileStream.open(filename);

    string str = filename;
    regex pattern("(\\d+)");
    smatch matches;
    regex_search(str, matches, pattern);
    int numV = stoi(matches[0].str());

    Graph g = Graph(numV);

    int u, v;
    string line;
    while (getline(fileStream, line))
    {
        if (line.length() < 2)
            break; //TODO: reset graph, add to a list
        stringstream ss(line);
        ss >> u >> v;
        g.addEdge(u, v, directed);
    }

    fileStream.close();

    return g;
}