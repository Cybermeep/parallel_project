#include "setup.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

Graph convertFileToGraph(char* filename)
{
    fstream fileStream;
    fileStream.open(filename);
    cout << filename << "\n";

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
    cout << "graph created with " << numV << " vertices\n";

    int u, v;
    while (fileStream >> u >> v)
    {
        //cout << u << " " << v << "\n";
        g.addEdge(u, v);
    }

    return g;
}