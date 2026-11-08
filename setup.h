#include "graph.h"
#include <fstream>
#include <string>

#ifndef SETUP
#define SETUP

Graph convertHCPFileToGraph(char* filename, bool directed);

Graph convertTxtFileToGraph(char* filename, bool directed);

#endif