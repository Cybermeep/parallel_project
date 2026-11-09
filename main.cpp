/*
	Usage: ./hamiltonian <input file> <directed (optional)>
*/

#include "graph.h"
#include "setup.h"
#include <iostream>
#include "solver.h"
#include <sys/time.h>

using namespace std;

//code for time from IT388 Lect 22 code
//returns ms
double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + t.tv_usec / 1000.0;

}

int main(int argc, char* argv[]) {
	double start, elapsed;

	char* input_file = argv[1];
	bool directed = true;
	if (argv[2] != NULL)
		directed = argv[2];
	//process file and convert to graph
	Graph g = convertHCPFileToGraph(input_file, directed);
	cout << "File converted to graph" << endl;

	//run the algorithm to get the path(s)
	vector<vector<int>> allPaths;
	start = getTime();
	//allPaths = dfsIterative(g, 1, false);
	allPaths = dfsRecursive(g, 1, false);
	elapsed = getTime() - start;
	cout << "Time: " << elapsed << " ms" << endl;

	//print path(s)
	if (allPaths.size() == 0) { //no hamiltonian paths were found
		cout << "No Hamiltonian paths." << endl;
	}
	else {
		for (int i = 0; i < allPaths.size(); i++)
		{
			cout << endl;
			cout << allPaths[i][0];
			for (int j = 1; j < allPaths[i].size(); j++)
				cout << " -> " << allPaths[i][j];
			cout << endl;
		}
	}

}