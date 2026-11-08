/*
	Usage: ./hamiltonian
*/

#include "graph.h"
#include "setup.h"
#include <iostream>
#include "solver.h"
#include <sys/time.h>

using namespace std;

//code for time from IT388 Lect 22 code
double getTime() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec + t.tv_usec / 1000000.0;

}

int main(int argc, char* argv[]) {
	double start, elapsed;

	//process file and convert to graph
	Graph g = convertHCPFileToGraph("TestData/tsphcp/SH_64.hcp", true);
	//Graph g = convertHCPFileToGraph("tsphcp/SSP_2007.hcp", true);
	//Graph g = convertHCPFileToGraph("tsphcp/SNm_124.hcp", true);
	//Graph g = convertHCPFileToGraph("test.txt", true);
	//Graph g = convertTxtFileToGraph("TestData/18_H/18_H.txt", true);
	//Graph g = convertTxtFileToGraph("TestData/20_NH/20_NH.txt", true);

	//run the algorithm to get the path(s)
	vector<vector<int>> allPaths;
	start = getTime();
	allPaths = dfs(g, 1, false);
	elapsed = getTime() - start;
	cout << "Time: " << elapsed << " seconds" << endl;

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