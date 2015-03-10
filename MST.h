#include "common.h"
#include "Minmatching/PerfectMatching.h"

#pragma once

class MST {
public:
	int** adjacentMatrix;
	int* parent; //Array to store constructed MST
	int* key; //Key values used to pick minimum weight edge in cut
	bool* mstSet; //To represent set of vertices not yet included in MST
	int N; //the size of pointset
	int last_vertex;
	int* tsp2;
	bool* visited; //Array for TSP implemeneation for 2-approximation

	MST(int** adjacentMatrix, int size);
	~MST();

	//deliverable a
	void makeTree();
	void printMST();

	//deliverable b
	void makeTSP2();
	int findValue(int padre);

	//deliverable c
	void makeTSP1_5(PerfectMatching * pm, float cost, int node_num, int * oddDegree);
	
	float calMean(int option);
	float calStd(int option);
	void combine(PerfectMatching * pm, float cost, int node_num, int * oddDegree);


private:
	void minimumMatching();
	int minKey(int key[], bool mstSet[]);

};