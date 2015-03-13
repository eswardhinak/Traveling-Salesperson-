#include "common.h"
#include "Minmatching/PerfectMatching.h"
#include <string>
#include <unordered_map>
#include "Node.hpp"
#pragma once

class MST {
public:
	float** adjacentMatrix;
	int* parent; //Array to store constructed MST
	int* key; //Key values used to pick minimum weight edge in cut
	bool* mstSet; //To represent set of vertices not yet included in MST
	int N; //the size of pointset
	int last_vertex;
	int* tsp2;
	int total_edges; 
	bool* visited; //Array for TSP implemeneation for 2-approximation
	std::unordered_map<int, Node*> nodes; //map from Node number to pointer to heap location
	float mst_total_weight;
	float tsp2_total_weight;
	float tsp1_5_total_weight;


	MST(float** adjacentMatrix, int size);
	~MST();
	void addEdge(int from, int to, unsigned int cost);
	unsigned int totalEdgeCost() const;
	void constructTree();

	//deliverable a
	void makeTree();
	void printMST();

	//deliverable b
	void makeTSP2();
	int findValue(int padre);

	//deliverable c
	void makeTSP1_5(int total_edges);

	float calMean(int option);
	float calStd(int option);
	void combine(PerfectMatching * pm, float cost, int node_num, int * oddDegree);
	int DFSCount(int v, bool * encountered);
	int isValidEdge(int u, int v);
	void unHideEdge(int u, int v);
	void hideEdge(int u, int v);
	void eulerTraversal(int current, int *  ci,  int * traversal);




private:
	void minimumMatching();
	int minKey(int key[], bool mstSet[]);

};