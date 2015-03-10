#include "common.h"
#include "Point.h"
#include <ctime>
#include "MST.h"
#include "Minmatching/PerfectMatching.h"

/*
This project is a starter code and wrappers for CSE101W15 Implementation project.

point.h - uniform random pointset generator

MST.h - minimum spanning tree

PerfectMatching.h - interface to min cost perfect matching code 

-------------------------------------
PerfectMatching is from the paper:

Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

sourcecode : pub.ist.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz

*/

void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, int** adjacentMatrix, int N, int * oddDegree) {
	int e = 0;
	edge_num = node_num*(node_num-1)/2 ; //complete graph
	edges = new int[2*edge_num];
	weights = new int[edge_num];
	for(int i = 0; i < node_num ; ++i) {
		for(int j = i+1 ; j< node_num; ++j) {
			edges[2*e] = i;
			edges[2*e+1] = j;
			weights[e] = adjacentMatrix[oddDegree[i]][oddDegree[j]];
			//cerr << "Weight" << e << ": " << weights[e] << endl;
 			e++;
		}
	}
	if (e != edge_num) {
		cout<<"the number of edge is wrong"<<endl;
		exit(1);
	}
}

void PrintMatching(int node_num, PerfectMatching* pm, int * oddDegree, int ** adjacentMatrix) {
	int i, j;

	for (i=0; i<node_num; i++) {
		j = pm->GetMatch(i);
		if (i < j) printf("%d %d %d\n", oddDegree[i], oddDegree[j], adjacentMatrix[oddDegree[i]][oddDegree[j]] );
	}
}

void calculateDegree(int * degree, MST* mst, int * node_num, int N){
	for (int j = 0; j < N; j++){
		if (mst->parent[j] > -1)
			degree[mst->parent[j]]++;
	}
	for (int i=0; i<N; i++){
		if (i!=0)
			degree[i]++;
		if (degree[i] % 2 != 0)
			(*node_num)++;
		cerr << i << "--->" << degree[i] << endl;
	}

}


int main() {
	cerr << "Program start.";
	set< pair<int,int> > generatedPointset;
	int** adjacentMatrix;
	int W, H, N;
	Point pointset;

	W = 1000;
	H = 1000;
	N = 100;

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal

	MST * mst2 = new MST(adjacentMatrix, N);
	mst2->makeTree();
	//mst2->printMST();
	cout << "Mean: " << mst2->calMean(MST_1) << endl;
	cout << "Standard Deviation: " << mst2->calStd(MST_1)<<endl;

	//Deliverable B: Find TSP2 path from the constructed MST

	clock_t begin = clock();
	mst2->makeTSP2();
	clock_t end = clock();
	double elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
	cerr << "TSP_time:  " << elapsed_seconds << endl;
	cout << "Mean: "<< mst2->calMean(TSP2) << endl;
	cout << "Standard Deviation: " << mst2->calStd(TSP2) << endl;


	//Deliverable C: Find TSP1.5 path from the constructed MST
	struct PerfectMatching::Options options;
	int i, e, node_num = 0;
	int * nn_ptr = &node_num;
	int * degree =  new int[N];

	calculateDegree(degree, mst2, nn_ptr, N);
	int * oddDegree = new int[node_num];

	int j = 0;
	for (int i = 0; i<N;i++){
		if(degree[i]%2 !=0){
			oddDegree[j] = i;
			j++;
		}
	}
	int edge_num = node_num * (node_num-1) / 2;
	int* edges;
	int* weights;
	PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
	LoadInput(node_num, edge_num, edges, weights, adjacentMatrix, N, oddDegree);
	for (e=0; e<edge_num; e++) {
		pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
		//cerr << edges[2*e] << "   "<< edges[2*e+1]<< endl;
	}

	pm->options = options;
	pm->Solve();

	float cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
	printf("Total cost of the perfect min-weight matching = %.1f\n", cost);
	mst2->makeTSP1_5(pm, cost, node_num, oddDegree);

	//PrintMatching(node_num, pm, oddDegree ,adjacentMatrix);

	delete pm;
	delete [] oddDegree;
	delete [] degree;
	delete [] edges;
	delete [] weights;
	
	
	return 0;
}

