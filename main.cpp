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

void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, float** adjacentMatrix, int N, int * oddDegree) {
	int e = 0;
	edge_num = node_num*(node_num-1)/2 ; //complete graph
	edges = new int[2*edge_num];
	weights = new int[edge_num];
	for(int i = 0; i < node_num ; ++i) {
		for(int j = i+1 ; j< node_num; ++j) {
			edges[2*e] = i;
			edges[2*e+1] = j;
			weights[e] = adjacentMatrix[oddDegree[i]][oddDegree[j]];
 			e++;
		}
	}
	if (e != edge_num) {
		cout<<"the number of edge is wrong"<<endl;
		exit(1);
	}
}

void PrintMatching(int node_num, PerfectMatching* pm, int * oddDegree, float ** adjacentMatrix) {
	int i, j;

	for (i=0; i<node_num; i++) {
		j = pm->GetMatch(i);
		if (i < j) printf("%d %d %f\n", oddDegree[i], oddDegree[j], adjacentMatrix[oddDegree[i]][oddDegree[j]] );
	}
}

void calculateDegree(int * degree, MST* mst, int * node_num, int N){
/*	for (int j = 0; j < N; j++){
		if (mst->parent[j] > -1)
			degree[mst->parent[j]]++;
	}
	for (int i=0; i<N; i++){
		if (i!=0)
			degree[i]++;
		if (degree[i] % 2 != 0)
			(*node_num)++;
	}
*/
	for (int i = 0; i < N; i++){
		Node * node = mst->nodes.find(i)->second;
		degree[i] = node->getEdges().size();
		if (degree[i] % 2 != 0){
			(*node_num)++;
		}
	}


}


int main() {
	cerr << "Program start.";
	set< pair<int,int> > generatedPointset;
	float** adjacentMatrix;
	int W, H, N;
	Point pointset;

	W = 1000;
	H = 1000;
	N = 1000;

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal

	MST * mst2 = new MST(adjacentMatrix, N);
	mst2->makeTree();
	//mst2->printMST();
	cout<< "Just the MST:" << endl;
	cout << "Mean: " << mst2->calMean(MST_1) << endl;
	cout << "Standard Deviation: " << mst2->calStd(MST_1)<<endl;
	cout << "MST total weight: " << mst2->mst_total_weight << endl;

	cout << "-------------------" << endl;

	//Deliverable B: Find TSP2 path from the constructed MST

	clock_t begin = clock();
	mst2->makeTSP2();
	clock_t end = clock();
	double elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
	cout << "TSP2 results:" << endl;
	cout << "Mean: "<< mst2->calMean(TSP2) << endl;
	cout << "Standard Deviation: " << mst2->calStd(TSP2) << endl;
	cout << "Total Weight: " << mst2->tsp2_total_weight << endl;



	//Deliverable C: Find TSP1.5 path from the constructed MST
	//construct MST based off of arrays

	for (int i = 1; i < N; i++){
		mst2->addEdge(mst2->parent[i], i, mst2->key[i]);
	}

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
	cerr <<"Node num: "<<  node_num << " " << *nn_ptr<< endl;
	int* edges;
	int* weights;
	PerfectMatching *pm = new PerfectMatching(node_num, edge_num);
	LoadInput(node_num, edge_num, edges, weights, adjacentMatrix, N, oddDegree);
	for (e=0; e<edge_num; e++) {
		pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
	}

	pm->options = options;
	pm->Solve();

	float cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
	printf("Total cost of the perfect min-weight matching = %.1f\n", cost);

	//PrintMatching(node_num, pm, oddDegree ,adjacentMatrix);
	//Code to add Perfect matching edges to MST
	int count = 0;
	for (int i=0; i<node_num; i++) {
		j = pm->GetMatch(i);
		if (i < j){
			mst2->addEdge(oddDegree[i], oddDegree[j], adjacentMatrix[oddDegree[i]][oddDegree[j]]);
			count++; 
		}
	}
	int total_edges = (N-1) + count;
	cerr << "TotalEdges: " << total_edges  << endl;
	mst2->makeTSP1_5(total_edges);


	delete pm;
	delete [] oddDegree;
	delete [] degree;
	delete [] edges;
	delete [] weights;
	
	
	return 0;
}

