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
	edges = new int[2*edge_num]; //hold the edges
	weights = new int[edge_num]; //weights of the edges
	//create the list of edges
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

//Print the Minimum matching 
void PrintMatching(int node_num, PerfectMatching* pm, int * oddDegree, float ** adjacentMatrix) {
	int i, j;

	for (i=0; i<node_num; i++) {
		j = pm->GetMatch(i);
		if (i < j) printf("%d %d %f\n", oddDegree[i], oddDegree[j], adjacentMatrix[oddDegree[i]][oddDegree[j]] );
	}
}
//calculate the degree. node_num is a pointer to the number of vertices with odd
//degrees. 
void calculateDegree(int * degree, MST* mst, int * node_num, int N){

	for (int i = 0; i < N; i++){
		Node * node = mst->nodes.find(i)->second;
		degree[i] = node->getEdges().size();
		if (degree[i] % 2 != 0){
			(*node_num)++;
		}
	}

}


int main(int argc, char ** argv) {
	set< pair<int,int> > generatedPointset;
	float** adjacentMatrix;
	int W, H, N;
	Point pointset;

	//convert command line arguments to W, H, N
	W = atoi(argv[1]); 
	H = atoi(argv[2]);
	N = atoi(argv[3]);

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl; 

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
	MST * mst2 = new MST(adjacentMatrix, N);
	mst2->makeTree();
	cout << "-------------------" << endl;
	cout<< "MST:" << endl;
	cout << "Mean: " << mst2->calMean(MST_1) << endl;
	cout << "Standard Deviation: " << mst2->calStd(MST_1)<<endl;
	cout << "MST total weight: " << mst2->mst_total_weight << endl;

	cout << "-------------------" << endl;

	//Deliverable B: Find TSP2 path from the constructed MST
	mst2->makeTSP2();
	double elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
	cout << "TSP-2 Approximation:" << endl;
	cout << "Mean: "<< mst2->calMean(TSP2) << endl;
	cout << "Standard Deviation: " << mst2->calStd(TSP2) << endl;
	cout << "Total Weight: " << mst2->tsp2_total_weight << endl;
	cout << "------------------" << endl;

	//Deliverable C: Find TSP1.5 path from the constructed MST
	//construct MST based off of arrays
	for (int i = 1; i < N; i++){
		mst2->addEdge(mst2->parent[i], i, mst2->key[i]); //add edges from existing MST representation
	}

	struct PerfectMatching::Options options; //algorithm for PerfectMatching
	int i, e, node_num = 0; 
	int * nn_ptr = &node_num;
	int * degree =  new int[N];

	calculateDegree(degree, mst2, nn_ptr, N); //get odd degrees 
	int * oddDegree = new int[node_num];

	int j = 0; //populate oddDegree for use by PerfectMatching code
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
	//add odd degree vertices to the Perfect Matching code.
	for (e=0; e<edge_num; e++) {
		pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
	}

	pm->options = options;
	pm->Solve(); //solve perfect matching

	float cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
	cerr << "-------------------" << endl;
	printf("Total cost of the perfect min-weight matching = %.1f\n", cost);
	cerr << "-------------------" << endl;

	int count = 0;
	//add all the edges from perfect matching to the graph representation
	for (int i=0; i<node_num; i++) {
		j = pm->GetMatch(i);
		if (i < j){
			mst2->addEdge(oddDegree[i], oddDegree[j], adjacentMatrix[oddDegree[i]][oddDegree[j]]);
			count++; 
		}
	}

	int total_edges = (N-1) + count; //total edges is number of edges in MST plus min matching
	cerr << "TotalEdges: " << total_edges  << endl;
	mst2->makeTSP1_5(total_edges);

	cout << "TSP-1.5 Approximation:" << endl;
	cout << "Mean: "<< mst2->calMean(TSP1_5) << endl;
	cout << "Standard Deviation: " << mst2->calStd(TSP1_5) << endl;
	cout << "Total Weight: " << mst2->tsp2_total_weight << endl;
	cout << "------------------" << endl;
	delete pm;
	delete [] oddDegree;
	delete [] degree;
	delete [] edges;
	delete [] weights;
	
	
	return 0;
}

