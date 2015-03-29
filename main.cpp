#include "common.h"
#include "Point.h"
#include <ctime>
#include <math.h>
#include "MST.h"
#include "MinMatching/PerfectMatching.h"

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

double calSTD(double mean, double * values, int trials, bool debug){
	double distance = 0;
	double meanDist = 0;
	for (int i =0 ; i < trials; i++){
		if (debug)
			cerr << values[i]<< endl;
		double dist = values[i] - mean;
		distance = dist*dist;
		//if (debug)
			//cerr << distance;
		meanDist += distance;
	}
	if (debug)
		cerr <<meanDist;
	meanDist /= trials;
	if (debug)
	cerr<<meanDist;
	double std = sqrt(meanDist);
	return std;
}
int main(int argc, char ** argv) {

	if (argc != 6){
		cerr << "Please enter the five command line arguments in this order:\n\tW - width\n\tH - height\n\tN - number of vertices\n\tT - number of trials \n\tE - extra credit 20PT-E and 20PT-V \n\t\t(1 to run EC and 0 for don't run it)\n";
		cerr << "\n\nNote: If number of vertices is very high (near 10000), EC1 will take\n\t exceptionally long time. Expect 30 minutes for 20PT-E \n\tand 20PT-V each per trial" << endl;
		return 0;
	}
	set< pair<int,int> > generatedPointset;
	float** adjacentMatrix;
	int W, H, N;

	//convert command line arguments to W, H, N
	W = atoi(argv[1]); 
	H = atoi(argv[2]);
	N = atoi(argv[3]);
	int trials = atoi(argv[4]);
	int flag = atoi(argv[5]);
	double * mstTrials = new double[trials];
	memset(mstTrials, 0, trials);
	double * timeMSTTrials = new double[trials];
	memset(timeMSTTrials, 0, trials);
	double * tsp2Trials = new double[trials];
	memset(tsp2Trials, 0, trials);
	double * timeTSP2Trials = new double[trials];
	memset(timeTSP2Trials, 0, trials);
	double * tsp1_5Trials = new double[trials];
	memset(tsp1_5Trials, 0, trials);
	double * timeTSP1_5Trials  = new double[trials];
	memset(timeTSP1_5Trials, 0, trials);
	double * ECEdges = new double[trials];
	memset(ECEdges, 0, trials);
	double * ECVertices = new double[trials];
	memset(ECVertices, 0, trials);
	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl; 
	for (int z = 0 ; z < trials; z++){
		float OptimumEdgeRemoval = 0;
		float OptimumVertexRemoval = 0;
		Point pointset;
		pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation

		generatedPointset = pointset.getPointset();
		adjacentMatrix = pointset.getAdjacentMatrix();

		//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
		clock_t begin = clock();
		MST * mst2 = new MST(adjacentMatrix, N);
		mst2->makeTree();
		clock_t end = clock();
		double elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
		//cout << "-------------------" << endl;
		//cout<< "MST:" << endl;
		//cout << "MST total weight: " << mst2->mst_total_weight << endl;
		//cout << "Time taken: "<< elapsed_seconds << endl;
		//cout << "-------------------" << endl;
		mstTrials[z] = mst2->mst_total_weight;
		timeMSTTrials[z] = elapsed_seconds;

		//Deliverable B: Find TSP2 path from the constructed MST
		begin = clock();
		mst2->makeTSP2();
		end = clock();
		elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
		//cout << "TSP-2 Approximation:" << endl;
		//cout << "Total Weight: " << mst2->tsp2_total_weight << endl;
		//cout << "Time taken: " << elapsed_seconds << endl;
		//cout << "------------------" << endl;
		tsp2Trials[z] = mst2->tsp2_total_weight;
		timeTSP2Trials[z] = elapsed_seconds;

		//Deliverable C: Find TSP1.5 path from the constructed MST
		//construct MST based off of arrays
		
		begin = clock();
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
		//cerr << "TotalEdges: " << total_edges  << endl;
		mst2->makeTSP1_5(total_edges);
		end = clock();
		elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
		//cout << "TSP-1.5 Approximation:" << endl;
		printf("Total cost of the perfect min-weight matching = %.1f\n", cost);
		//cout << "Total Weight: " << mst2->tsp1_5_total_weight << endl;
		//cout << "Time taken: " << elapsed_seconds << endl;
		//cout << "------------------" << endl;
		tsp1_5Trials[z] = mst2->tsp1_5_total_weight;
		timeTSP1_5Trials[z] = elapsed_seconds;
		if (flag == 1){
			cout << "20PT-E EC ...." <<endl;
			mst2->edgeRemoval(total_edges);
			OptimumEdgeRemoval = mst2->er_total_weight;
			cout << "20PT_V EC ... " <<endl;
			mst2->vertexRemoval(total_edges);
			OptimumVertexRemoval = mst2->vr_total_weight;
		}
		ECEdges[z] = OptimumEdgeRemoval;
		ECVertices[z] = OptimumVertexRemoval;
		delete pm;
		delete [] oddDegree;
		delete [] degree;
		delete [] edges;
		delete [] weights;
		delete mst2;
	}

	//Mean calculation
	cout << "\n\n\nWidth: " << W << "   Height: " << H << "   Vertices: " << N << "    Trials: " << trials << endl;
	cout << "___________________________" << endl;
	double meanMST = 0.0;
	double meanMSTTime = 0.0;
	double meanTSP2 = 0.0;
	double meanTSP2Time = 0.0;
	double meanTSP1_5 = 0.0;
	double meanTSP1_5Time = 0.0;
	double meanECEdges = 0.0;
	double meanECVertices = 0.0;
	for (int i = 0; i < trials; i++){
		meanMST += mstTrials[i];
		meanMSTTime += timeMSTTrials[i];
		meanTSP2 += tsp2Trials[i];
		meanTSP2Time += timeTSP2Trials[i];
		meanTSP1_5 += tsp1_5Trials[i];
		meanTSP1_5Time += timeTSP1_5Trials[i]; 
		meanECEdges += ECEdges[i];
		meanECVertices += ECVertices[i];
	}
	meanMST /= trials;
	meanMSTTime /= trials;
	meanTSP2 /= trials;
	meanTSP2Time /= trials;
	meanTSP1_5 /= trials;
	meanTSP1_5Time /= trials;
	meanECEdges /= trials;
	meanECVertices /= trials;

	double stdMST = calSTD(meanMST, mstTrials, trials, false);
	double stdMSTTime = calSTD(meanMSTTime, timeMSTTrials, trials, false);
	double stdTSP2 = calSTD(meanTSP2, tsp2Trials, trials, false);
	double stdTSP2Time = calSTD(meanTSP2Time, timeTSP2Trials, trials, false);
	double stdTSP1_5 = calSTD(meanTSP1_5, tsp1_5Trials, trials, false);
	double stdTSP1_5Time = calSTD(meanTSP1_5Time, timeTSP1_5Trials, trials, false);
	double stdECEEdges = calSTD(meanECEdges, ECEdges, trials, false);
	double stdECVertices = calSTD(meanECVertices, ECVertices, trials, false);

	delete [] mstTrials;
	delete [] timeMSTTrials;	
	delete [] tsp2Trials;
	delete [] timeTSP2Trials;
	delete [] tsp1_5Trials;
	delete [] timeTSP1_5Trials;
	delete [] ECEdges;
	delete [] ECVertices;
	
	cout << "Costs" << endl;
	cout << "MST --> Mean: " << meanMST << "  Std: " << stdMST <<endl;
	cout << "TSP-2 --> Mean: " << meanTSP2 << "  Std: " <<stdTSP2 <<endl;
	cout << "TSP-1.5 --> Mean: " <<meanTSP1_5 << " Std: " <<stdTSP1_5 << endl;
	cout << "ECE 20PT-E --> Mean: " << meanECEdges << " Std: " <<stdECEEdges << endl;
	cout << "ECE 20PT-V ---> Mean: " << meanECVertices << " Std: " <<stdECVertices <<endl;

	/*cout << "Time To Calculate" << endl;
	cout << "MST --> Mean: " << meanMSTTime << "  Std: " << stdMSTTime << endl;
	cout << "TSP-2 -->  Mean: " <<meanTSP2Time << "  Std: " <<stdTSP2Time <<endl;
	cout << "TSP-1.5--> Mean: " <<meanTSP1_5Time << " Std: " <<stdTSP1_5Time <<endl;*/
	return 0;
}

