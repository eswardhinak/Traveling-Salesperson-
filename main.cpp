#include "common.h"
#include "Point.h"
#include <ctime>
#include "MST.h"


int main() {
	cerr << "Program start.";
	set< pair<int,int> > generatedPointset;
	int** adjacentMatrix;
	int W, H, N;
	Point pointset;

	W = 20000;
	H = 20000;
	N = 10000;

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation
	//pointset.printPointset();

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
		//create a MST creator class that will return the MST
		//construct an MST
	MST * mst2 = new MST(adjacentMatrix, N);
	mst2->makeTree();
	//mst2->printMST();
	cout << "Mean: " << mst2->calMean(MST_1) << endl;
	cout << "Standard Deviation: " << mst2->calStd(MST_1)<<endl;
	clock_t begin = clock();
	mst2->makeTSP2();
	clock_t end = clock();
	double elapsed_seconds = double(end - begin) / CLOCKS_PER_SEC;
	cerr << "TSP_time:  " << elapsed_seconds << endl;
	cout << "Mean: "<< mst2->calMean(TSP2) << endl;
	cout << "Standard Deviation: " << mst2->calStd(TSP2) << endl;



	//Deliverable B: Find TSP2 path from the constructed MST

	//Deliverable C: Find TSP1.5 path from the constructed MST

	return 0;
}

