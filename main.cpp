#include "common.h"
#include "Point.h"

int ** calculateMST(int ** adjacentMatrix, int size, int max){
	//Will use Prim's Algorithm to calculate MST
	int source = 0;
	int current = source;
	int min_source;
	int * visited;
	visited = (int*)calloc(size, sizeof(int));
	for (int k =0 ; k<size; k++){
		visited[k]=0;
	}
	visited[source] = 1;
	int ** minST;
	minST = (int**)calloc(size, sizeof(int *));

	for(int i=0; i<size ; ++i) {
		minST[i] = (int*)calloc(size, sizeof(int));
	}
	cerr << "Passed array allocation" << endl;
	for (int j = 0; j < (size-1); j++){
		cerr << "current: " << current << endl;
		int min_edge=0;
		int min_len = max;
		int i = 0;
		for (i=0; i < size; i++){
			int curr_len = adjacentMatrix[current][i];
			if (curr_len < min_len && !visited[i]){
				min_edge = i;
				min_len = curr_len;
			}
		}
		minST[current][min_edge] = min_len;
		minST[min_edge][current] = min_len;
		current = min_edge;
		visited[current] = 1;
		
	}
	return minST;

}

int main() {
	set< pair<int,int> > generatedPointset;
	int** adjacentMatrix;
	int W, H, N;
	Point pointset;

	W = 10000;
	H = 10000;
	N = 1000;

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation
	pointset.printPointset();

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
		//create a MST creator class that will return the MST
		//construct an MST
	int max;
	if (W>=H) max = W; 
	else max = H;
	int ** min_span_tree = calculateMST(adjacentMatrix, N, max+1);
	
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			cout << min_span_tree[i][j] << "\t ";
		}
		cerr << endl;
	}

	//Deliverable B: Find TSP2 path from the constructed MST

	//Deliverable C: Find TSP1.5 path from the constructed MST

	return 0;
}

