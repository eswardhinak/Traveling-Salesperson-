#include "MST.h"
#include <math.h>
MST::MST(int** input, int size) {
	adjacentMatrix = input;
	key = new int[size];   
    mstSet = new bool[size];  
	parent = new int[size];
	visited = new bool[size];
	tsp2 = new int[size];
	N = size;
}

MST::~MST() {

}
void MST::makeTree() { 
     // Initialize all keys as INFINITE
     for (int i = 0; i < N; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 	
     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST 
 
     // The MST will have V vertices
     for (int count = 0; count < N-1; count++)
     {
        //cerr << count << endl;
        // Pick thd minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet);
 
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
 
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < N; v++)
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if adjacentMatrix[u][v] is smaller than key[v]
          if (adjacentMatrix[u][v] && mstSet[v] == false && adjacentMatrix[u][v] <  key[v]){
             parent[v]  = u, key[v] = adjacentMatrix[u][v];
          }
     }
     /**
     cout << "key array:" << endl;
     for (int h=0; h < N; h++){
     	cout << key[h] << ' ';
     }
     cout << endl << "mstSet array" << endl;
     for (int j=0; j < N; j++){
     	cout  << mstSet[j] << ' ';
     }
     cout << endl << "parent array" << endl;
     for (int l = 0; l < N ; l++){
     	cout << parent[l] << " ";
     }**/
}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int MST::minKey(int key[], bool mstSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;
 
   for (int v = 0; v < N; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;
   //cerr << "Min: " << min << "  Min_Index: " << min_index << endl;
   return min_index;
}

// A utility function to print the constructed MST stored in parent[]
void MST::printMST() {
	cout<<endl;
	cout<<"Minimum spanning tree from the adjacency matrix"<<endl;
	cout<<"Edge   Weight"<<endl;
	for (int i = 1; i < N; i++) {
		cout<<parent[i]<<" - "<<i<<"  "<<adjacentMatrix[i][parent[i]]<<endl;
	}
}

//calculate mean of all edges in the MST
float MST::calMean(int option) {
	float mean = 0.0;
 	if(option == MST_1) {
		for (int j=1; j < N; j++){
			mean += key[j];
		}
		mean /= (N-1);
	}else if(option == TSP2) {
		for (int j=0; j < N; j++){
			mean += tsp2[j];
		}
		mean /= N;
	} else if(option == TSP1_5) {

	}

	return mean;
}

//calculate standard deviation of all edges in the MST
float MST::calStd(int option) {
	float std = 0.0;

	if(option == MST_1) {
		float mean = calMean(option);
		float dist_from_mean = 0.0;
		float mean_dist = 0.0;
		for (int j=1; j < N; j++){
			float dist = key[j] - mean;
			dist_from_mean= dist * dist;
			mean_dist += dist_from_mean;
		}
		mean_dist /= (N-1);
		std = sqrt(mean_dist);


	}else if(option == TSP2) {
		float mean = calMean(option);
		float dist_from_mean = 0.0;
		float mean_dist = 0.0;
		for (int j = 0; j < N; j++){
			float dist = key[j]  - mean;
			dist_from_mean = dist*dist;
			mean_dist += dist_from_mean;
		}
		mean_dist /= N;
		std = sqrt(mean_dist);
	} else if(option == TSP1_5) {

	}

	return std;
}

void MST::makeTSP2() {
	//make a Eulerian tour by DFS
	//http://www.geeksforgeeks.org/travelling-salesman-problem-set-2-approximate-using-mst/
	cerr << "Starting 2-approximation TSP" << endl;
	int source=0;
	int prev_vertex =source;
	int count = 0;
	int curr_vertex = source;
	int tsp_cost = 0;
	int new_edge=0;
	//make a tour using DFS
	while (count < N-1){
		curr_vertex = findValue(curr_vertex);
		//cerr << prev_vertex << "--" << curr_vertex << ' ';
		if (curr_vertex == -1)
			return;
		new_edge = adjacentMatrix[curr_vertex][prev_vertex];
		tsp_cost += new_edge;
		//cerr << adjacentMatrix[curr_vertex][prev_vertex] << endl;
		prev_vertex = curr_vertex;
		tsp2[count] = new_edge;
		count++;
	}
	new_edge = adjacentMatrix[curr_vertex][source];
	tsp_cost += new_edge;
	tsp2[count] = new_edge;
	cerr << tsp_cost << endl;

}
int MST::findValue(int padre){
	if (padre == -1){
		return -1;
	}
	for (int i = 1; i < N; i++){
		if (parent[i] == padre && !visited[i]){
			visited[i] = true;
			return i;
		}
	}
	return findValue(parent[padre]);
}

void MST::makeTSP1_5() {
	
	//construct minimum-weight-matching for the given MST
	minimumMatching();

	//make all edges has even degree by combining mimimum-weight matching and MST
	combine();

	//calculate heuristic TSP cost 
}

void MST::minimumMatching() { //if you choose O(n^2)
	//add edges to MST using minimum matching code provided

}

void MST::combine() {
	//combine minimum-weight matching with the MST to get a multigraph which has vertices with even degree
	
}