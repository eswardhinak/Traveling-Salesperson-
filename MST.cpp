#include "MST.h"
#include <math.h>
#include "Edge.hpp"
#include "Node.hpp"
#include <iostream>
#include <array>

MST::MST(float** input, int size) {
	adjacentMatrix = input;
	key = new float[size];   
	memset(key, 0, size);
    mstSet = new bool[size]; 
    memset(mstSet, false, size); 
	parent = new int[size];
	memset(parent, -1, size);
	visited = new bool[size];
	memset(visited,false,size);
	tsp2 = new float[size];
	memset(tsp2, 0, size);
	N = size;
	tsp2_total_weight = 0;
	mst_total_weight = 0;
	tsp1_5_total_weight = 0;
}

MST::~MST() {
	//go through and delete all the dynamically allocated arrays. 
	delete [] key;
	delete [] mstSet;
	delete [] parent;
	delete [] visited;
	delete [] tsp2;
	auto it = nodes.begin();
	for (; it!= nodes.end(); it++){
		delete(it->second); //delete node
		it->second = NULL; //set node pointer to null
	}
	for (int i = 0; i < N; i++){
		free (adjacentMatrix[i]);
	}
	free(adjacentMatrix);
	adjacentMatrix = nullptr;

}

/*Method to add edge to the MST*/

void MST::addEdge(int from, int to, unsigned int cost){
	Node * node1;
	Node * node2;

	std::unordered_map<int, Node *>::const_iterator first_node = nodes.find(from); //first vertex
	//if no vertex with that name was found
	if (first_node  == nodes.end()){
		node1 = new Node(from);
		std::pair<int, Node*> first_node(from, node1);
		nodes.insert(first_node);
	}
	else{
		node1 = first_node->second;
	}
	std::unordered_map<int, Node*>::const_iterator second_node = nodes.find(to);
	//no node with that name was found
	if (second_node == nodes.end()){
		node2 = new Node(to); // create new node
		std::pair<int, Node*> second_node (to, node2);
		nodes.insert(second_node);
	}
	else{
		node2 = second_node->second; //find vertex inside our graph
	}

	node1->addEdge(node2, cost);
	node2->addEdge(node1, cost);
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
     for (int j = 1; j < N; j++){
     	mst_total_weight += key[j];
     }
}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int MST::minKey(float key[], bool mstSet[])
{
   // Initialize min value
   float min = std::numeric_limits<float>::max();
   int min_index;
 
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
		mst_total_weight = mean;
		mean /= (N-1);
	}else if(option == TSP2) {
		for (int j=0; j < N; j++){
			mean += tsp2[j];
		}
		tsp2_total_weight = mean;
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
		float mean = calMean(option);
		//float dist_from_mean = 
	}

	return std;
}

void MST::makeTSP2() {
	//make a Eulerian tour by DFS
	int source=0; //"root" vertex to start DFS from. 
	int prev_vertex =source;
	int count = 0;
	int curr_vertex = source;
	int tsp_cost = 0;
	float new_edge=0;
	//make a tour using DFS
	while (count < N-1){
		curr_vertex = findValue(curr_vertex);
		if (curr_vertex == -1)
			return;
		new_edge = adjacentMatrix[curr_vertex][prev_vertex];
		tsp_cost += new_edge;
		prev_vertex = curr_vertex;
		tsp2[count] = new_edge;
		count++;
	}
	new_edge = adjacentMatrix[curr_vertex][source];
	tsp_cost += new_edge;
	tsp2[count] = new_edge;
	for (int j = 0; j < N; j++){
		tsp2_total_weight += tsp2[j];
	}

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

//counts the number of vertices reachable in the present graph 
int MST::DFSCount(int v, bool * encountered){
	encountered[v] = true;
	int count = 1;

	Node * nodeOne = (nodes.find(v))->second;
	auto it = nodeOne->edges.begin();
	for (; it != nodeOne->edges.end(); it++){
		Edge *  edgeOne = *it;
		int node_id = (edgeOne->getTo())->getName();
		if (!edgeOne->getTraversed() && !encountered[node_id])
			count += DFSCount(node_id, encountered);
	}
	return count;
}
int MST::isValidEdge(int u, int v){
	//cerr << "Is Valid Edge: " << u << " " << v << endl;
	int count = 0; // stores count of adjacent vertices
	Node * nodeOne = (nodes.find(u))->second;
	auto it = (nodeOne->edges).begin();
	for (; it != nodeOne->edges.end(); it++){
		Edge * edgeOne = *it;
		if (!edgeOne->getTraversed()){
			count++;
		}
	} 
	//only path that we can take 
	if (count == 1) return true;

	bool * encountered = new bool[total_edges]; //to be used in DFSCount
	memset(encountered, false, total_edges); //reset locations to false
	int count1 = DFSCount(u, encountered); //count using the edge

	hideEdge(u, v); //remove the edge from the graph
	memset(encountered, false, total_edges); //reset locations to false
	int count2 = DFSCount(u, encountered); //recount without the edge

	unHideEdge(u, v); //add the edge back to the graph
	//if count1 is greater than count2 then the edge is a bridge.
	delete [] encountered;
	return (count1 > count2) ? false:true; 
}
/*Unhide both instances of an edge with vertices u and v*/
void MST::unHideEdge(int u, int v){
	//get the Node pointers
	Node * nodeOne = (nodes.find(u))->second; 
	Node * nodeTwo = (nodes.find(v))->second;
	//get the Edge pointers
	Edge * edgeOne;
	auto it1 = (nodeOne->edges).begin();
	for (; it1 != (nodeOne->edges).end(); it1++){
		Edge * currEdge = *it1;
		if (currEdge->getTo()->getName() == v && currEdge->getTraversed() == true){
			edgeOne = currEdge;
			break;
		}
	}
	Edge * edgeTwo;
	auto it2 = (nodeTwo->edges).begin();
	for (; it2 != (nodeTwo->edges).end(); it2++){
		Edge * currEdge = *it2;
		if (currEdge->getTo()->getName() == u && currEdge->getTraversed() == true){
			edgeTwo = currEdge;
			break;
		}
	}
	//set traversal flag to false
	edgeOne->setTraversed(false);
	edgeTwo->setTraversed(false);
}

/*Hide both instances of an edge with vertices u and v*/
void MST::hideEdge(int u, int v){
	//get Node pointers
	Node * nodeOne = (nodes.find(u))->second;
	Node * nodeTwo = (nodes.find(v))->second;
	//get Edge pointers
	Edge * edgeOne;
	auto it1 = (nodeOne->edges).begin();
	for (; it1 != (nodeOne->edges).end(); it1++){
		Edge * currEdge = *it1;
		if (currEdge->getTo()->getName() == v && currEdge->getTraversed() == false){
			edgeOne = currEdge;
			break;
		}
	}
	Edge * edgeTwo;
	auto it2 = (nodeTwo->edges).begin();
	for (; it2 != (nodeTwo->edges).end(); it2++){
		Edge * currEdge = *it2;
		if (currEdge->getTo()->getName() == u && currEdge->getTraversed() == false){
			edgeTwo = currEdge;
			break;
		}
	}
	//set traversal flag to true
	edgeOne->setTraversed(true);
	edgeTwo->setTraversed(true);	
}

//method that starts the 1.5 approximation of TSP
void MST::makeTSP1_5(int total_edges) {
	//total edges in MST + minimum matching
	this->total_edges = total_edges;
 	int * traversal = new int[total_edges];
	int ci = 0;
	int * ci_ptr = &ci;
	memset(traversal, 0, total_edges);
	traversal[0] = 0;
	eulerTraversal(0, ci_ptr, traversal, total_edges);
	bool * duplicateSearch = new bool[N];
	noDuplicates = new int[N]; //for tsp+shortcut(1.5)
	memset(duplicateSearch, false, N);
	int j=0;
	for (int i = 0; i < total_edges; i++){
		int curr_vertex = traversal[i];
		if (duplicateSearch[curr_vertex] == false){
			noDuplicates[j] = curr_vertex;
			duplicateSearch[curr_vertex] = true;
			j++;
		}	
	}
	tsp1_5_total_weight = 0.0;
	int k;
	for (k = 0; k < (N-1); k++){
		tsp1_5_total_weight += adjacentMatrix[noDuplicates[k]][noDuplicates[k+1]];
		//cerr << noDuplicates[k] << " " << noDuplicates[k+1] << endl;
	}
	tsp1_5_total_weight += adjacentMatrix[noDuplicates[k]][noDuplicates[0]];
	delete [] traversal;
	delete [] duplicateSearch;
	//delete [] noDuplicates;
}

void MST::eulerTraversal(int current, int * ci,  int * traversal, int total_edges){
	Node * nodeOne = (nodes.find(current))->second;
	//cerr << "Current Node: " << nodeOne->getName() << endl;
	if (*ci > -1 && *ci < total_edges)
		traversal[*ci] = current;
	auto it = (nodeOne->edges).begin();
	for (; it != (nodeOne->edges).end(); it++){
		Edge * edgeOne = *it;
		int toNode = edgeOne->getTo()->getName();
		if (!edgeOne->getTraversed() && isValidEdge(current, toNode)){
			hideEdge(current,toNode);
			(*ci)++;
			eulerTraversal(toNode, ci , traversal, total_edges);
		}
	}
}
void MST::vertexRemoval(int numEdges){
	int * traversalNew;
	traversalNew = new int[N];
	//cerr << "Original:";
	for (int j = 0; j < N; j++){
		//cerr << noDuplicates[j];
		traversalNew[j] = noDuplicates[j];
	}
	//cerr << endl;
	vr_total_weight = 0.0;
	float minChange = 0;
	int finalI, finalJ;
	do {
		minChange = 0;
		for (int i = 1; i < N-2; i++){
			for (int j = i+1; j < N-1; j++){
				int vert1 = traversalNew[i];
				int vert1Prev = traversalNew[i-1];
				int vert1Next = traversalNew[i+1];
				int vert2 = traversalNew[j];
				int vert2Prev = traversalNew[j-1];
				int vert2Next = traversalNew[j+1];
				float totalPrev = 0.0;
				float newTotal = 0.0;
				if (j - i  == 1){
					 totalPrev = adjacentMatrix[vert1Prev][vert1] + adjacentMatrix[vert2][vert2Next] + adjacentMatrix[vert1][vert2];
			 		 newTotal = adjacentMatrix[vert1Prev][vert2] + adjacentMatrix[vert1][vert2Next]  + adjacentMatrix[vert1][vert2];
				}
				else{
					float prev_edge1_combo = adjacentMatrix[vert1Prev][vert1] + adjacentMatrix[vert1][vert1Next];
					float prev_edge2_combo = adjacentMatrix[vert2Prev][vert2] + adjacentMatrix[vert2][vert2Next];
					totalPrev = prev_edge2_combo + prev_edge1_combo;
					newTotal = adjacentMatrix[vert1Prev][vert2]+ adjacentMatrix[vert2][vert1Next] + adjacentMatrix[vert2Prev][vert1] + adjacentMatrix[vert1][vert2Next];
				}
				//cerr << "Previous : " << totalPrev << "Potential: " << newTotal << endl;
				float change = newTotal  - totalPrev;
				if (minChange > change){
					//	cerr << "First Vertex Combo:  " << prev_edge1_combo;

					//cerr << "   Second vertex combo: " << prev_edge2_combo << endl;

					//cerr << change << endl;
					minChange = change;
					finalI = i;
					finalJ = j;
				}

			}
		}
		if (minChange < 0){
		//	cerr << finalI << " " << finalJ << endl;
			int swap = traversalNew[finalJ];
			traversalNew[finalJ] = traversalNew[finalI];
			traversalNew[finalI] = swap;
		}
		//cerr << "now: " <<endl;
		for (int i = 0; i < N; i++){
			//cerr << traversalNew[i];
		}
		//cerr << endl;
	}while(minChange < 0);
	for (int i = 0; i < N; i++){
		vr_total_weight += adjacentMatrix[traversalNew[i]][traversalNew[(i+1)%N]];
	}
	delete [] traversalNew;
}
void MST::edgeRemoval(int numEdges){
	int * traversalNew;
	traversalNew = new int[N];
	//cerr << "OriginalL :";
	for (int k = 0; k < N; k++){
		//cerr << noDuplicates[k] << "  ";
		traversalNew[k] = noDuplicates[k];
	}
	er_total_weight = 0.0;
	int tc = 0;
	float minChange = 0;
	do {
		minChange = 0;
		int finalI, finalJ; //best indices for swap 
		for (int i = 0; i < N-2; i++){
			for (int j = i+2; j < N; j++){
				int vert1 = traversalNew[i];
				int vert2 = traversalNew[(i+1)%N];
				int vert3 = traversalNew[j];
				int vert4 = traversalNew[(j+1)%N];
				float edge1 = adjacentMatrix[vert1][vert2];
				float edge2 = adjacentMatrix[vert3][vert4];
				float total_edgecost = edge1 + edge2;
				float edge1New = adjacentMatrix[vert1][vert3];
				float edge2New = adjacentMatrix[vert2][vert4];
				float newTotal = edge1New + edge2New;
				float change = newTotal - total_edgecost;
				if (minChange > change){
					minChange = change;
					finalI = i;
					finalJ = j;
				}
			}
		}
		//do the best swap 
		if (minChange < 0){
			int midpoint = finalI+1+(finalJ-(finalI+1))/2;
			int count = 0;
			//cerr << "Replacing: " << traversalNew[finalI] << "-->" << traversalNew[(finalI+1)%N] << " and " << traversalNew[finalJ] << "-->" << traversalNew[(finalJ+1)%N] << endl;
			for (int k = finalI+1; k <= midpoint; k++){
				int swap = traversalNew[k];
				traversalNew[k] = traversalNew[finalJ - count];
				traversalNew[finalJ-count] = swap;
				count++;
			}
			/*cerr << "Now:" <<endl;
			for (int s = 0; s < N; s++){
				cerr << traversalNew[s];
			}
			cerr << endl;*/
		}
	}while(minChange < 0);
	for (int i = 0; i < N; i++){
		er_total_weight += adjacentMatrix[traversalNew[i]][traversalNew[(i+1)%N]];
	}
	delete [] traversalNew;
}
