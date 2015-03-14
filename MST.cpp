#include "MST.h"
#include <math.h>
#include "Edge.hpp"
#include "Node.hpp"
#include <iostream>

MST::MST(float** input, int size) {
	adjacentMatrix = input;
	key = new int[size];   
    mstSet = new bool[size];  
	parent = new int[size];
	visited = new bool[size];
	tsp2 = new int[size];
	N = size;
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
	//cerr << tsp_cost << endl;

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
}/*
int MST::recursiveFindNextAvailableNode(bool * visited, int current){
	auto it = nodes.find(current);
	Node * currentNode = it->second;

}*/
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
	//
	int * traversal = new int[total_edges];
	int ci = 0;
	int * ci_ptr = &ci;
	memset(traversal, 0, total_edges);
	traversal[0] = 0;
	eulerTraversal(0, ci_ptr, traversal);
	bool * duplicateSearch = new bool[N];
	int * noDuplicates = new int[N];
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
	tsp1_5_total_weight = 0;
	for (int k = 0; k < (N-1); k++){
		tsp1_5_total_weight += adjacentMatrix[noDuplicates[k]][noDuplicates[k+1]];
	}
}

void MST::eulerTraversal(int current, int * ci,  int * traversal){
	Node * nodeOne = (nodes.find(current))->second;
	//cerr << "Current Node: " << nodeOne->getName() << endl;
	traversal[*ci] = current;
	auto it = (nodeOne->edges).begin();
	for (; it != (nodeOne->edges).end(); it++){
		Edge * edgeOne = *it;
		int toNode = edgeOne->getTo()->getName();
		//cerr << toNode << endl;
		if (!edgeOne->getTraversed() && isValidEdge(current, toNode)){
			hideEdge(current,toNode);
			(*ci)++;
			eulerTraversal(toNode, ci , traversal);
		}
	}
}


void MST::minimumMatching() { //if you choose O(n^2)
	//add edges to MST using minimum matching code provided

}

void MST::combine(PerfectMatching * pm, float cost, int node_num, int * oddDegree) {
	//combine minimum-weight matching with the MST to get a multigraph which has vertices with even degree
	//build a tree first. 	
}