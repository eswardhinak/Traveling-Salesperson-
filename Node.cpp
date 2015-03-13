/* File Name: Node.cpp
 * Author: Eswar Dhinakaran
 * Description: This contains the class description for a Node. The Node
 * 		has member variables of distance (for Dijkstra), and name
 * 		It also contains a unordered map of edges that are connected
 * 		to this Node.
 * Last Modified: 28 November 2014
 *
*/


#include "Node.hpp"
#include <iostream>
using namespace std;

/*ew Constructor */
Node::Node(int name){
	this->name = name; //sets the name variable 
}
int Node::getName() const{
	return name;
}
bool Node::addEdge(Node *to, unsigned int cost){
	Edge * e= new Edge(this, to, cost); //creates a RTS edge
	int adj_name = to->getName(); //gets name of connected Node
	//creates pair
	std::pair<int, Edge*> edge_pair (adj_name, e);
	edges.push_back(edge_pair); //inserts into vector
	return true;
}



/*Method to calculate total cost of all edges connected to this Node*/
unsigned int Node::totalEdgeCost() const{
	auto it = edges.begin();//iterator to edges map
	unsigned int total_cost = 0; //initialize total cost
	//iterate through edges and get cost
	for (; it != edges.end(); ++it){
		total_cost += (it->second)->getCost();
	}
	return total_cost; 
}

/*Accessor method to get the Distance of this Node*/
unsigned int Node::getDistance() const{
	return distance;
}

/*Setter method to set the distance of this Node*/
void Node::setDistance(unsigned int distance){
	this->distance = distance;
}

/*Accessor method used in Prim/Dijkstra to check the state of this Node*/
bool Node::wasVisited() const {
	return visited;
}

/*Accessor method used in Prim/Dijkstra to set the state of this Node*/
void Node::setVisited(bool visited){
	this->visited = visited;
}

/*Method to clear all the edges*/
void Node::clearEdges(){
	edges.clear();
}

/*Accessor method to get the edges map from Node*/
const std::unordered_map<int, Edge*>& Node::getEdges() const{
	const unordered_map<int, Edge*>& map_ref = edges;
	return map_ref;
}
