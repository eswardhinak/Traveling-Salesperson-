#include "Edge.hpp"

/*Constructor*/
Edge::Edge(Node * from, Node * to, unsigned int cost){
	this->from = from;
	this->to = to;
	this->cost = cost;
	this->traversed = false;
}

//Accessor Method for To Vertex
Node * Edge::getTo() const{
	return to;
}

Node * Edge::getFrom() const{
	return from;
}

void Edge::setCost(unsigned int cost){
	this->cost = cost;
}


unsigned int Edge::getCost()const{
	return cost;
}

bool Edge::operator<(const Edge &right) const {
	if (cost > right.getCost())
		return true;
	return false;
}

bool Edge::getTraversed() const {
	return traversed;
}

void Edge::setTraversed(bool x) {
	this->traversed = x;
}


bool Edge::getHide() const {
	return hide;
}

void Edge::setHide(bool x) {
	this->hide = x;
}