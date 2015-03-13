#ifndef EDGE_HPP
#define EDGE_HPP

class Node;

/**
 * Represents an edge in a graph.
 *
 * Maintains pointers to the vertices that the edge originates 
 * from and terminates at. Edges have both a cost and a length,
 * which are both non-negative integers.
 *
 * Follows value semantics, so can be copy constructed.
 */
class Edge {
  public:
    /**
     * Constructs an Edge from the given parameters.
     */
    Edge(Node *from, Node *to,
            unsigned int cost);

    /**
     * Returns a pointer to the Vertex that this Edge originates
     * from.
     */
    Node *getFrom() const;
    
    /**
     * Returns a pointer to the Vertex that this Edge terminates
     * at.
     */
    Node *getTo() const;

    /**
     * Sets the cost of this Edge.
     */
    void setCost(unsigned int cost);

    /**
     * Returns the cost of this Edge.
     */
    unsigned int getCost() const;
    bool getTraversed() const;
    void setTraversed(bool x);
    void setHide(bool x);
    bool getHide() const;
  


    /*
     * Compares this Edge to another Edge. Suitable for
     * use with a priority queue where Edges with the lowest
     * weight have the highest priority.
     *
     * Returns true if this Edge's cost is more than
     * right's cost.
     */
    bool operator<(const Edge &right) const;
    
  private:
    /**
     * Vertex that this Edge originates from.
     */
    Node *from;

    /**
     * Vertex that this Edge terminates at.
     */
    Node *to;

    /**
     * Cost of this Edge.
     */
    unsigned int cost;
    bool hide;
    bool traversed;

};
#endif