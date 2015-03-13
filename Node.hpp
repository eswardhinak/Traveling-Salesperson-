#ifndef NODE_HPP
#define NODE_HPP

#include <string>
#include <vector>
#include "Edge.hpp"

/**
 * Represents a Node in a graph.
 *
 * Vertices are connected to other Vertices via Edges. Each Node
 * maintains a collection of all Edges that originate from it.
 */
class Node {
    // Graph needs access to Edge map for Dijkstra/Prim algorithms.
    friend class MST;
    
  public:
    /**
     * Initialize the Node with the given name.
     */
    Node(int name);

    /**
     * Add an edge to this Node. If an edge already exists to the given
     * Node, updates the cost and length of the edge to match the
     * passed parameters.
     */
    int getName() const;
    bool addEdge(Node *to, unsigned int cost);
    
    /**
     * Returns the Node's name.
     */
    
    /**
     * Gets the Node's distance value.
     */
    unsigned int getDistance() const;

    /**
     * Sets the Node's distance value.
     */
    void setDistance(unsigned int distance);
    
    /**
     * Gets the Node's visited state.
     */
    bool wasVisited() const;

    /**
     * Sets the Node's visited state.
     */
    void setVisited(bool visited);

    /**
     * Clears all edges from this Node.
     */
    void clearEdges();

    /**
     * Gets total cost of all edges terminating at this Node.
     */
    unsigned int totalEdgeCost() const;
    const std::vector<std::pair<int, Edge*> &getEdges() const;

    
  private:
    /**
     * Returns a reference to the internal map of Edges.
     * Used by UndirectedGraph for Dijkstra/Prim algorithms.
     */

    /**
     * Name of this Node.
     */
    int name;
    
    /**
     * Distance of this Node from initial Node.
     * Used by Dijkstra's algorithm.
     */
    unsigned int distance;
    
    /**
     * Whether this node has been visited.
     * Used by Prim's algorithm.
     */
    bool visited;

    /**
     * Map of adjacent Node name to Edge describing the adjacency.
     */
    std::vector<std::pair<int, Edge*> edges;
};

#endif
