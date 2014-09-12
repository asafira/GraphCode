#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CS207/Util.hpp"
#include "Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  /////////////////////////////
  // PUBLIC TYPE DEFINITIONS //
  /////////////////////////////

  /** Type of this graph. */
  typedef Graph graph_type;



  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  typedef unsigned size_type;

  ////////////////////////////////
  // CONSTRUCTOR AND DESTRUCTOR //
  ////////////////////////////////

  /** @brief Construct an empty graph. */
  Graph() 
    : size_(0), num_edges_(0),  nodes_() {
  }
  /** Default destructor */
  ~Graph() = default;

  /////////////
  // General //
  /////////////

  /** @brief Return the number of nodes in the graph.
   *  @post return number of nodes in Graph
   *  Complexity: O(1).
   */
  size_type size() const {
   
    return size_;
  }

  /** @brief Remove all nodes and edges from this graph.
   *  @post num_nodes() == 0 && num_edges() == 0
   *
   *  Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    num_edges_ = 0;
    size_ = 0;
    nodes_.clear();
    edges_.clear();
    index_lookup_.clear();
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   public:
    /** @brief Construct an invalid node.
     */
    Node() {
    }

    /** @brief Return the position of this node
     *  @post Return 3d Point result stored in node
     */
    const Point& position() const {
      return graph_->nodes_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      
      return index_;
    }

    /** Test whether this node and @a x are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& x) const {
      
      (void) x;          // Quiet compiler warning
      return (this->graph_ == x.graph_) && (this->index() == x.index());
    }

    /** Test whether this node is less than @a x in the global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& x) const {
     
      (void) x;           // Quiet compiler warning
      return this->index() < x.index();
    }

   private:

    Node (const Graph* graph, size_type index) 
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
    
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Define a 
    Graph* graph_;
    size_type index_;
  };

  /** @brief Synonym for size()
   *  @post returns the number of nodes in the graph
   */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {

    nodes_[size_] = position;
    
    size_++;
    (void) position;      // Quiet compiler warning
    return Node(this, size_ - 1);        // Invalid node
  }

  /** @brief Return the node with index @a i.
   *  @pre 0 <= @a i < num_nodes()
   *  @post result_node.index() == i
   *
   *  This class has complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
    
    (void) i;             // Quiet compiler warning
    return Node(this, i);        // Invalid node
  }

  /////////////////
  // GRAPH EDGES //
  /////////////////

  /** @class Graph::Edge
   *  @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */

  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      size_type node_1_index = graph_->edges_[edge_index_].index_1;
      return Node(graph_, node_1_index);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      size_type node_2_index = graph_->edges_[edge_index_].index_2;
      return Node(graph_, node_2_index);      // Invalid Node
    }

    /** @brief Tests whether this edge and @a x are equal.
     *  @pre 0 <= @x <= num_edges()
     *  @pre @x is size_type
     *  @post result is true if this edge == edge x and false otherwise
     */
    bool operator==(const Edge& x) const {
 
      if (this->graph_==x.graph_) {

        if (this->node1() == x.node1() && this->node2() == x.node2())
          return true;

        else if (this->node2() == x.node1() && this->node1() == x.node2())
          return true;
      }

      (void) x;          // Quiet compiler warning
      return false;
    }

    /** @brief Defines the less than operator for comparing edges
     *  @post Result is a boolean value such that the edges obey trichotomy
     */
    bool operator<(const Edge& x) const {
      size_type this_index_1 = graph_->edges_[edge_index_].index_1;
      size_type x_index_1 = x.graph_->edges_[edge_index_].index_1;
      (void) x;           // Quiet compiler warning
      return this_index_1 < x_index_1; 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    Edge(const Graph* graph, size_type index)
       : edge_index_(index), graph_(const_cast<Graph*>(graph)) {
    } 

    friend class Graph;

    // Keep edge index and graph pointer private variables
    size_type edge_index_;
    Graph* graph_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** @brief Add an edge to the graph, or return the current edge 
   *  if it already exists.
   *  @pre @a a and @a b are distinct valid nodes of this graph
   *  @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   *  @post has_edge(@a a, @a b) == true
   *  @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   *  Can invalidate edge indexes -- in other words, old edge(@a i) might not
   *  equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   *  Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    
    // Create internal edge_nodes struct
    edge_nodes_* edge_ab = new edge_nodes_();
    edge_ab->index_1 = a.index();
    edge_ab->index_2 = b.index();

    // Check if edge already exists
    if (index_lookup_.find(*edge_ab) != index_lookup_.end())
      return Edge(this, index_lookup_[*edge_ab]);

    // Create edge with reversed nodes
    edge_nodes_* edge_ba = new edge_nodes_();
    edge_ba->index_1 = b.index();
    edge_ba->index_2 = a.index();

    // Check if reverse edge exists
    if (index_lookup_.find(*edge_ba) != index_lookup_.end()) {
      
      // Substitute reversed edge into graph, return it
      edges_[index_lookup_[*edge_ba]] = *edge_ab;
      return Edge(this, index_lookup_[*edge_ba]);
    }
   
    // Create new edge, incrementin num_edges and next_edge_index appropriately
    num_edges_++;
    edges_[num_edges_ - 1] = *edge_ab;
   
    (void) a, (void) b;   // Quiet compiler warning
    return Edge(this, num_edges_ - 1);        // Invalid Edge
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    (void) i;             // Quiet compiler warning
    return Edge(this, i);        // Invalid Edge
  }


 private:

  // Use internal struct of two node indices as edge identifier
  struct edge_nodes_ {
    size_type index_1;
    size_type index_2;

    // Define < for stl::map
    bool operator<(edge_nodes_ other) const {

        if ( index_1 != other.index_1)
          return index_1 < other.index_1; 

        else
          return index_2 < other.index_2;
    }
  };

  // Keep track of size and num_edges
  size_type size_;
  size_type num_edges_;

  // Keep maps of identifying size_types to Points and edge_nodes
  std::map<size_type, Point> nodes_;
  std::map<size_type, edge_nodes_> edges_;

  // Keep revese-lookup map for add_edge speedup
  std::map<edge_nodes_, size_type> index_lookup_;

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
};

#endif
