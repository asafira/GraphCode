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
template <typename V>
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


  typedef V node_value_type;
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

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  typedef NodeIterator node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  typedef EdgeIterator edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  typedef IncidentIterator incident_iterator;

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
  }

  /////////////////
  // GRAPH NODES //
  /////////////////

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node>{
   public:
    /** @brief Construct an invalid node.
     */
    Node() {
    }

    /** @brief Return the position of this node
     *  @post Return 3d Point result stored in node
     */
    const Point& position() const {
      return graph_->nodes_[index_].point;
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

    node_value_type& value() {

      return graph_->nodes_[index_].value;
    }

    const node_value_type& value() const {

      return graph_->nodes_[index_].value;
      
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
 
    std::vector<size_type> adj_list;
    
    nodes_.push_back({position, node_value_type(), adj_list});
    
    size_++;
    (void) position;      // Quiet compiler warning
    return Node(this, size_ - 1);        // Invalid node
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW1: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @brief Return the node with index @a i.
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
  

  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }
 
  NodeIterator node_end() const {
    return NodeIterator(this, size());
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

  class Edge : private totally_ordered<Edge>{
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

    size_type index() {

      return edge_index_;
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

    assert (a.index() < size_ && b.index() < size_);
    assert (a.index() != b.index());

    edge_type_ input_edge;

    if (a < b)
      input_edge = {a.index(), b.index()};
    
    else
      input_edge = {b.index(), a.index()};

    size_type uid = insert_into_vector(edges_, input_edge);
    insert_into_vector(nodes_[a.index()].adj_list, b.index());
    insert_into_vector(nodes_[b.index()].adj_list, a.index());

    num_edges_ = edges_.size();
 
    (void) a, (void) b;
    return Edge(this, uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW1: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    return false;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    (void) i;             // Quiet compiler warning
    return Edge(this, i);        // Invalid Edge
  }

  

  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
 
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  ///////////////
  // Iterators //
  ///////////////

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid NodeIterator. */
    NodeIterator() : p_(nullptr) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

     Node* p_;

     Node operator*() const {
        return *p_;
     }

     NodeIterator& operator++() {

       size_type curr_index = p_->index(); 
       curr_index++;

       if (curr_index == graph_->size()) 
         p_ = NULL;
       
       else  
         *p_ = graph_->node(curr_index);
       
       return *this;
     }

     bool operator==(const NodeIterator& other_iter) const {
       return p_ == other_iter.p_;
     }

   private:

    NodeIterator(const Graph* graph, size_type index) 
      : graph_(const_cast<Graph*>(graph)) {

      if ( index < graph->size() ) {
        Node current_node = graph->node(index);
        p_ = &current_node;
      }

      else 
        p_ = NULL;

    }

    friend class Graph;   

    Graph* graph_;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : p_(nullptr) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
   
     Edge* p_;

     Edge operator*() const {
        return *p_;
     }

     EdgeIterator& operator++() {

       size_type curr_index = p_->index(); 
       curr_index++;

       if (curr_index == graph_->num_edges()) 
         p_ = NULL;
       
       else  
         *p_ = graph_->edge(curr_index);
       
       return *this;
     }

     bool operator==(const EdgeIterator& other_iter) const {
       return p_ == other_iter.p_;
     }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    EdgeIterator(const Graph* graph, size_type index) 
      : graph_(const_cast<Graph*>(graph)) {
      
      if ( index < graph->num_edges() ) {
        Edge current_edge = graph->edge(index);
        p_ = &current_edge;
      }

      else
        p_ = NULL;
 
    }

    Graph* graph_;
  };

  // HW1 #3: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef std::ptrdiff_t difference_type;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

 private:

  // Use internal struct of two node indices as edge identifier
  struct edge_type_ {
    size_type index_1;
    size_type index_2;


    // Define < for std::map
    bool operator<(edge_type_ other) const {

        if ( index_1 != other.index_1)
          return index_1 < other.index_1; 

        else
          return index_2 < other.index_2;
    }
  };

  
  // Inspiration for this came from lafstern.org/matt/col1.pdf
  template <class Vector, class T>
  size_type insert_into_vector(Vector& v,  const T& elem) {

    typename Vector::iterator high = std::lower_bound(v.begin(), v.end(), elem);
    if (high == v.end() || elem < *high) {
      v.insert(high, elem);
    }
    return (&(*high) - &v[0]);
  } 

  struct node_type_ {
    Point point;
    node_value_type value;
    std::vector<size_type> adj_list;
  };

  // Keep track of size and num_edges
  size_type size_;
  size_type num_edges_;

  std::vector<node_type_> nodes_;
  std::vector<edge_type_> edges_;

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
};

#endif
