#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP
/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <queue>
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


  /** @brief Performs a Breadth-first search on an input
   *  graph @a g, assigning the distance each node is
   *  from the node closest to @a point and outputting
   *  the max distance.
   *  @post All nodes have a value v s.t. it takes at least
   *  v edges to traverse before reaching the node closest
   *  to @a point.
  */
  int shortest_path_lengths(Graph<int>& g, const Point& point) {
  
    // Find the node closest point and give it value 0
    Node root = *std::min_element(g.node_begin(), g.node_end(), MyComparator(point));
    root.value() = 0;
    
    // Initialize a queue for a graph BFS
    std::queue<Graph::Node> nodes_to_visit;
    nodes_to_visit.push(root);
    
    // Initialize a boolean array to keep track of what nodes we visited during BFS
    std::vector<bool> visited(g.size(), false);
    visited[(root).index()] = true;

    Node curr_node;
    int max = 0;

    // Loop through nodes until queue is empty  
    while (!nodes_to_visit.empty()) {
      
      // Pop a node from the queue
      curr_node = nodes_to_visit.front();
      nodes_to_visit.pop();
      
      // Update max accordingly
      if (max < curr_node.value())
        max = curr_node.value();
     
      // Iterate through all adjacent nodes near the current node
      for( incident_iterator it = curr_node.edge_begin(); it != curr_node.edge_end(); ++it) {
        if (!visited[(*it).node2().index()]){
          
          // Add neighboring nodes to queue, update visisted vector, assign value
          nodes_to_visit.push((*it).node2());
          visited[(*it).node2().index()] = true;
          (*it).node2().value() = curr_node.value() + 1; 

        }
      }
    }

    // Assign a value of -1 for each unvinisted node.
    for ( auto node_index : visited) {

      if( !visited[node_index])
        g.node(node_index).value() = -1;
    } 

    // Return the max separation in the graph
    return max;
  }

  /* @brief Compares two node objects and returns
   * one closest to point @a p
   */
  struct MyComparator {
    Point p_;
    MyComparator(const Point& p) : p_(p) {
    };

    template <typename NODE>
    bool operator()(const NODE& node1, const NODE& node2) const {
      return (normSq(p_ - node1.position()) < normSq(p_ - node2.position()));
      }
    };

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
      graph_ = NULL;
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
      
      return this->index() < x.index();
    }
     
    // Return value of current node
    node_value_type& value() {
      
      return graph_->nodes_[index_].value;
    }

    // Return immutable value of this node
    const node_value_type& value() const {
      
      return graph_->nodes_[index_].value;
    }

    // Return the degree of this node
    size_type degree() const {
      return graph_->nodes_[index_].adj_list.size();
    }
    
    // Return the first iterator over all nodes adjacent to this node
    incident_iterator edge_begin() const {
    
      return IncidentIterator(this, 0);
    }

    // Return the last iterator over all nodes adjacent to this node
    incident_iterator edge_end() const {

      return IncidentIterator(this, degree());
    }
   
  private:

    // Construct a valid Node
    Node (const Graph* graph, size_type index) 
      : graph_(const_cast<Graph*>(graph)), index_(index) {
    }
    
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // Class variables for Node
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
  Node add_node(const Point& position, const node_value_type& myvalue= node_value_type()) {
 
    std::vector<size_type> adj_list;
    
    nodes_.push_back({position, myvalue, adj_list});
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
    
    if (n == node(n.index()))
      return true;

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
    return Node(this, i);        // Invalid node
  }
  

  // Return an iterator pointing to the first node
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
 
  // Return the iterator pointing after all nodes
  node_iterator node_end() const {
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
      graph_ = NULL;
    }

    /** Return a node of this Edge */
    Node node1() const {
      
      size_type index_1;

      // Checks which node needs to be outputted as node1()
      if (nodes_flipped_) 
        index_1  = graph_->edges_[edge_index_].index_2;

      else
        index_1  = graph_->edges_[edge_index_].index_1;

      return Node(graph_, index_1);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
 
      size_type index_2;

      // Checks which node needs to be outputted as node2()
      if (nodes_flipped_)
        index_2 = graph_->edges_[edge_index_].index_1;

      else
        index_2 = graph_->edges_[edge_index_].index_2;

      return Node(graph_, index_2);      // Invalid Node
    }

    /** @brief Tests whether this edge and @a x are equal.
     *  @pre 0 <= @x <= num_edges()
     *  @pre @x is size_type
     *  @post result is true if this edge == edge x and false otherwise
     */
    bool operator==(const Edge& x) const {

      // Check if part of the same graph 
      if (this->graph_==x.graph_) {

        // Check if has same nodes
        if (this->node1() == x.node1() && this->node2() == x.node2())
          return true;

        else if (this->node2() == x.node1() && this->node1() == x.node2())
          return true;
      }

      return false;
    }

    /** @brief Defines the less than operator for comparing edges
     *  @post Result is a boolean value such that the edges obey trichotomy
     */
    bool operator<(const Edge& x) const {
      size_type this_index_1 = graph_->edges_[edge_index_].index_1;
      size_type x_index_1 = x.graph_->edges_[edge_index_].index_1;
      
      return this_index_1 < x_index_1; 
    }

   private:
    
    // Allow Graph to access Edge's private member data and functions.
    Edge(const Graph* graph, size_type index)
       : edge_index_(index), graph_(const_cast<Graph*>(graph)), nodes_flipped_(false) { 
    }

    // Output the index of the current edge
    size_type index() {

      return edge_index_;
    }

    // "flip" the node output of the current edge
    void flip_nodes() {

      nodes_flipped_ = !nodes_flipped_;
    }

    // Allow graph to access Edge privates
    friend class Graph;

    // Keep edge index and graph pointer private variables
    size_type edge_index_;
    Graph* graph_;
    bool nodes_flipped_;
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
 
    // Check that nodes are valid
    assert (a.index() < size() && b.index() < size());
    assert (a.index() != b.index());

    // Iterate through adjacency list of a to check if Edge already exists
    std::vector<size_type> a_adj_list = nodes_[a.index()].adj_list; 
    for (std::vector<size_type>::iterator it = a_adj_list.begin(); it != a_adj_list.end();  it++) {

      if (Edge(this, *it).node2() == b) 
        return Edge(this, *it);

      if (Edge(this, *it).node1() == b) {
        Edge new_edge = Edge(this, *it);
        // Flip edges before returning abide by specifications
        new_edge.flip_nodes();
        return new_edge;
      }
    }

    // Create uid for new edge
    size_type uid = num_edges_;

    //Add nodes into adjacency lists as new edge in edges vector
    nodes_[b.index()].adj_list.push_back(uid);
    nodes_[a.index()].adj_list.push_back(uid);
    edges_.push_back({a.index(), b.index()});

    num_edges_++;

    return Edge(this, uid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    
    // Iterate through adjacency list for a to find if edge exists 
    std::vector<size_type> a_adj_list = nodes_[a.index()].adj_list; 
    for (std::vector<size_type>::iterator it = a_adj_list.begin(); it != a_adj_list.end();  it++) {

      if (Edge(this, *it).node2() == b || Edge(this, *it).node1()==b) 
        return true;
      }

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

  
  // Returns iterator to first edge in this graph
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
 
  // Returns iterator after last edge in this graph
  edge_iterator edge_end() const {
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
    NodeIterator() : graph_(NULL) {
    }

     /* @brief Dereferences current iterator, outputs Node
      * @pre iterator must point towards 
      */
     Node operator*() const {
        return graph_->node(curr_node_index_);
     }

     /* @brief
      *
      */
     node_iterator& operator++() {
        curr_node_index_++;
        if (curr_node_index_ >= graph_->size()) {
          curr_node_index_ = graph_->size();
        }
        return *this;

     }

     bool operator==(const node_iterator& other_iter) const {
       if (graph_ == other_iter.graph_)
         if (curr_node_index_ == other_iter.curr_node_index_)
           return true;
       return false;
     }


   private:
 
    /* @brief Constructs new iterator for nodes in
     * Graph @a *graph, pointing at the node with 
     * index @a index.
     * @pre 0 <= index <= graph->size();
     */
    NodeIterator(const Graph* graph, size_type index) 
      : graph_(const_cast<Graph*>(graph)) {

      assert( index <= graph_->size() );

      curr_node_index_ = index;
    }


    Graph* graph_;
    size_type  curr_node_index_;
    
    friend class Graph;   
  };

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
    EdgeIterator() : graph_(NULL) {
    }

     /* @brief Returns an Edge the iterator is pointing
      * to.
      * @pre current index may not exceed number of edges
      */
     Edge operator*() const {
        
        assert (curr_edge_index_ < graph_->num_edges_);
        return Edge(graph_, curr_edge_index_);
     }

     /* @brief Incrememnt operator, moves iterator to next
      * edge. 
      * @post If *it == edge_end(), it++ == edge_end()
      */
     edge_iterator& operator++() {
       if (curr_edge_index_ != graph_->num_edges()) {
         curr_edge_index_++;
       }
       return *this;
     }
   
     /* @brief Checks if two iterators point to the same Edge
      */
     bool operator==(const edge_iterator& other_iter) const {
       return (graph_ == other_iter.graph_ && 
                curr_edge_index_ == other_iter.curr_edge_index_);
     }

   private:

    /* @brief Constructs an iterator over the edges
     * in @a *graph pointing at the edge with index
     * @index
     * @pre 0 <= index <=  graph->num_edges_;
     */
    EdgeIterator(const Graph* graph, size_type index) 
      : graph_(const_cast<Graph*>(graph)), curr_edge_index_(index) {
      
      // Check if valid index
      assert( index <= graph->num_edges() );
    }

    Graph* graph_;
    size_type curr_edge_index_;
    
    friend class Graph;
  };

  // HW1 #3: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const


  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    IncidentIterator() : node_(NULL), graph_(NULL){
    }

     /* @brief returns the value the iterator is pointing to.
      * @pre iterator must be pointing towards an edge with
      * index @curr_edge_index_, 0 <= curr_edge_index_ <= 
      * graph_->num_edges()
      */
     Edge operator*() const {
        assert ( curr_edge_index_ < node_->degree() );

        Edge curr_edge = Edge(graph_, graph_->nodes_[node_->index()].adj_list[curr_edge_index_]);
        if (curr_edge.node1() != *node_)
          curr_edge.flip_nodes();
        return curr_edge;
     }

     // @brief Returns an iterator to the next edge
     incident_iterator& operator++() {
       curr_edge_index_++; 
       if (curr_edge_index_ >= node_->degree())  
         curr_edge_index_ = node_->degree();

       return *this;
     }
    
     // @brief Checks if graph, nodes, and edge matches for given iterators
     bool operator==(const incident_iterator& other_iter) const {
       if (graph_ == other_iter.graph_)
          if(node_ == other_iter.node_)
            if (curr_edge_index_ == other_iter.curr_edge_index_)
              return true;
     
       return false;
     }

   private:
    friend class Graph;

    /* @brief Constructs an iterator pointing at the edge
     * with position @a index in the adjacency list of node.
     * @ pre 0 <= index <= node.degree()
     */
    IncidentIterator(const Node* node, size_type index) 
      : node_(const_cast<Node*>(node)) {
      
      graph_ = const_cast<Graph*>(node_->graph_);

      assert( index <= node_->degree() );
 
      curr_edge_index_ = index;
 
    }

    Graph* graph_;
    Node* node_;
    size_type curr_edge_index_;
  };

 private:

  // Use internal struct of two node indices as edge identifier
  struct edge_type_ {
    size_type index_1;
    size_type index_2;


    // Define < for std::map
    bool operator<(edge_type_ e2) const {
      return std::tie(index_1, index_2) < std::tie(e2.index_1, e2.index_2);
    }
  };

  
 /* // Inspiration for this came from lafstern.org/matt/col1.pdf
  template <class Vector, class T>
  void insert(Vector& v,  const T& elem) {

    typename Vector::iterator high = std::lower_bound(v.begin(), v.end(), elem);

    if (high == v.end() || elem < *high) {
      v.insert(high, elem);
      //return curr_pos;
      //return (&(*curr_pos) - &v[0])i;
    }
    
  } */

  // Struct to hold node data
  struct node_type_ {
    Point point;
    node_value_type value;
    std::vector<size_type> adj_list;
  };

  // Keep track of size and num_edges
  size_type size_;
  size_type num_edges_;

  // Vectors to hold nodes and edges
  std::vector<node_type_> nodes_;
  std::vector<edge_type_> edges_;

  // Disable copy and assignment of a Graph
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
};

#endif
