/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Eges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

#include "Graph.hpp"
#include "Point.hpp"
#include "BoundingBox.hpp"
#include <fstream>

// HW3: YOUR CODE HERE
// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
typedef Graph<bool,bool> GraphType;  //<  DUMMY Placeholder

/** Remove all the nodes in graph @a g whose posiiton is contained within
 * BoundingBox @a bb
 * @post For all i, 0 <= i < @a g.num_nodes(),
 *        not bb.contains(g.node(i).position())
 */
void remove_box(GraphType& g, const BoundingBox& bb) {
  
  for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter) {

    if (bb.contains((*node_iter).position())) {

      node_iter = g.remove_node(node_iter);
    }
  }
  return;
}

class GraphSymmetricMatrix {
public:

  GraphSymmetricMatrix(const GraphType* g) : matrix_params_({g->num_nodes(), g->num_nodes()}), g_(g){};

  template <typename VectorIn, typename VectorOut, typename Assign>
  void mult( const VectorIn& v, VectorOut& w, Assign) const {

    std::size_t node_num = 0;

    double curr_tot;
    while (node_num < g_->num_nodes()) {
      curr_tot = 0.0;
      for (auto adj_iter = g_->node(node_num).edge_begin(); adj_iter != g_->node(node_num).edge_end(); ++adj_iter) {
        
        curr_tot += ((double) v[(*adj_iter).node2().index()])*a_ij(node_num, (*adj_iter).node2().index());
      }

      Assign::apply(w[node_num], curr_tot);
      ++node_num; 
    }
  }
  
  template<typename Vector>
  mtl::vec::mat_cvec_multiplier<GraphSymmetricMatrix, Vector>
  operator*(const Vector& v) const {
    return mtl::vec::mat_cvec_multiplier
                           <GraphSymmetricMatrix, Vector>(*this, v);
  }

  std::size_t num_rows() const {
    return matrix_params_.num_rows;
  }

  std::size_t num_cols() const {
    return matrix_params_.num_cols;
  }

  private:

  int l_ij(std::size_t i, std::size_t j) const {

    if (i == j)
      return -1 * g_->node(i).degree();

    else if (g_->has_edge(g_->node(i), g_->node(j)))
      return 1;

    else
      return 0;
  } 

  int a_ij(std::size_t i, std::size_t j) const {

    if (i == j && g_->node(i).value())
      return 1;
    else if (g_->node(i).value() && g_->node(j).value())
      return 0;
    else
      return l_ij(i,j);
  } 

  struct matrix_params {

    std::size_t num_rows;
    std::size_t num_cols;
  };

  matrix_params matrix_params_;
  const GraphType* g_;
};


  namespace mtl {
    namespace ashape {


      template<>
      struct ashape_aux<GraphSymmetricMatrix> {
        typedef nonscal type;
      };
    }

    template <>
    struct Collection<GraphSymmetricMatrix> {
      typedef double value_type;
      typedef unsigned size_type;
    };
  }

/*  
  template <typename Vector>
  Vector operator*(const Vector& x) const {
    return x;
  }
*/
  inline std::size_t size(const GraphSymmetricMatrix& A) {
    return A.num_rows()*A.num_cols();
  }

  inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
    return A.num_rows();
  }

  inline std::size_t num_cols(const GraphSymmetricMatrix& A) {
    return A.num_cols();
  }


int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }

  // Define an empty Graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<typename GraphType::node_type> node_vec;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    node_vec.push_back(graph.add_node(2*p - Point(1,1,0)));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    graph.add_edge(node_vec[t[0]], node_vec[t[1]]);
    graph.add_edge(node_vec[t[0]], node_vec[t[2]]);
    graph.add_edge(node_vec[t[1]], node_vec[t[3]]);
    graph.add_edge(node_vec[t[2]], node_vec[t[3]]);
  }

  // Get the edge length, should be the same for each edge
  double h = graph.edge(0).length();

  // Make holes in our Graph
  remove_box(graph, BoundingBox(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // HW3: YOUR CODE HERE
  // Define b using the graph, f, and g.
  // Construct the GraphSymmetricMatrix A using the graph
  // Solve Au = b using MTL.

  Point p_p = Point({0.6, 0.6, 0});
  Point p_m = Point({-0.6, -0.6, 0});

  Point bb_p1 = Point({-0.6, -0.2, -1.0});
  Point bb_p2 = Point({0.6, 0.2, 1});
  BoundingBox bb = BoundingBox(bb_p1, bb_p2);

  for (auto node_iter = graph.node_begin(); node_iter != graph.node_end(); ++node_iter) {

    auto node = *node_iter;

    if (norm_inf(node.position()) == 1)
      node.value() = true;

    else if (norm_inf(node.position() - p_p) < 0.2 || norm_inf(node.position() - p_m) < 0.2)
      node.value() = true;
 
    else if (bb.contains(node.position()))
      node.value() = true;

    else
      node.value() = false;
  }

  mtl::dense_vector<double> b(graph.size());
  
  
  for (auto node_iter = graph.node_begin(); node_iter != graph.node_end(); ++node_iter) {

    auto node = *node_iter;

    if (norm_inf(node.position()) == 1)
      b[node.index()] = 0;

    else if (norm_inf(node.position() - p_p) < 0.2 || norm_inf(node.position() - p_m) < 0.2)
      b[node.index()] = -0.2;
 
    else if (bb.contains(node.position()))
      b[node.index()] = 1;

    else {
        
      b[node.index()] = h*h*5*std::cos(norm_1(node.position()));
      for (auto inc_iter = node.edge_begin(); inc_iter != node.edge_end(); ++inc_iter) {

        if (norm_inf(node.position() - p_p) < 0.2 || norm_inf(node.position() - p_m) < 0.2)
          b[node.index()] -= -0.2;
 
        if (bb.contains(node.position()))
          b[node.index()] -= 1;
      }
    }
  }

  GraphSymmetricMatrix A = GraphSymmetricMatrix(&graph);
  itl::pc::identity<GraphSymmetricMatrix> pre_con(A);
  mtl::dense_vector<double> x(graph.size(), 1.0);

  itl::noisy_iteration<double> iter(b, 500, 1.e-6);

  itl::cg(A, x, b, pre_con, iter);

  return 0;
}
