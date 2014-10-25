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


    double curr_tot;
    for (auto node_iter = g_->node_begin(); node_iter != g_->node_end(); ++node_iter){
      curr_tot = 0.0;
      auto node = *node_iter;
      curr_tot += v[node.index()] * a_ij(node.index(), node.index());

      for (auto adj_iter = node.edge_begin(); adj_iter != node.edge_end(); ++adj_iter) {
        
        curr_tot += (v[(*adj_iter).node2().index()])*a_ij(node.index(), (*adj_iter).node2().index());
      }

      Assign::apply(w[node.index()], curr_tot);
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

  double l_ij(std::size_t i, std::size_t j) const {

    if (i == j)
      return -1.0 * (double) g_->node(i).degree();

    else if (g_->has_edge(g_->node(i), g_->node(j)))
      return 1.0;

    else
      return 0.0;
  } 

  double a_ij(std::size_t i, std::size_t j) const {

    if (i == j && g_->node(i).value())
      return 1.0;
    else if (g_->node(i).value() || g_->node(j).value())
      return 0.0;
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

  inline std::size_t size(const GraphSymmetricMatrix& A) {
    return A.num_rows()*A.num_cols();
  }

  inline std::size_t num_rows(const GraphSymmetricMatrix& A) {
    return A.num_rows();
  }

  inline std::size_t num_cols(const GraphSymmetricMatrix& A) {
    return A.num_cols();
  }



template <class Real, class OStream = std::ostream>
  class visual_iteration : public itl::cyclic_iteration<Real> 
  {
      typedef itl::cyclic_iteration<Real> super;
      typedef visual_iteration self;

    public:
  
      template <class Vector>
      visual_iteration(const Vector& r0, int max_iter_, Real tol_, GraphType* graph_, 
                       mtl::dense_vector<double>* u_, Real atol_ = Real(0), 
                       int cycle_ = 100, OStream& out = std::cout)
        : super(r0, max_iter_, tol_, atol_, cycle_, out), graph(graph_), u(u_) {

        viewer.launch();
        auto node_map = viewer.empty_node_map(*graph);
       
        viewer.add_nodes(graph->node_begin(), graph->node_end(), CS207::NodeColor(*u), CS207::VectorZPosition(*u), node_map);
        viewer.add_edges(graph->edge_begin(), graph->edge_end(), node_map);
        viewer.center_view();
        CS207::sleep(0.1);
      }

      bool finished() { return super::finished(); }

      template <typename T>
      bool finished(const T& r) 
      {
         bool ret= super::finished(r);
         viewer.clear(); 
         auto node_map = viewer.empty_node_map(*graph);

         viewer.add_nodes(graph->node_begin(), graph->node_end(), CS207::NodeColor(*u), CS207::VectorZPosition(*u), node_map);
         viewer.add_edges(graph->edge_begin(), graph->edge_end(), node_map);
         viewer.center_view();
         viewer.set_label(this->i);  
         CS207::sleep(0.1);
         return ret;
      }

    protected:
      GraphType* graph; 
      CS207::SDLViewer viewer;
      mtl::dense_vector<double>* u;

  };
int main(int argc, char** argv) {

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

  Point p_pp = Point({0.6, 0.6, 0});
  Point p_mm = Point({-0.6, -0.6, 0});
  Point p_pm = Point({0.6, -0.6, 0});
  Point p_mp = Point({-0.6, 0.6, 0});

  Point bb_p1 = Point({-0.6, -0.2, -1.0});
  Point bb_p2 = Point({0.6, 0.2, 1});
  BoundingBox bb = BoundingBox(bb_p1, bb_p2);

  for (auto node_iter = graph.node_begin(); node_iter != graph.node_end(); ++node_iter) {

    auto node = *node_iter;

    if (norm_inf(node.position()) == 1)
      node.value() = true;

    else if (norm_inf(node.position() - p_pp) < 0.2 || norm_inf(node.position() - p_mm) < 0.2)
      node.value() = true;
 
    else if (norm_inf(node.position() - p_pm) < 0.2 || norm_inf(node.position() - p_mp) < 0.2)
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

    else if (norm_inf(node.position() - p_pp) < 0.2 || norm_inf(node.position() - p_mm) < 0.2)
      b[node.index()] = -0.2;
 
    else if (norm_inf(node.position() - p_pm) < 0.2 || norm_inf(node.position() - p_mp) < 0.2)
      b[node.index()] = -0.2;
 
    else if (bb.contains(node.position()))
      b[node.index()] = 1;

    else {
        
      b[node.index()] = h*h*5*std::cos(norm_1(node.position()));
      for (auto inc_iter = node.edge_begin(); inc_iter != node.edge_end(); ++inc_iter) {

        auto edge = *inc_iter;
        if (norm_inf(edge.node2().position() - p_pp) < 0.2 || norm_inf(edge.node2().position() - p_mm) < 0.2)
          b[node.index()] -= -0.2;
 
        else if (norm_inf(edge.node2().position() - p_pm) < 0.2 || norm_inf(edge.node2().position() - p_mp) < 0.2)
          b[node.index()] -= -0.2;
        
        else if (bb.contains(edge.node2().position()))
          b[node.index()] -= 1;

      }
    }
  }

  GraphSymmetricMatrix A = GraphSymmetricMatrix(&graph);


  mtl::dense_vector<double> u(graph.size(), 0.0);

  visual_iteration<double> iter(b, 500, 1.e-10, &graph, &u, 0, 50);
  itl::cg(A, u, b, iter);

  return 0;
}
