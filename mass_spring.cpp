/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;
double c;


/** Custom structure of data to store with Nodes */
struct NodeData {
  Point velocity;  //< Node velocity
  double mass;     //< Node mass
};

struct EdgeData {
  double length;
};

// HW2 #1 YOUR CODE HERE
// Define your Graph type
typedef Graph<NodeData, EdgeData> GraphType;
typedef typename GraphType::node_type Node;
typedef typename GraphType::edge_type Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on Node
 *           at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().velocity * dt;
  }

  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) 
      n.value().velocity += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force being applied to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double) {
    
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
     
      return Point(0,0,0);
    }

    else {
  
      double k = 100.0;

      Point grav_force = {0,0,(-1)*n.value().mass*grav};

      Point spring_force = {0,0,0};
      
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {

        Point disp = ((*it).node1().position() - (*it).node2().position());
        spring_force += -1 * k * (disp/norm(disp)) * (norm(disp) - (*it).value().length);
      }
   
    Point total_force = grav_force + spring_force;
      
    return total_force;
    }
  }
};

struct GravityForce {
  
  Point operator()(Node n, double) {
    
    Point grav_force = {0,0,(-1)*n.value().mass*grav};
      
    return grav_force;
  }
};



struct MassSpringForce {
  
  Point operator()(Node n, double) {
    
    double k = 100.0;
    
    Point spring_force = {0,0,0};
      
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {

      Point disp = ((*it).node1().position() - (*it).node2().position());
      spring_force += -1 * k * (disp/norm(disp)) * (norm(disp) - (*it).value().length);
    }
   
  return spring_force;
  }
};


struct DampingForce {
  
  Point operator()(Node n, double) {
    
    
    Point damping_force = -1.0 * c * n.value().velocity;
   
    return damping_force;
  }
};


struct NoForce {

  Point operator()(Node , double) {

    return {0,0,0};
  }
};

template<typename F1, typename F2, typename F3 = NoForce>
struct make_combined_force {

  F1 force1_;
  F2 force2_;
  F3 force3_;

  make_combined_force(F1 force1, F2 force2) : force1_(force1), force2_(force2) {
  }

  make_combined_force(F1 force1, F2 force2, F3 force3) 
                       : force1_(force1), force2_(force2), force3_(force3) {
  }

  Point operator()(Node n, double t) {

    return (force1_(n, t) + force2_(n, t) + force3_(n, t));
  }
};

int main(int argc, char** argv) {
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a graph
  GraphType graph;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  std::vector<Node> nodes;
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);

      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);

      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // HW2 #1 YOUR CODE HERE
  // Set initial conditions for your nodes, if necessary.
  // Construct Forces/Constraints

  double mass = 1.0/ (double) graph.num_nodes();
  Point initial_velocity = Point(0,0,0);
  c = 1.0 / (double) graph.num_nodes();
  
  for (auto it = graph.node_begin(); it != graph.node_end(); ++it)
    (*it).value() = {initial_velocity, mass};

  for (auto it = graph.edge_begin(); it != graph.edge_end(); ++it)
    (*it).value() = {norm((*it).node1().position()- (*it).node2().position())};

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0.0;
  double t_end   = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;

    symp_euler_step(graph, t, dt, make_combined_force<GravityForce, MassSpringForce, DampingForce>(GravityForce(), MassSpringForce(), DampingForce()));

    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}
