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

// NoCons is the default (or null) constraint that does not constrain the points in graph.
struct NoCons {

  void operator()(GraphType&, double) {return;}
};

template <typename G, typename F, typename C = NoCons>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the {n+1} node positions
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().velocity * dt;
  }
  
  // Enforce constraints
  constraint(g, t);

  // Compute the {n+1} node velocities
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
//    if (n.position() != Point(0,0,0) && n.position() != Point(1,0,0)) 
    n.value().velocity += force(n, t) * (dt / n.value().mass);

  }

  return t + dt;
}

/* @brief StickyPlaneCons is a planar constraint on a set of nodes.
   @param[in] @a coeffs is a Point object describing the coefficients of a plane
              in standard form.
   @param[in] @a offset is the Plane's offset from the origin
   See operator() description in struct for usage info
*/
struct StickyPlaneCons {
 

  StickyPlaneCons(Point coeffs, double offset)
                  : coeffs_(coeffs), offset_(offset) {}

  /* @brief Reinforces the positions of nodes in the graph g stay on one side of the 
            plane given by dot(Point(x,y,z), coeffs_) - offset_ = 0 
 
     @param[in] Graph g is the graph to 
   
     @pre plane must pass through the z axis
 
     @post nodes in graph g that have passed the plane described by dot(n.position(),
           coeffs_) - offset_ = 0 will have their positions reset to the position on the 
           plane with correct z coordinate


   */
  void operator()(GraphType& graph, double) {

    for(auto it =  graph.node_begin(); it != graph.node_end(); ++it) {
    
      auto n = *it;
      if (dot(n.position(), coeffs_) < offset_) {

        // Finds the point on plane closest to n.position()
        Point point_on_plane = Point(0,0,offset_/coeffs_.z);
        Point v = n.position() - point_on_plane;
        double dist = dot(v, coeffs_)/norm(coeffs_); 
        n.position() += -1.0 * dist * coeffs_/norm(coeffs_);
        n.value().velocity = Point(n.value().velocity.x, n.value().velocity.y, 0);
      } 
    }
  }

  private:
    Point coeffs_;
    double offset_;
};

/* @brief SphereCons is a functor that constrains a graph g to only contain points outside of a 
          sphere centered at @a center and with radius @a r. It moves all within the sphere 
          to the closest points on the sphere.

   @param[in] Point center defines the center of the sphere
   @param[in] double r defines the radius of the sphere

 */
struct SphereCons {

  SphereCons(Point center, double radius) : c(center), r(radius) {}

  /* @brief the operator() method will take in a graph and the time and impose
            the sphere constraint

     @post ALl nodes in graph with positions that were outside or on of the sphere 
           will be unmodified,
     @post All nodes in the graph within the sphere will move the point on the sphere that is closest
           to its position
   */
  void operator()(GraphType& graph, double) {

    for( auto it = graph.node_begin(); it != graph.node_end(); ++it) {
    
      auto n = *it;
      if (norm(n.position() - c) < r) {

        Point force_direc = (n.position() - c)/norm(n.position()-c);
        n.position() = c + r*force_direc;
        n.value().velocity -= dot(n.value().velocity, force_direc)*force_direc;
      }
    }
  }

  private:
    Point c;
    double r;
};


// Similar to SphereCons, except that positions within the sphere will have their nodes deleted
struct CuttingSphereCons {

  CuttingSphereCons(Point center, double radius) : c(center), r(radius) {}

  void operator()(GraphType& graph, double) {

    auto it = graph.node_begin();
    while(it != graph.node_end()) {

      if (norm((*it).position() - c) < r)
        it = graph.remove_node(it);

      else
        ++it;
    }
  }

  private:
    Point c;
    double r;
};


 
/* Combined_cons combines the constraints in the graph */
template<typename C1 = NoCons, typename C2 = NoCons, typename C3 = NoCons>
struct combined_cons {

  // Has constructors for 1, 2, or 3 constraints

  combined_cons(C1 cons1) : cons1_(cons1) {
  }
  
  combined_cons(C1 cons1, C2 cons2) : cons1_(cons1), cons2_(cons2) {
  }

  combined_cons(C1 cons1, C2 cons2, C3 cons3) 
                       : cons1_(cons1), cons2_(cons2), cons3_(cons3) {
  }

  //performs the constraints sequentially
  void operator()(GraphType& g, double t) {

    cons1_(g,t);
    cons2_(g,t);
    cons3_(g,t);

    return;
  }

  private:
    C1 cons1_;
    C2 cons2_;
    C3 cons3_;
};



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


/* @brief GravityForce takes in a node n and a time and returns the force on that node
          given a defined mass in its value() struct

   @param[in] node n is a valid node
   @param[in] t is the current timestep
   @pre n.value() must contain a parameter mass of type double.

   @post returns the 3D Point that is the gravitational force on that node.
*/
struct GravityForce {
  
  Point operator()(Node n, double) {
    
    Point grav_force = {0,0,(-1)*n.value().mass*grav};
      
    return grav_force;
  }
};


/* @brief MassSpringForce returns the spring force on a node in the graph

   @param[in] node n is a node in the graph
   @param[in] time t is the current timestep
   @param[out] Point that is the spring force on the node
   
   @pre All edges in the graph connected to n have had their rest lengths 
        put into a struct variable length stored in edge.value()
*/
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



/* @brief DampingForce returns the damping force on a node in the graph

   @param[in] node n is a node in the graph
   @param[in] time t is the current timestep
   @param[out] Point that is the damping force on the node
   
*/
struct DampingForce {
  
  Point operator()(Node n, double) {
    
    
    Point damping_force = -1.0 * c * n.value().velocity;
   
    return damping_force;
  }
};

/* @brief The null force to act as the default force in make_combined_force */
struct NoForce {

  Point operator()(Node , double) {

    return {0,0,0};
  }
};

/* @brief Combines forces in its arguements (up to 3) to return a single total force
          on a node

*/
template<typename F1, typename F2 = NoForce, typename F3 = NoForce>
struct make_combined_force {

  F1 force1_;
  F2 force2_;
  F3 force3_;

  make_combined_force(F1 force1) : force1_(force1) {
  }
  
  make_combined_force(F1 force1, F2 force2) : force1_(force1), force2_(force2) {
  }

  make_combined_force(F1 force1, F2 force2, F3 force3) 
                       : force1_(force1), force2_(force2), force3_(force3) {
  }

  // Return the sum of the forces
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

    // Note that we need to include the types of forces for the combined functors
    symp_euler_step(graph, t, dt, 
                    make_combined_force<GravityForce, MassSpringForce, DampingForce>
                    (GravityForce(), MassSpringForce(), DampingForce()),
                    combined_cons<StickyPlaneCons, NoCons>
                    (StickyPlaneCons(Point(0,0,1), -0.75), NoCons()));

    // Clear the viewer's nodes and edges 
    viewer.clear();
    node_map.clear();

    // Update viewer with nodes' new positions 
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      CS207::sleep(0.001);
  }

  return 0;
}; 
