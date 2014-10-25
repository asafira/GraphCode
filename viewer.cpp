/**
 * @file viewer.cpp
 * Test script for the SDLViewer and Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point list
 *
 * Prints
 * A B
 * where A = number of nodes
 *       B = number of edges
 * and launches an SDLViewer to visualize the system
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "BoundingBox.hpp"

#include "Graph.hpp"


typedef Graph<bool,bool> GraphType;

void remove_box(GraphType& g, const BoundingBox& bb) {
  
  for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter) {

    if (bb.contains((*node_iter).position())) {

      node_iter = g.remove_node(node_iter);
    }
  }
  return;
}




int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  using GraphType = Graph<bool, bool>;
  GraphType graph;
  std::vector<typename GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(2*p - Point(1,1,0)));
    //nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t; 
 /*
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);
  */
  while (CS207::getline_parsed(tets_file, t)) {
    graph.add_edge(nodes[t[0]], nodes[t[1]]);
    graph.add_edge(nodes[t[0]], nodes[t[2]]);
    graph.add_edge(nodes[t[1]], nodes[t[3]]);
    graph.add_edge(nodes[t[2]], nodes[t[3]]);
  }

  double h = graph.edge(0).length();

 
  remove_box(graph, BoundingBox(Point(-0.8+h,-0.8+h,-1), Point(-0.4-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h,-0.8+h,-1), Point( 0.8-h,-0.4-h,1)));
  remove_box(graph, BoundingBox(Point(-0.8+h, 0.4+h,-1), Point(-0.4-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point( 0.4+h, 0.4+h,-1), Point( 0.8-h, 0.8-h,1)));
  remove_box(graph, BoundingBox(Point(-0.6+h,-0.2+h,-1), Point( 0.6-h, 0.2-h,1)));

  // Print number of nodes and edges
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

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


  // Launch a viewer
  CS207::SDLViewer viewer;
  viewer.launch();
  Point source = Point(-1,0,1);
  //int max_degree_of_separation = graph.shortest_path_lengths(graph, source);

  // Set the viewer
  //viewer.draw_graph(graph);
  auto node_map = viewer.empty_node_map(graph);
  
  viewer.add_nodes(graph.node_begin(), graph.node_end(), CS207::HW3Color(), node_map);
  //viewer.add_nodes(graph.node_begin(), graph.node_end(), CS207::DazzlingColor(max_degree_of_separation), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
  viewer.center_view();

  return 0;
}
