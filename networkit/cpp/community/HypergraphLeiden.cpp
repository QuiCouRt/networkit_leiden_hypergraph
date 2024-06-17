/*
 * HypergraphLeiden.cpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/


#include <networkit/community/HypergraphLeiden.hpp>

namespace NetworKit {
HypergraphLeiden::HypergraphLeiden(const Graph &graph, int iterations, bool randomize, double gamma)
    : CommunityDetectionAlgorithm(graph), gamma(gamma), numberOfIterations(iterations),
      random(randomize) {
    this->result = Partition(graph.numberOfNodes());
    this->result.allToSingletons();
}


void HypergraphLeiden::run() {
    std::cout << "Thanks for viewing my code!";
}

} // namespace NetworKit