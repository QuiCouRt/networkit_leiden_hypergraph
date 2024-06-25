/*
 * HypergraphLeiden.cpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/

#include <networkit/community/HypergraphLeiden.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {
HypergraphLeiden::HypergraphLeiden(const Hypergraph &graph, int iterations, bool randomize, double gamma)
    : CommunityDetectionAlgorithm(graph), gamma(gamma), numberOfIterations(iterations), random(randomize) {
    this->result = Partition(graph.numberOfNodes());
    this->result.allToSingletons();
}

void HypergraphLeiden::run() {
    handler.assureRunning();
    const Hypergraph *currentGraph = HG;

    Partition zeta((*currentGraph).upperNodeIdBound());
    zeta.allToSingletons();
    Partition refined;

    for (int i = 0; i < numberOfIterations; ++i) {
        calculateVolumesHypergraph(*currentGraph);
        parallelMoveHypergraph(*currentGraph);
        flattenPartitionHypergraph();
        refined = parallelRefineHypergraph(*currentGraph);
    }

    result = refined;
    hasRun = true;
}

void HypergraphLeiden::flattenPartitionHypergraph(){

}

void HypergraphLeiden::calculateVolumesHypergraph(const Hypergraph &HG){

}

void HypergraphLeiden::parallelMoveHypergraph(const Hypergraph &HG){

}

Partition HypergraphLeiden::parallelRefineHypergraph(const Hypergraph &HG){
  return Partition(HG.numberOfNodes());
}

} // namespace NetworKit