/*
 * HypergraphConverterForModeling.cpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/
#include <fstream>
#include <networkit/community/HypergraphConverterForModeling.hpp>
#include <networkit/graph/Hypergraph.hpp>


namespace NetworKit {

HypergraphModeling::HypergraphModeling(const Hypergraph &graph, const Partition &zeta)
    : CommunityDetectionAlgorithm(graph), zeta(zeta) {
    this->result = Partition(graph.numberOfNodes());
    this->result.allToSingletons();
}

void HypergraphModeling::run() {

    handler.assureRunning();
    const Hypergraph *currentGraph = HG;

    FILE * f;
    f = fopen ("output/Hypergraph.txt", "wt");
    if (f == NULL)
        DEBUG("Unable to open file for writing!");
    else
    {
    //fprintf (f, "%d", currentGraph->numEdges);
    (*currentGraph).forEdges([&](edgeid eid, edgeweight ew) {
        for (node nid: (*currentGraph).getEdgeIncidence(eid)){
            fprintf (f, "%ld ", nid);
        }
        fprintf (f, "\n");
    });
    fclose (f);

    hasRun = true;
    }

    FILE * p;
    p = fopen ("output/Partition.txt", "wt");
    if (p == NULL)
        DEBUG("Unable to open file for writing!");
    else
    {
    //fprintf (f, "%d", currentGraph->numEdges);
    (*currentGraph).forNodes([&](node nid, nodeweight nw) {
      fprintf (p, "%ld\n", zeta[nid]); 
    });

    fclose (p);

    hasRun = true;
    }
}


} // namespace NetworKit