/*
 * HypergraphLeiden.hpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/


#ifndef NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_
#define NETWORKIT_COMMUNITY_HYPERGRAPH_LEIDEN_HPP_

#include <atomic>
#include <condition_variable>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <thread>
#include <tlx/unused.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class HypergraphLeiden final : public CommunityDetectionAlgorithm {
public:
    /**
     *
     * @param graph A networkit hypergraph
     * @param iterations Number of Leiden Iterations to be run
     * @param randomize Randomize node order?
     * @param gamma Resolution parameter
     */
    HypergraphLeiden(const Hypergraph &graph, int iterations = 3, bool randomize = false, double gamma = 1, int type_contribution = 1, double gamma_cut = 0.035);

    void run() override;

    int VECTOR_OVERSIZE = 10000;

private:
    /*inline double modularityDelta(double cutD, double degreeV, double volD) const {
        return cutD - gamma * degreeV * volD * inverseGraphVolume;
    };

    inline double modularityThreshold(double cutC, double volC, double degreeV) const {
        return cutC - gamma * (volC - degreeV) * degreeV * inverseGraphVolume;
    }*/

    static inline void lockLowerFirst(index a, index b, std::vector<std::mutex> &locks) {
        if (a < b) {
            locks[a].lock();
            locks[b].lock();
        } else {
            locks[b].lock();
            locks[a].lock();
        }
    }

    void calculateVolumesHypergraph(const Hypergraph &graph);

    void parallelMoveHypergraph(const Hypergraph &graph, const Partition &zeta);

    double deltaModHypergraph(const Hypergraph &graph, const Partition &zeta, index S, const Partition &p, index c, index target_c);

    Partition parallelRefineHypergraph(const Hypergraph &graph,const Partition &zeta);

    double HypergraphCut(const Hypergraph &graph, const Partition &zeta, index S);

    double GraphVolume; // vol(V)

    std::vector<double> communityVolumes_1; //volume of each community 
    std::vector<double> communityVolumes_2;

    int d; //max size of hyperedge

    std::vector<double> EdgeSizeWeight; // vector of sum of weight of edges of size i (i=0 to d)

    double totalEdgeWeight;

    static constexpr int WORKING_SIZE = 1000;

    double gamma; // Resolution parameter

    double gamma_cut;

    int type_contribution;

    bool changed;

    int numberOfIterations;

    Aux::SignalHandler handler;

    bool random;

    std::vector<std::set<index>> NeighborComm;

    int step=0.0;

};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_PARALLEL_LEIDEN_HPP_