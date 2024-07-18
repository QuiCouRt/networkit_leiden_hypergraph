/*
 * HypergraphLeiden.hpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/


#ifndef NETWORKIT_COMMUNITY_HYPERGRAPH_MODELING_HPP_
#define NETWORKIT_COMMUNITY_HYPERGRAPH_MODELING_HPP_

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

class HypergraphModeling final : public CommunityDetectionAlgorithm {
public:
    /**
     *
     * @param graph A networkit hypergraph
     * @param graph A networkit partition
     */
    HypergraphModeling(const Hypergraph &graph, const Partition &zeta);

    void run() override;

    int VECTOR_OVERSIZE = 10000;

private:

    Aux::SignalHandler handler;

    Partition zeta;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_HYPERGRAPH_MODELING_HPP_