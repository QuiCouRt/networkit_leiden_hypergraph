/*
 * ModularityHypergraph.hpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/



#ifndef NETWORKIT_COMMUNITY_MODULARITY_HYPERGRAPH_HPP_
#define NETWORKIT_COMMUNITY_MODULARITY_HYPERGRAPH_HPP_

#include <networkit/community/QualityMeasureHypergraph.hpp>
#include <networkit/graph/Hypergraph.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Modularity is a quality index for community detection. It assigns a quality value in [-0.5, 1.0]
 * to a partition of a graph which is higher for more modular networks and partitions which better
 * capture the modular structure.
 *
 * Modularity is defined as:
 *
 * $$mod(\\zeta) := \\frac{\\sum_{C \\in \\zeta} \\sum_{ e \\in E(C) } \\omega(e)}{\\sum_{e \\in E}
 * \\omega(e)}
 * \\frac{ \\sum_{C \\in \\zeta}( \\sum_{v \\in C} \\omega(v) )^2 }{4( \\sum_{e \\in E} \\omega(e)
 * )^2 }$$
 */
class ModularityHypergraph final : public QualityMeasureHypergraph {

private:
    double gTotalEdgeWeight;

public:
    /** Default constructor */
    ModularityHypergraph();

    /**
     * Returns the Modularity of the given clustering with respect to the graph @a G.
     *
     * @param zeta The clustering.
     * @param G The graph.
     * @return The modularity.
     */
    double getQualityHypergraph(const Partition &zeta, const Hypergraph &G, int type_contribution);

    /**
     * Returns the gain of Modularity obtained by moving the set of nodes S in community c, 
     *
     * @param zeta The clustering.
     * @param G The graph.
     * @param S The set of nodes that we want to move.
     * @param c The target community.
     * @param type_contribution The type of modularity : 0=strict or 1=majority.
     * @return The modularity gain.
     */
    double deltaModularityHypergraph(const Partition &zeta, const Hypergraph &G, std::vector<node> S, index c, int type_contribution);

    /**
     * @param totalEdgeWeight Sum of all edge weights in @a G. If specified, it does not
     *        have to be computed.
     */
    void setTotalEdgeWeight(double totalEdgeWeight);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_MODULARITY_HPP_