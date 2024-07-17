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
 * $$mod(zeta) := \frac{1}{number of edge} \sum_{C \in zeta} (count^{Hypergraph}(C) - \mathds{E}(count^{\tilde{Hypergraph}}(C))) $$
 * So modularity = coverange - expected coverange
 */
class ModularityHypergraph final : public QualityMeasureHypergraph {

private:
    double gTotalEdgeWeight;

public:
    /** Default constructor */
    ModularityHypergraph();

    /**
     * Returns the Modularity of the given clustering with respect to the hypergraph @a G.
     *
     * @param zeta The clustering.
     * @param G The hypergraph.
     * @param gamma Resolution parameter
     * @param type_contribution The type of modularity : 00=strict , or 01=majority , or 10=strict/partially weighted , or 11=majority/partially weighted.
     * @return The modularity.
     */
    double getQualityHypergraph(const Partition &zeta, const Hypergraph &G, double gamma = 1.0, int type_contribution=1);

    /**
     * Returns the gain of Modularity obtained by moving the set of nodes S in community c.
     * Warning : all nodes of S muss belong to the same community. 
     *
     * @param zeta The clustering.
     * @param G The hypergraph.
     * @param S The set of nodes that we want to move.
     * @param c The target community.
     * @param gamma Resolution parameter
     * @param type_contribution The type of modularity : 00=strict , or 01=majority , or 10=strict/partially weighted , or 11=majority/partially weighted.
     * @return The modularity gain.
     */
    double deltaModularityHypergraph(const Partition &zeta, const Hypergraph &G, std::set<node> S, index c, double gamma=1.0, int type_contribution=1);

    /**
     * @param totalEdgeWeight Sum of all edge weights in @a G. If specified, it does not
     *        have to be computed.
     */
    void setTotalEdgeWeight(double totalEdgeWeight);
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_MODULARITY_HPP_