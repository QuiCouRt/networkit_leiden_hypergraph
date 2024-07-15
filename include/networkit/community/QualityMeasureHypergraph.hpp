/*
 * ModularityHypergraph.hpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/

#ifndef NETWORKIT_COMMUNITY_QUALITY_HYPERGRAPH_MEASURE_HPP_
#define NETWORKIT_COMMUNITY_QUALITY_HYPERGRAPH_MEASURE_HPP_

#include <networkit/graph/Hypergraph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for all clustering quality measures.
 */
class QualityMeasureHypergraph {

public:
    virtual double getQualityHypergraph(const Partition &zeta, const Hypergraph &G, double gamma, int type_contribution) = 0;
    virtual double deltaModularityHypergraph(const Partition &zeta, const Hypergraph &G, std::set<node> S, index c, double gamma, int type_contribution) = 0;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_QUALITY_MEASURE_HPP_
