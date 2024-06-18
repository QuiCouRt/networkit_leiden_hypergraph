/*
 * ModularityHypergraph.cpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/


#include <cmath>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Coverage.hpp>
#include <networkit/graph/Hypergraph.hpp>
#include <networkit/community/ModularityHypergraph.hpp>

namespace NetworKit {

ModularityHypergraph::ModularityHypergraph() : QualityMeasureHypergraph(), gTotalEdgeWeight(0.0) {}

void ModularityHypergraph::setTotalEdgeWeight(double totalEdgeWeight) {
    gTotalEdgeWeight = totalEdgeWeight;
}

double ModularityHypergraph::getQualityHypergraph(const Partition &zeta, const Hypergraph &G, int type_contribution) {
    assert(G.numberOfNodes() <= zeta.numberOfElements());
    double cov = 0.0;
    if (type_contribution==0) {// >>> STRICT EDGE CONTRIBUTION
  
        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add edge weigh

        //For loop to obtain the max size edge
        int d= 0;
        for (edgeid eid = 0; eid < G.upperEdgeIdBound(); ++eid) {
            if (G.getEdgeIncidence(eid).size()>d){
                d=G.getEdgeIncidence(eid).size();
            }
        }

        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of edges of size i (i=0 to d)


#ifdef DEBUG
        if (!G.getEdgeIncidence(0).empty()) {
            ERROR("First edge error");
        }
#endif


        // compute intra-cluster edge weights per cluster
        bool isFirst=true;
        bool edgeBelongs=true;
        index c, c_prime;
        for (edgeid eid = 0; eid < G.upperEdgeIdBound(); ++eid) {
            isFirst=true;
            edgeBelongs=true;

            for (auto nid: G.getEdgeIncidence(eid)) {
                if(isFirst){
                    c = zeta[nid];
                    isFirst=false;
                }
                c_prime = zeta[nid];
                if (c!=c_prime){
                    edgeBelongs=false;
                    break;
                }
            }
            double w = G.getEdgeWeight(eid);
            if (edgeBelongs){
                intraEdgeWeightSum += w;
            }
            EdgeSizeWeight[G.getEdgeIncidence(eid).size()] += w;
            totalEdgeWeight+= w;
        }

        cov = intraEdgeWeightSum / totalEdgeWeight;
        cov = EdgeSizeWeight[2];

/*
    Coverage coverage;
    // deprecated: intraEdgeWeightSum / gTotalEdgeWeight;
    double cov = coverage.getQuality(zeta, G);
    //// term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C}
    ///\omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
    double expCov;
    double modularity; // mod = coverage - expected coverage
    if (gTotalEdgeWeight == 0.0) {
        gTotalEdgeWeight = G.totalEdgeWeight(); // compute total edge weight in G
        DEBUG("computed total edge weight: ", gTotalEdgeWeight);
    }

    if (gTotalEdgeWeight == 0.0) {
        ERROR("G: m=", G.numberOfEdges(), "n=", G.numberOfNodes());
        throw std::invalid_argument(
            "Modularity is undefined for graphs without edges (including self-loops).");
    }

    //!< cluster -> sum of the weights of incident edges for all nodes
    std::vector<double> incidentWeightSum(zeta.upperBound(), 0.0);

    // compute volume of each cluster
    G.parallelForNodes([&](node v) {
        // add to cluster weight
        index c = zeta[v];
        assert(zeta.lowerBound() <= c);
        assert(c < zeta.upperBound());
        // account for self-loops a second time
#pragma omp atomic
        incidentWeightSum[c] += G.weightedDegree(v) + G.weight(v, v);
    });

    // compute sum of squared cluster volumes and divide by squared graph volume
    expCov = 0.0;

#pragma omp parallel for reduction(+ : expCov)
    for (omp_index c = static_cast<omp_index>(zeta.lowe/isaline/Bureau/testFolder/Stage_2024_Hypergraph/code/networkit/include/networkit/community/QualityMeasure.hpp:23:20: note:     ‘virtual double NetworKit::QualityMeasure::getQuality(const NetworKit::Partition&, const NetworKit::Graph&)’
[build]    23 |     virtual double getQuality(const Partition &zeta, const Graph &G) = 0rBound());
         c < static_cast<omp_index>(zeta.upperBound()); ++c) {
        // squared
        expCov +=
            ((incidentWeightSum[c] / gTotalEdgeWeight) * (incidentWeightSum[c] / gTotalEdgeWeight))
            / 4;
    }

    DEBUG("expected coverage: ", expCov);

    // assert ranges of coverage
    assert(cov <= 1.0);
    assert(cov >= 0.0);
    assert(expCov <= 1.0);
    assert(expCov >= 0.0);

    modularity = cov - expCov;
    DEBUG("modularity = ", modularity);

    // reset totalEdgeWeight
    gTotalEdgeWeight = 0.0;

    assert(!std::isnan(modularity)); // do not return NaN
    // do not return anything not in the range of modularity values
    assert(modularity >= -0.5);
    assert(modularity <= 1);
    return modularity;
    */
}
    return cov;
}// End if for strict edge contribution : type_contribution == 0 

} /* namespace NetworKit */