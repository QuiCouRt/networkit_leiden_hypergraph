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

// TODO check if binom exist in networkit
double BinomialCoefficient(const double n, const double k) {
  std::vector<double> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (double i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

namespace NetworKit {

ModularityHypergraph::ModularityHypergraph() : QualityMeasureHypergraph(), gTotalEdgeWeight(0.0) {}

void ModularityHypergraph::setTotalEdgeWeight(double totalEdgeWeight) {
    gTotalEdgeWeight = totalEdgeWeight;
}

double ModularityHypergraph::getQualityHypergraph(const Partition &zeta, const Hypergraph &G, int type_contribution) {
    assert(G.numberOfNodes() <= zeta.numberOfElements());
    double cov = 0.0;
    double expCov = 0.0;
    double modularity =0.0; // mod = coverage - expected coverage

    //For loop to obtain the max size edge
    double d= 0;
    for (edgeid eid = 0; eid < G.upperEdgeIdBound(); ++eid) {
        if (G.getEdgeIncidence(eid).size()>d){
            d=G.getEdgeIncidence(eid).size();
        }
    }

    if (type_contribution==0) {// >>>>> STRICT EDGE CONTRIBUTION

        // >>> Firstly, we compute the coverage : 
        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add all edge weights
        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of weight of edges of size i (i=0 to d)

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
            // an edge belongs to a community if all nodes in this edge belong to the same community 
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

        // >>> Secondly, we compute the expected coverage : 
        std::set<std::set<index>> subsets =zeta.getSubsets(); //obtain all partition 
        std::vector<double> intraVolume(zeta.upperBound(), 0.0);
        double volume=0.0;

        // Compute volume of each community 
        double i=0;
        for (std::set c : subsets){
            for (index n : c){
                for (edgeid eid: G.getNodeIncidence(n)){
                    intraVolume[i]+=G.getEdgeWeight(eid);
                }
            }
            i++;
        }

        //Compute volume of all hypergraph
        G.forNodes([&](node nid, nodeweight nw) {
            for (edgeid eid: G.getNodeIncidence(nid)){
                volume+=G.getEdgeWeight(eid);
            }
        });

        double sum_vol =0.0;
        // Compute expCov : \sum_{d} (\frac{|E_d|}{vol(V)^d} \sum_{c} vol(c)^d)
        for (double j=0 ; j < d+1; j++){
            sum_vol=0.0;
            for (i=0; i < zeta.upperBound(); i++){
                sum_vol += pow(intraVolume[i],j);
            }
            expCov += (EdgeSizeWeight[j] / pow(volume,j)) * sum_vol;
        }
        
        expCov = expCov / totalEdgeWeight; 
        modularity = cov - expCov;

    }// End if for strict edge contribution : type_contribution == 0 


    if (type_contribution==1) {// >>>>> MAJORITY EDGE CONTRIBUTION
        
        //>>> Firstly, we compute the coverage : 

        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add all edge weights

        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of weight of edges of size i (i=0 to d)

        // compute intra-cluster edge weights per cluster
        bool isFirst=true;
        bool edgeBelongs=true;
        index c, c_prime;
        for (edgeid eid = 0; eid < G.upperEdgeIdBound(); ++eid) {
            isFirst=true;
            edgeBelongs=true;
            // an edge belongs to a community if the majority of its nodes belongs to the same community 
            for (auto nid: G.getEdgeIncidence(eid)) {
                // TODO max occurrence  >? d/2 = G.getEdgeIncidence(eid).size() / 2 (diision entiere si tu veux)
            }
            double w = G.getEdgeWeight(eid);
            if (edgeBelongs){
                intraEdgeWeightSum += w;
            }
            EdgeSizeWeight[G.getEdgeIncidence(eid).size()] += w;
            totalEdgeWeight+= w;
        }

        cov = intraEdgeWeightSum / totalEdgeWeight;


        //>>> Secondly, we compute the expected coverage : 

        std::set<std::set<index>> subsets =zeta.getSubsets(); //obtain all partition 
        std::vector<double> intraVolume(zeta.upperBound(), 0.0);
        double volume=0.0;

        // Compute volume of each community 
        double i=0;
        for (std::set c : subsets){
            for (index n : c){
                for (edgeid eid: G.getNodeIncidence(n)){
                    intraVolume[i]+=G.getEdgeWeight(eid);
                }
            }
            i++;
        }

        //Compute volume of all hypergraph
        G.forNodes([&](node nid, nodeweight nw) {
            for (edgeid eid: G.getNodeIncidence(nid)){
                volume+=G.getEdgeWeight(eid);
            }
        });

        double sum_vol =0.0;
        double sum =0.0;
        // Compute expCov : \sum_d \frac{|E_d|}{vol(V)^d} \sum_{i= \frac{d}{2} +1}^{d}  \binom{d}{i} \sum_{c \in C} (vol(c))^i (vol(V)- vol(c))^{d-i}
        for (double j=0 ; j < d+1; j++){
            sum=0.0;
            for (i= (j / 2) + 1; i <= j; i++){
                sum_vol=0.0;
                //TODO !!!
                for (double c=0; c < zeta.upperBound(); c++){
                    sum_vol += pow(intraVolume[i],j);
                }
                sum_vol= BinomialCoefficient(j,i)* sum_vol;
                sum += sum_vol;
            }
            expCov += (EdgeSizeWeight[j] / pow(volume,j)) * sum;
        }

        expCov = expCov / totalEdgeWeight; 
        
        modularity = cov - expCov;
    }
    return modularity;
}

} /* namespace NetworKit */