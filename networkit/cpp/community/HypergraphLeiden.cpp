/*
 * HypergraphLeiden.cpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/

#include <networkit/community/HypergraphLeiden.hpp>
#include <networkit/graph/Hypergraph.hpp>

double BinomialCoef(const double n, const double k) {
  std::vector<double> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (double i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}


namespace NetworKit {
HypergraphLeiden::HypergraphLeiden(const Hypergraph &graph, int iterations, bool randomize, double gamma, int type_contribution)
    : CommunityDetectionAlgorithm(graph), gamma(gamma), numberOfIterations(iterations), random(randomize), type_contribution(type_contribution) {
    this->result = Partition(graph.numberOfNodes());
    this->result.allToSingletons();
}

void HypergraphLeiden::run() {
    handler.assureRunning();
    const Hypergraph *currentGraph = HG;
    
    //For loop to obtain the max size edge
    d= 0;
    (*currentGraph).forEdges([&](edgeid eid, edgeweight ew) {
        if ((*currentGraph).getEdgeIncidence(eid).size()>d){
            d=(*currentGraph).getEdgeIncidence(eid).size();
        }
    });

    // vector of sum of weight of edges of size i (i=0 to d)
    EdgeSizeWeight.resize(d+1 + VECTOR_OVERSIZE);
    totalEdgeWeight=0.0;
    (*currentGraph).forEdges([&](edgeid eid, edgeweight ew) {
        EdgeSizeWeight[(*currentGraph).getEdgeIncidence(eid).size()] += ew;
        totalEdgeWeight += ew;
    });

    // at initialization each node is in a singleton community
    Partition zeta((*currentGraph).upperNodeIdBound());
    zeta.allToSingletons();
    // For each community, we obtain the set of neighboring communities
    NeighborComm.resize((*currentGraph).upperNodeIdBound());
    (*currentGraph).forNodes([&](node nid, nodeweight nw) {
      index c = zeta[nid];
      for (edgeid eid: (*currentGraph).getNodeIncidence(nid)){
        for (node nid_prime: (*currentGraph).getEdgeIncidence(eid)){
          NeighborComm[nid].insert(nid_prime);
        }
      }
    });

    Partition refined;

    for (int i = 0; i < numberOfIterations; ++i) {
        calculateVolumesHypergraph(*currentGraph);
        parallelMoveHypergraph((*currentGraph), zeta);
        refined = parallelRefineHypergraph(*currentGraph);
    }

    result = refined;
    hasRun = true;
}


void HypergraphLeiden::calculateVolumesHypergraph(const Hypergraph &graph){
  auto timer = Aux::Timer();
  timer.start();
  communityVolumes_1.clear();
  communityVolumes_1.resize(result.upperBound() + VECTOR_OVERSIZE);
  //if (graph.isWeighted()) {
  if (true) {
    std::vector<double> threadVolumes(omp_get_max_threads());
    graph.parallelForNodes([&](node nid) {
      {
      for (edgeid eid: graph.getNodeIncidence(nid)){
        edgeweight ew =graph.getEdgeWeight(eid);
#pragma omp atomic
        communityVolumes_1[result[nid]] += ew;
        threadVolumes[omp_get_thread_num()] += ew;
      }}
    });
    for (const auto vol : threadVolumes) {
      GraphVolume += vol;
    }
    GraphVolume = GraphVolume;
  } 
  TRACE("Calculating Volumes took " + timer.elapsedTag());
}

double HypergraphLeiden::deltaModHypergraph(const Hypergraph &graph, const Partition &zeta, index S, index c, index target_c){
  double modularityGain =0.0; // modularityGain = covGain - expCovGain
  double covGain =0.0;
  double expCovGain = 0.0;

  //The volume of S
  double vol_S=communityVolumes_1[S];
  double vol_c = communityVolumes_2[c];
  double vol_target_c = communityVolumes_2[target_c];

  if (type_contribution==0) {// >>>>> STRICT EDGE CONTRIBUTION
    //covGain = \frac{In_w(S,c')-In_w(S,c)}{|E|_w} 
    //with $In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq A ~ \wedge ~ e \nsubseteq B ~ \wedge ~ e \subseteq A \cup B} w_e$
    double in_c_S = 0.0;
    double in_target_c_S = 0.0;
    bool edgeBelongs_c=true;
    bool edgeBelongs_target_c=true;
    std::set<node> set_S = zeta.getMembers(S);
    std::set<edgeid> border_edge;
    for (node n : set_S){
      for (edgeid eid: graph.getNodeIncidence(n)){
        border_edge.insert(eid);
      }
    }

    for (edgeid eid: border_edge){
      edgeBelongs_c =false;
      for (node nid: graph.getEdgeIncidence(eid)){
        if (!edgeBelongs_c && zeta[nid]!=S){ // eid is not strictly included in S, so it's possible to take eid into account
          edgeBelongs_c = true;
        }
        if (result[nid] !=c){
          edgeBelongs_c = false;
          break; 
        }
      }
      if (edgeBelongs_c){
        in_c_S +=graph.getEdgeWeight(eid);
      }
                
      edgeBelongs_target_c=false;
      for (node nid: graph.getEdgeIncidence(eid)){
        if (!edgeBelongs_target_c && zeta[nid]!=S){ // eid is not strictly included in S, so it's possible to take eid into account
          edgeBelongs_target_c = true;
        }
        if (result[nid] !=target_c && zeta[nid]!=S){
          edgeBelongs_target_c = false;
          break; 
        }
      }
      if (edgeBelongs_target_c){
        in_target_c_S += graph.getEdgeWeight(eid);
      }
    }
    covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;

    //expCovGain = \sum_{d\geq 2} \frac{|E_d|_w}{|E|_w vol_w(V)^d}(vol_w(c'+S)^d -vol_w(c')^d -vol_w(c)^d + vol_w(c-S)^d)
    for (double j=0 ; j < d+1; j++){
      expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(GraphVolume,j))) * ( pow(vol_target_c + vol_S,j)-pow(vol_target_c,j)  - pow(vol_c,j) + pow(vol_c - vol_S,j));
    }
  }


  if (type_contribution==1) {// >>>>> MAJORITY EDGE CONTRIBUTION
    // covGain = \frac{In(S,c')-In(S,c)}{|E|} 
    // with In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq_{maj} A ~ \wedge ~ e \nsubseteq_{maj} B ~ \wedge ~ e \subseteq_{maj} A \cup B} w_e
    double in_c_S = 0.0;
    double in_target_c_S = 0.0;
    bool edgeBelongs_c=true;
    bool edgeBelongs_target_c=true;
    std::set<node> set_S = zeta.getMembers(S);
    std::set<edgeid> border_edge;
    for (node n : set_S){
      for (edgeid eid: graph.getNodeIncidence(n)){
        border_edge.insert(eid);
      }
    }

    for (edgeid eid: border_edge){
      edgeBelongs_c =false;
      edgeBelongs_target_c=false;
      std::vector<double> CommEdge(3, 0.0);
      for (node nid: graph.getEdgeIncidence(eid)){
        if (result[nid]==c){CommEdge[0]++;}
        if (result[nid]==target_c){CommEdge[1]++;}
        if (zeta[nid]!=S){CommEdge[2]++;}
      }
                
      if (!(CommEdge[1] >= int(graph.getEdgeIncidence(eid).size()/ 2. +1)) && !(CommEdge[0]- CommEdge[2] >= int(graph.getEdgeIncidence(eid).size()/ 2. +1)) ){
        if (CommEdge[0] >= int(graph.getEdgeIncidence(eid).size()/ 2. +1)){
          edgeBelongs_c =true;
        }
        if (CommEdge[1] + CommEdge[2] >= int(graph.getEdgeIncidence(eid).size()/ 2. +1)){
          edgeBelongs_target_c =true;
        }
      }
            
      if (edgeBelongs_c){
        in_c_S += graph.getEdgeWeight(eid);
      }

      if (edgeBelongs_target_c){
        in_target_c_S += graph.getEdgeWeight(eid);
      }
    }
    covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;


    // expCovGain = \sum_{d\geq 2} \frac{|E_d|}{|E| vol(V)^d} \sum_{i=\frac{d}{2}+1}^{d} \binom{d}{i} (vol(c'+S)^i(vol(V)-vol(c'+S))^{d-i} -vol(c')^i(vol(V)-vol(c'))^{d-i} -vol(c)^i(vol(V)-vol(c))^{d-i} + vol(c-S)^i(vol(V)-vol(c-S))^{d-i})
    double sum = 0.0;
    for (double j=0 ; j < d+1; j++){
      sum = 0.0;
      for (double i= int ((j / 2.) + 1); i <= j; i++){
        sum += BinomialCoef(j,i) *(pow(vol_target_c + vol_S,i)*pow(GraphVolume - vol_target_c - vol_S,j-i) -  pow(vol_target_c,i)*pow(GraphVolume - vol_target_c,j-i) - pow(vol_c ,i)*pow(GraphVolume- vol_c,j-i) + pow(vol_c - vol_S,i)*pow(GraphVolume - vol_c + vol_S,j-i));
      }
      expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(GraphVolume,j))) * sum;
    }
  } 

  modularityGain= covGain - expCovGain; 

  return modularityGain;
}

void HypergraphLeiden::parallelMoveHypergraph(const Hypergraph &graph, const Partition &zeta){
  //TODO
  //double maxDelta = std::numeric_limits<double>::lowest();
  //DEBUG("Local Moving : ", graph.numberOfNodes(), " Nodes ");

  //The zeta partition corresponds to our block of nodes that we will move in this phase of greedy move
  //Initialy zeta = result 
  //After we transform result
  count moved = 0;
  count totalNodes = 0;
  int singleton = 0;

  std::vector<bool> inQueue(zeta.upperBound(), false); //Boolean array allowing you not to add the same node twice in the queue
  std::queue<index> queue; //nodes that remain to be processed
  
  communityVolumes_2.resize(communityVolumes_1.size());
  copy(communityVolumes_1.begin(), communityVolumes_1.end(), communityVolumes_2.begin());
  bool resize = false;
  int waitingForResize = 0;
  int waitingForNodes = 0;

  uint64_t vectorSize = communityVolumes_1.capacity();
  int upperBound = result.upperBound();

  std::vector<index> currentBlocks; // block to cover during a loop
  std::vector<index> newNodes;
  newNodes.reserve(WORKING_SIZE);
  //std::vector<double> cutWeights(communityVolumes.capacity());
  std::vector<index> pointers;

  // We add all node in queue 
  std::unique_ptr<std::atomic<bool>[]> exists(new std::atomic<bool>[zeta.upperBound()] {});
  zeta.parallelForEntries([&](index, index s) {
    if (s != none) {
      exists[s] = true;
    }
  });
  count k = 0; // number of actually existing clusters
//#pragma omp parallel for reduction(+ : k)
  for (index i = 0; i < static_cast<index>(zeta.upperBound()); ++i) {
    if (exists[i]) {
      k++;
      currentBlocks.push_back(i);
      inQueue[i] = true;
    }
  }

  // We randomize the order of the blocks
  if (random) {
    auto &mt = Aux::Random::getURNG();
    std::shuffle(currentBlocks.begin(), currentBlocks.end(), mt);
  }

  do {
    handler.assureRunning();
    for (index u : currentBlocks) {
      // Resize if necessary
      if (resize) {
        waitingForResize++;
        while (resize) {
          std::this_thread::yield();
        }
        waitingForResize--;
      }
      //cutWeights.resize(vectorSize);

      //Check u not empty
      std::set<node> S= zeta.getMembers(u);
      if (S.empty()){
        continue;
      }
      index current_comm = result[*(S.begin())];
      assert(inQueue[u]); // index in currentBlocks must be in the queue
      double maxDelta = std::numeric_limits<double>::lowest();
      index bestCommunity = none;
      //double degree = 0;
        //for (auto z : pointers) {
        //    cutWeights[z] = 0;
        //}
      pointers.clear();

      for (index neighborCommunity : NeighborComm[u]){
        index target_comm = result[(zeta.giveOne(neighborCommunity))];
        if (target_comm != current_comm) {
          double delta = deltaModHypergraph(graph, result, u, target_comm, neighborCommunity);
          if (delta > maxDelta) {
            maxDelta = delta;
            bestCommunity = target_comm;
          }
        }
      }
        
      moved++;
      // On laisse la communauté dans un singleton 
      /*if (maxDelta < 0) {
        singleton++;
        bestCommunity = upperBound++;
        if (bestCommunity >= communityVolumes_1.capacity()) {
          resize = true;
          vectorSize += VECTOR_OVERSIZE;
          communityVolumes_1.resize(vectorSize);
          resize = false;
        }
      }*/
      // maxDelta > 0 on 
      for (node nid : zeta.getMembers(u)){
        result[nid] = bestCommunity;
      }
      communityVolumes_2[bestCommunity] += communityVolumes_1[u];
      communityVolumes_2[current_comm]-= communityVolumes_1[u];
      changed = true;
      inQueue[u] = false;

      for (index neighborCommunity : NeighborComm[u]){
        index target_comm = result[(zeta.giveOne(neighborCommunity))];
        if (result[target_comm] != bestCommunity &&  neighborCommunity !=u) {
          if (!inQueue[neighborCommunity]) {
            newNodes.push_back(neighborCommunity);
            inQueue[neighborCommunity] = true;
            if (newNodes.size() == WORKING_SIZE) {
              for (node v : newNodes) {
                queue.push(v);
              }
              newNodes.clear();
              newNodes.reserve(WORKING_SIZE);
            }
          }
        }
      }  
    }

    totalNodes += currentBlocks.size();
    if (!newNodes.empty()) {
      currentBlocks = std::move(newNodes);
      newNodes.clear();
    } else if (!queue.empty()) {
      currentBlocks.clear();
      while (!queue.empty()) {
        currentBlocks.push_back(queue.front());
        queue.pop();
      }
    } else {
      break;
    }
    } while (true);

    result.setUpperBound(upperBound);
    if (Aux::Log::isLogLevelEnabled(Aux::Log::LogLevel::DEBUG)) {
        DEBUG("Total worked: ", totalNodes, " Total moved: ", moved, " moved to singleton community: ", singleton);
    }
  
}


Partition HypergraphLeiden::parallelRefineHypergraph(const Hypergraph &HG){
  //TODO
  return Partition(HG.numberOfNodes());
}

} // namespace NetworKit