# Isaline's Codes

Several files have been implemented for the execution of the Leiden algorithm on hypergraphs.

# Role of different files



## ModularityHypergraph.cpp

This file contains modularity calculation functions. This file has two functions:

```cpp

double getQualityHypergraph(const Partition &zeta, const Hypergraph &G, double gamma = 1.0, int type_contribution=1);
```

This function returns the Modularity of the given clustering `zeta` with respect to the hypergraph `G`. Parameters are 
- `zeta` the clustering
- `G` our hypergraph
- `gamma` the resolution parameter: 1 (standard modularity), 0 (one community), 2 * numberOfEdge (singleton communities)
- `type_contribution` the type of modularity : 00=strict , or 01=majority , or 10=strict and partially weighted , or 11=majority and partially weighted.

```cpp

double deltaModularityHypergraph(const Partition &zeta, const Hypergraph &G, std::set<node> S, index c, double gamma=1.0, int type_contribution=1);
```

This function returns the gain of Modularity obtained by moving the set of nodes `S` in community `c`. Warning : all nodes of `S` muss belong to the same community.

Parameters are 
- `zeta` the clustering
- `G` our hypergraph
- `S` the set of nodes that wae xant to move
- `c` the target community
- `gamma` the resolution parameter: 1 (standard modularity), 0 (one community), 2 * numberOfEdge (singleton communities)
- `type_contribution` the type of modularity : 00=strict , or 01=majority , or 10=strict and partially weighted , or 11=majority and partially weighted.

-> **As there are 4 types of modularity, each function is in 4 parts managed with if conditions on the `type_contribution`.**

We can make the code more beautiful with fewer conditions. Indeed, the **modularity** is **equal** to the coverage of the graph by the clustering `cov` **minus** expectation of the coverage `expCov`. 
But the `cov` does not change between strict modularity and partially weighted strict modularity. The same for majority modularity. There would be a way to make the code more readable by making a less crude if condition. 
This was done in calculating the modularity in the file **HypergraphLeiden.cpp**.



## HypergraphLeiden.cpp

This file contains the implementation of the Leiden algorithm. This file has six functions. The main function is :

```cpp

HypergraphLeiden(const Hypergraph &graph, int iterations = 3, bool randomize = false, double gamma = 1, double gamma_cut = none, int type_contribution = 1, double (*weightFun)(double)= &weight_1);
```

This function returns a clustering `result` of the hypergraph `graph` given by the Leiden method. Parameters are :
- `graph` our hypergraph
- `iterations` the maximal number of Leiden iterations (the algorithm may stop before)
- `randomize` boolean to choose if the node traversal order is random
- `gamma` the resolution parameter: 1 (standard modularity), 0 (one community), 2 * numberOfEdge (singleton communities)
- `gamma_cut` the resolution parameter for refinement phase (to determine if a node is well-connected with its community), by default = 1 / vol(`graph`)
- `type_contribution` the type of modularity : 00=strict , or 01=majority , or 10=strict and partially weighted , or 11=majority and partially weighted.
- `weightFun` the weight function chosen to weight the edges. This function calculates the weight of an edge based on its order/size.

-> **I had only passed a function as an argument to another function in C++. The used method for weightFun is probably not the right one.**

The auxiliary functions are : 

```cpp

void calculateVolumesHypergraph(const Hypergraph &graph);
```

This function initializes the volume of each community in `communityVolumes_1`. These communities are defined by the partition `result`. Note: this is the only function where I tried to do a little parallelization before abandoning this idea to concentrate on having a Leiden algorithm that works. 

Note : there are `communityVolumes_1` and `communityVolumes_1_unweighted` because we need both types of volume depending on the modularity chosen (totally or partially weighted). 

Note :  there are `communityVolumes_1` and `communityVolumes_2`.
- `communityVolumes_1` corresponds to the volume of our current nodes/blocks.
- `communityVolumes_2` corresponds to the volume of the communities we are building (by greedy move of blocks)


```cpp

void MoveHypergraph(const Hypergraph &graph, const Partition &zeta);
```

This function does the greedy move phase. 
In this function, `zeta` represents our blocks obtained in the previous step which we will consider as indivisible nodes/blocks. These blocks will be moved according to the greedy move method to give a larger community (saved in `result` partition). 


```cpp

double deltaModHypergraph(const Hypergraph &graph, const Partition &zeta, index S, const Partition &p, index c, index target_c);
```

This function computates the gain of modularity of moving a set of nodes `S` from the community `c` to the community `target_c`.
- `zeta` is the partition of our blocks that we manipulate and of which `S` is a part. 
- `p` is the community partition that we are currently modifying and of which `c` and `target_c` are a part. 


```cpp

Partition RefineHypergraph(const Hypergraph &graph,const Partition &zeta);
```

This function does the refinement phase. 
In this function we have three different partitions:
- `zeta` (that didn't change during the greedy move phase) which is the partition of our nodes/blocks that we manipulate. The volumes are saved in `communityVolumes_1`.
- `result` which corresponds to the big community obtained during the greedy move phase. It is in these big communities that we will once again associate our blocks together except the no well-connected ones. The volumes are saved in `bigCommunityVolumes`.
- `refined` which is the partition that we are building and which initially is equal to `zeta`, each block is in a singleton community. The volumes are saved in `communityVolumes_2`.

```cpp

double HypergraphCut(const Hypergraph &graph, const Partition &zeta, index S);
```

This function computes the cut between a set (given by `S` and `zeta`) and the big community to which it belongs (given by `result`). 



## HypergraphLouvain.cpp

This file contains the implementation of the Louvain algorithm. This file has four functions. 

-> **Overall this file is a copy/paste of the file HypergraphLeiden.cpp**. Only the tools and functions managing the refinement phase were removed. 

-> So if an optimization is done in the file HypergraphLeiden.cpp, it must then also be done in the file HypergraphLouvain.cpp. This was prevented by manipulating subclasses but I wanted to code quickly.

 The main function is :

```cpp

HypergraphLouvain(const Hypergraph &graph, int iterations = 3, bool randomize = false, double gamma = 1, int type_contribution = 1, double (*weightFun)(double)= &weight_1bis);
```

This function returns a clustering `result` of the hypergraph `graph` given by the Louvain method. Parameters are :
- `graph` our hypergraph
- `iterations` the maximal number of Leiden iterations (the algorithm may stop before)
- `randomize` boolean to choose if the node traversal order is random
- `gamma` the resolution parameter: 1 (standard modularity), 0 (one community), 2 * numberOfEdge (singleton communities)
- `type_contribution` the type of modularity : 00=strict , or 01=majority , or 10=strict and partially weighted , or 11=majority and partially weighted.
- `weightFun` the weight function chosen to weight the edges. This function calculates the weight of an edge based on its order/size.




## HypergraphConverterForModeling.cpp

This file contains a function which allows us to write our graph and a partition of this graph in a `.txt` file which can simply be read to have our graph in **hypernetx format**.

The function is : 
```cpp

HypergraphModeling(const Hypergraph &graph, const Partition &zeta);
```

Parameters are :
- `graph` our hypergraph
- `zeta` a clustering of this graph



-> I am completely aware that this file is ugly and should not be kept. But I chose the simplest and quickest solution to have a modeling.



## Hypergraph_visualization.ipynb
This file read `.txt` files to give a modeling of the hypergraph with **hypernetx**.

### Example of use : 

It is enough to have in our code the following lines for example if we want to see the result of Leiden on our `hypergraph`. 

```cpp
    HypergraphLeiden pl(hypergraph);
    pl.run();
    Partition zeta = pl.getPartition();

    HypergraphModeling m(hypergraph,zeta);
    m.run();
```

This will create the `.txt` files in networkit/output.
Then is enough to run the file **Hypergraph_visualization.ipynb** to obtain a modeling of our hypergraph and these communities. 



# Ideas for Leiden code optimization

Here are some ideas to optimize the implementation of the Leiden algorithm. We work only on HypergraphLeiden.cpp file.



## Parallelization 

Parallelization is simple to set up. Just follow the same method as Leiden in parallel for the graphs. 



## Neighborhood

In Leiden, we manipulate our nodes/blocks. For the greedy move phase, we move these blocks into the communities of their neighboring blocks. So we need to know who are the neighboring communities of a block. This is saved in `NeighborComm` which is the vector of index of neighboring community for a block. 

`NeighborComm` is recalculated at each iteration. I'm sure there is a way to update it while the algorithm is running. 



## Edge in partition

The gain in modularity depends on the number of edges that we lose or gain when a set of nodes leaves or joins our community. 
For this, it was necessary to have knowledge of the edges which are at the border of our set (the edges which are likely to belong to our community if we gain nodes).  

It would be nice to be able to keep and update this information (the edge border of our communities) during execution.



## Big partition access

During execution, we manipulate several partitions. The partition of our blocks and the partition of the largest communities that we are building. Our blocks are contained in these largest communities. But to obtain the index of the large community in which a block is, we must go back to the node level.

I added a method in partition.cpp to obtain a node of a block and thus be able to obtain index of the large community to which this block belongs. 

```cpp
index Partition::giveOne(index s) const {
    assert(s <= omega);
    index i= std::numeric_limits<uint64_t>::max();
    for (index e = 0; e < this->z; ++e) {
        if (data[e] == s) {
            i = e;
            break;
        }
    }
    return i;
}
```

This method returns either the index of a node belonging to the block, or `std::numeric_limits<uint64_t>::max()`. 

There is clearly a way to be more efficient. For example, trying to make the large partition a block partition.

# License

[MIT](https://choosealicense.com/licenses/mit/)