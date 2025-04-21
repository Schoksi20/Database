#ifndef HELPER_ALGORITHMS_H
#define HELPER_ALGORITHMS_H

#include "Graph.h"
#include <vector>
#include <chrono>

typedef std::vector<VertexID> Clique;

class HelperAlgorithms
{
public:
    static Clique runMCDD(const Graph &graph);

    static Clique runMCEGO(const Graph &graph);

    static Clique runImprovedMCEGO(const Graph &graph);

private:
    static Clique findMaximalCliqueByDegeneracy(const Graph &graph);
    static int colorByDegeneracy(const Graph &graph, const std::vector<VertexID> &vertices);

    static Clique exploreEgoNetwork(const Graph &graph, VertexID center);
    static Clique improvedExploreEgoNetwork(const Graph &graph, VertexID center);
    static Clique findMaximalCliqueInDenseNeighborhood(const Graph &graph, VertexID center,
                                                       const std::vector<VertexID> &neighbors);
    static int computeUpperBound(const Graph &graph, const std::vector<VertexID> &candidates);
};

#endif