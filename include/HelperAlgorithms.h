#ifndef HELPER_ALGORITHMS_H
#define HELPER_ALGORITHMS_H

#include "Graph.h"
#include <vector>
#include <chrono>

// A clique is represented as a vector of vertex IDs
typedef std::vector<VertexID> Clique;

// Static class for helper algorithms
class HelperAlgorithms
{
public:
    // MC-DD: Degeneracy-based maximum clique algorithm
    static Clique runMCDD(const Graph &graph);

    // MC-EGO: Ego-centric maximum clique algorithm
    static Clique runMCEGO(const Graph &graph);

    // Improved MC-EGO: Enhanced ego-centric maximum clique algorithm
    static Clique runImprovedMCEGO(const Graph &graph);

private:
    // Helper methods for MC-DD
    static Clique findMaximalCliqueByDegeneracy(const Graph &graph);
    static int colorByDegeneracy(const Graph &graph, const std::vector<VertexID> &vertices);

    // Helper methods for MC-EGO
    static Clique exploreEgoNetwork(const Graph &graph, VertexID center);
    static Clique improvedExploreEgoNetwork(const Graph &graph, VertexID center);
    static Clique findMaximalCliqueInDenseNeighborhood(const Graph &graph, VertexID center,
                                                       const std::vector<VertexID> &neighbors);
    static int computeUpperBound(const Graph &graph, const std::vector<VertexID> &candidates);
};

#endif // HELPER_ALGORITHMS_H