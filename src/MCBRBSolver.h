#ifndef MCBRB_SOLVER_H
#define MCBRB_SOLVER_H

#include "Graph.h"
#include "HelperAlgorithms.h"
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <set>

class MCBRBSolver
{
public:
    MCBRBSolver(const Graph &graph);
    ~MCBRBSolver() = default;

    // Initialize using helper algorithms (MC-DD and MC-EGO)
    void initialize();

    // Main branch-and-bound search routine
    void branchAndBound();

    // Return the best clique found
    std::vector<VertexID> getBestClique() const;

    // Status reporting
    void reportStatistics() const;

private:
    // Core algorithm components
    void search(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);

    // Candidate set management
    void reduceCandidates(std::vector<VertexID> &candidates, int currentCliqueSize);
    std::vector<VertexID> intersectCandidates(const std::vector<VertexID> &candidates, VertexID v);
    std::vector<VertexID> fastIntersectCandidates(const std::vector<VertexID> &candidates, VertexID v);

    // Advanced kernelization techniques
    void advancedKernelization(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);
    void advancedKernelizationWithDominance(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);
    void kernelizeByDegreeReduction(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);

    // Bitset operations for performance
    void constructBitsetForCandidates(const std::vector<VertexID> &candidates);
    int countCommonNeighbors(VertexID v, VertexID w);
    bool isClique(const std::vector<VertexID> &vertices);
    int getBitsetIntersectionSize(VertexID v, const std::vector<VertexID> &vertices);

    // Bit matrix operations
    bool shouldUseBitMatrix(const std::vector<VertexID> &candidates);
    void constructBitMatrix(const std::vector<VertexID> &candidates);
    bool areBitMatrixAdjacent(VertexID u, VertexID v);

    // Folding operation tracking
    struct FoldingOperation
    {
        VertexID center;    // The degree-2 vertex
        VertexID neighbor1; // First neighbor
        VertexID neighbor2; // Second neighbor
    };
    std::vector<FoldingOperation> foldingOperations;
    std::set<std::pair<VertexID, VertexID>> virtualEdges;
    std::unordered_set<VertexID> foldedVertices;

    // Enhanced folding function
    void enhancedFolding(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates,
                         VertexID v, VertexID n1, VertexID n2, std::vector<int> &localDegree);

    // Dynamic representation selection
    bool shouldSwitchToMatrix(const std::vector<VertexID> &candidates);
    bool shouldUseMatrix(const std::vector<VertexID> &candidates) const;

    // Branch ordering strategies
    std::vector<VertexID> determineBranchOrder(const std::vector<VertexID> &candidates);
    std::vector<VertexID> computeDegeneracyOrdering(const std::vector<VertexID> &candidates);

    // Upper bound computation methods
    int greedyColoring(const std::vector<VertexID> &candidates);
    int advancedColoring(const std::vector<VertexID> &candidates);
    int dsaturColoring(const std::vector<VertexID> &candidates);

    // Special cases
    void searchTriangle(const std::vector<VertexID> &currentClique, const std::vector<VertexID> &candidates);

    // Matrix handling for dense subgraphs
    void constructMatrix(const std::vector<VertexID> &candidates);
    bool matrixAdjacent(VertexID u, VertexID v) const;
    void setMatrixBit(int row, int col);
    bool getMatrixBit(int row, int col) const;

    // Reference to the input graph
    const Graph &graphRef;

    // State of the search
    std::vector<VertexID> bestClique;

    // Statistics tracking
    long long branches;
    int maxDepth;

    // Matrix representation
    bool usingMatrix;
    std::vector<VertexID> matrixMapping;        // Map original vertices to matrix indices
    std::vector<VertexID> inverseMapping;       // Map matrix indices to original vertices
    std::vector<unsigned char> adjacencyMatrix; // Bit-packed adjacency matrix
    int matrixSize;

    // Bitset representation
    bool usingBitsets;
    std::vector<std::vector<uint32_t>> vertexBitsets;
    std::unordered_map<VertexID, int> vertexToConsecutiveIdx;
    std::vector<VertexID> consecutiveIdxToVertex;

    // Bit matrix representation variables
    bool usingBitMatrix;
    std::unordered_map<VertexID, int> vertexToMatrixIndex;
    std::vector<VertexID> matrixToVertex;
    std::vector<uint32_t> bitMatrix;
};

#endif // MCBRB_SOLVER_H