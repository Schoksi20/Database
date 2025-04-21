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

    void initialize();

    void branchAndBound();

    std::vector<VertexID> getBestClique() const;

    void reportStatistics() const;

private:
    void search(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);

    void reduceCandidates(std::vector<VertexID> &candidates, int currentCliqueSize);
    std::vector<VertexID> intersectCandidates(const std::vector<VertexID> &candidates, VertexID v);
    std::vector<VertexID> fastIntersectCandidates(const std::vector<VertexID> &candidates, VertexID v);

    void advancedKernelization(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);
    void advancedKernelizationWithDominance(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);
    void kernelizeByDegreeReduction(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates);

    void constructBitsetForCandidates(const std::vector<VertexID> &candidates);
    int countCommonNeighbors(VertexID v, VertexID w);
    bool isClique(const std::vector<VertexID> &vertices);
    int getBitsetIntersectionSize(VertexID v, const std::vector<VertexID> &vertices);

    bool shouldUseBitMatrix(const std::vector<VertexID> &candidates);
    void constructBitMatrix(const std::vector<VertexID> &candidates);
    bool areBitMatrixAdjacent(VertexID u, VertexID v);

    struct FoldingOperation
    {
        VertexID center;
        VertexID neighbor1;
        VertexID neighbor2;
    };
    std::vector<FoldingOperation> foldingOperations;
    std::set<std::pair<VertexID, VertexID>> virtualEdges;
    std::unordered_set<VertexID> foldedVertices;

    void enhancedFolding(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates,
                         VertexID v, VertexID n1, VertexID n2, std::vector<int> &localDegree);

    bool shouldSwitchToMatrix(const std::vector<VertexID> &candidates);
    bool shouldUseMatrix(const std::vector<VertexID> &candidates) const;

    std::vector<VertexID> determineBranchOrder(const std::vector<VertexID> &candidates);
    std::vector<VertexID> computeDegeneracyOrdering(const std::vector<VertexID> &candidates);

    int greedyColoring(const std::vector<VertexID> &candidates);
    int advancedColoring(const std::vector<VertexID> &candidates);
    int dsaturColoring(const std::vector<VertexID> &candidates);

    void searchTriangle(const std::vector<VertexID> &currentClique, const std::vector<VertexID> &candidates);

    void constructMatrix(const std::vector<VertexID> &candidates);
    bool matrixAdjacent(VertexID u, VertexID v) const;
    void setMatrixBit(int row, int col);
    bool getMatrixBit(int row, int col) const;

    const Graph &graphRef;

    std::vector<VertexID> bestClique;

    long long branches;
    int maxDepth;

    bool usingMatrix;
    std::vector<VertexID> matrixMapping;
    std::vector<VertexID> inverseMapping;
    std::vector<unsigned char> adjacencyMatrix;
    int matrixSize;

    bool usingBitsets;
    std::vector<std::vector<uint32_t>> vertexBitsets;
    std::unordered_map<VertexID, int> vertexToConsecutiveIdx;
    std::vector<VertexID> consecutiveIdxToVertex;

    bool usingBitMatrix;
    std::unordered_map<VertexID, int> vertexToMatrixIndex;
    std::vector<VertexID> matrixToVertex;
    std::vector<uint32_t> bitMatrix;
};

#endif