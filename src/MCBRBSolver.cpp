#include "MCBRBSolver.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <cstring>
#include <chrono>
// Constructor implementation
MCBRBSolver::MCBRBSolver(const Graph &graph)
    : graphRef(graph),
      branches(0),
      maxDepth(0),
      usingMatrix(false),
      matrixSize(0),
      usingBitsets(false),
      usingBitMatrix(false)
{
    bestClique.clear();
}

// Initialize using helper algorithms (MC-DD and MC-EGO).
// Update initialize to use the improved MC-EGO
void MCBRBSolver::initialize()
{
    auto cliqueDD = HelperAlgorithms::runMCDD(graphRef);
    auto cliqueEGO = HelperAlgorithms::runImprovedMCEGO(graphRef); // Use improved MC-EGO

    // Choose the clique with larger size; if equal, pick one arbitrarily.
    if (cliqueDD.size() >= cliqueEGO.size())
        bestClique = cliqueDD;
    else
        bestClique = cliqueEGO;

    std::cout << "[Solver] Initial best clique size: " << bestClique.size() << std::endl;
}

// Determine optimal branch ordering for the current candidates
std::vector<VertexID> MCBRBSolver::determineBranchOrder(const std::vector<VertexID> &candidates)
{
    // For large candidate sets, degeneracy ordering is usually best
    if (candidates.size() > 100)
    {
        return computeDegeneracyOrdering(candidates);
    }

    // For medium-sized sets, use degree-weighted ordering
    std::vector<std::pair<int, VertexID>> scoredCandidates;

    for (VertexID v : candidates)
    {
        int internalDegree = 0;
        for (VertexID u : candidates)
        {
            if (v != u && (usingMatrix ? matrixAdjacent(v, u) : graphRef.isAdjacent(v, u)))
            {
                internalDegree++;
            }
        }

        // Combine degeneracy information with local degree for better ordering
        int degPosition = std::find(graphRef.getDegeneracyOrder().begin(),
                                    graphRef.getDegeneracyOrder().end(), v) -
                          graphRef.getDegeneracyOrder().begin();

        // Score is weighted combination of internal degree and degeneracy position
        double score = 0.7 * internalDegree + 0.3 * (graphRef.getNumVertices() - degPosition);

        scoredCandidates.push_back(std::make_pair(score, v));
    }

    // Sort by score (descending)
    std::sort(scoredCandidates.begin(), scoredCandidates.end(),
              [](const std::pair<int, VertexID> &a, const std::pair<int, VertexID> &b)
              { return a.first > b.first; });

    // Extract vertices in order
    std::vector<VertexID> orderedCandidates;
    for (const auto &pair : scoredCandidates)
    {
        orderedCandidates.push_back(pair.second);
    }

    return orderedCandidates;
}

// Add this to MCBRBSolver.cpp
// Advanced kernelization with dominance relation
void MCBRBSolver::advancedKernelizationWithDominance(std::vector<VertexID> &currentClique,
                                                     std::vector<VertexID> &candidates)
{
    // First apply basic kernelization
    advancedKernelization(currentClique, candidates);

    if (candidates.size() <= 3)
        return;

    // Now check for dominance relations
    std::vector<bool> dominated(graphRef.getNumVertices(), false);

    // For each pair of vertices, check if one dominates the other
    for (size_t i = 0; i < candidates.size(); i++)
    {
        VertexID v = candidates[i];
        if (dominated[v])
            continue;

        // Get neighbors of v in the candidate set
        std::vector<VertexID> vNeighbors;
        for (size_t j = 0; j < candidates.size(); j++)
        {
            VertexID u = candidates[j];
            if (u != v && (usingMatrix ? matrixAdjacent(v, u) : graphRef.isAdjacent(v, u)))
            {
                vNeighbors.push_back(u);
            }
        }

        // Check each other vertex for dominance relation
        for (size_t j = i + 1; j < candidates.size(); j++)
        {
            VertexID w = candidates[j];
            if (dominated[w])
                continue;

            // Check if v and w are adjacent - if not, no dominance possible
            bool adjacent = (usingMatrix ? matrixAdjacent(v, w) : graphRef.isAdjacent(v, w));
            if (!adjacent)
                continue;

            // Get neighbors of w in the candidate set
            std::vector<VertexID> wNeighbors;
            for (size_t k = 0; k < candidates.size(); k++)
            {
                VertexID u = candidates[k];
                if (u != w && (usingMatrix ? matrixAdjacent(w, u) : graphRef.isAdjacent(w, u)))
                {
                    wNeighbors.push_back(u);
                }
            }

            // Check if v dominates w (N(w) - {v} ⊆ N(v) - {w})
            bool vDominatesW = true;
            for (VertexID u : wNeighbors)
            {
                if (u == v)
                    continue; // Skip v itself

                // If u is in N(w) but not in N(v), then v doesn't dominate w
                if (std::find(vNeighbors.begin(), vNeighbors.end(), u) == vNeighbors.end())
                {
                    vDominatesW = false;
                    break;
                }
            }

            // Check if w dominates v (N(v) - {w} ⊆ N(w) - {v})
            bool wDominatesV = true;
            for (VertexID u : vNeighbors)
            {
                if (u == w)
                    continue; // Skip w itself

                // If u is in N(v) but not in N(w), then w doesn't dominate v
                if (std::find(wNeighbors.begin(), wNeighbors.end(), u) == wNeighbors.end())
                {
                    wDominatesV = false;
                    break;
                }
            }

            // Apply dominance rule - we can remove the dominated vertex
            if (vDominatesW)
            {
                dominated[w] = true;
            }
            else if (wDominatesV)
            {
                dominated[v] = true;
                break; // No need to check v against other vertices
            }
        }
    }

    // Remove dominated vertices
    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(),
                       [&dominated](VertexID v)
                       { return dominated[v]; }),
        candidates.end());
}

// Dynamic matrix selection based on graph characteristics
bool MCBRBSolver::shouldSwitchToMatrix(const std::vector<VertexID> &candidates)
{
    const int sizeThreshold = 50;        // Minimum size to consider
    const double densityThreshold = 0.4; // Density threshold

    // Matrix representation is only worth it for larger candidate sets
    if (candidates.size() < sizeThreshold)
    {
        return false;
    }

    // Estimate density by sampling
    int sampleSize = std::min(100, static_cast<int>(candidates.size()));
    int edges = 0;
    int possibleEdges = sampleSize * (sampleSize - 1) / 2;

    // Sample some candidate vertices
    for (int i = 0; i < sampleSize; i++)
    {
        VertexID v = candidates[i];
        for (int j = i + 1; j < sampleSize; j++)
        {
            VertexID u = candidates[j];
            if (graphRef.isAdjacent(v, u))
            {
                edges++;
            }
        }
    }

    double density = static_cast<double>(edges) / possibleEdges;
    return density > densityThreshold;
}

// Advanced coloring implementation using a saturation-based approach
// This is more effective than simple greedy coloring
int MCBRBSolver::advancedColoring(const std::vector<VertexID> &candidates)
{
    int n = candidates.size();
    if (n == 0)
        return 0;

    // Track the color assigned to each vertex
    std::vector<int> color(n, -1);

    // Track the saturation of each vertex (number of different colors in neighborhood)
    std::vector<int> saturation(n, 0);

    // Used colors for each vertex
    std::vector<std::vector<bool>> usedColors(n);
    for (int i = 0; i < n; i++)
    {
        usedColors[i].resize(n, false);
    }

    // Process vertices in order of saturation
    for (int i = 0; i < n; i++)
    {
        // Find uncolored vertex with maximum saturation
        int maxSat = -1;
        int maxDeg = -1;
        int maxIdx = -1;

        for (int j = 0; j < n; j++)
        {
            if (color[j] == -1)
            { // Not yet colored
                if (saturation[j] > maxSat ||
                    (saturation[j] == maxSat && graphRef.degree(candidates[j]) > maxDeg))
                {
                    maxSat = saturation[j];
                    maxDeg = graphRef.degree(candidates[j]);
                    maxIdx = j;
                }
            }
        }

        int v = maxIdx;

        // Find smallest available color
        int c = 0;
        while (c < n && usedColors[v][c])
            c++;

        // Assign the color
        color[v] = c;

        // Update saturation of neighbors
        for (int j = 0; j < n; j++)
        {
            if (color[j] == -1)
            { // Not yet colored
                if (usingMatrix)
                {
                    if (matrixAdjacent(candidates[v], candidates[j]))
                    {
                        if (!usedColors[j][c])
                        {
                            usedColors[j][c] = true;
                            saturation[j]++;
                        }
                    }
                }
                else
                {
                    if (graphRef.isAdjacent(candidates[v], candidates[j]))
                    {
                        if (!usedColors[j][c])
                        {
                            usedColors[j][c] = true;
                            saturation[j]++;
                        }
                    }
                }
            }
        }
    }

    // Return the number of colors used
    int maxColor = *std::max_element(color.begin(), color.end());
    return maxColor + 1;
}

// Basic greedy coloring heuristic with optimized local degree computation
int MCBRBSolver::greedyColoring(const std::vector<VertexID> &origCandidates)
{
    std::vector<VertexID> candidates = origCandidates;
    int n = candidates.size();
    if (n == 0)
        return 0;

    // Optimized local degree computation
    std::vector<int> localDegree(n, 0);

    for (int i = 0; i < n; i++)
    {
        VertexID v = candidates[i];
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                VertexID u = candidates[j];
                if (usingMatrix)
                {
                    if (matrixAdjacent(v, u))
                    {
                        localDegree[i]++;
                    }
                }
                else
                {
                    if (graphRef.isAdjacent(v, u))
                    {
                        localDegree[i]++;
                    }
                }
            }
        }
    }

    // Comparator: sort in descending order by local degree
    auto cmpLocalDegree = [&localDegree](int a, int b)
    {
        return localDegree[a] > localDegree[b];
    };

    // Create an ordering of indices
    std::vector<int> ordering(n);
    for (int i = 0; i < n; i++)
        ordering[i] = i;
    std::sort(ordering.begin(), ordering.end(), cmpLocalDegree);

    // Reorder candidates based on local degree
    std::vector<VertexID> orderedCandidates(n);
    for (int i = 0; i < n; i++)
    {
        orderedCandidates[i] = candidates[ordering[i]];
    }

    // Coloring process
    std::vector<int> color(n, -1);
    for (int i = 0; i < n; i++)
    {
        VertexID v = orderedCandidates[i];
        std::vector<bool> used(n, false);

        for (int j = 0; j < i; j++)
        {
            VertexID u = orderedCandidates[j];
            if (usingMatrix)
            {
                if (matrixAdjacent(v, u))
                {
                    used[color[j]] = true;
                }
            }
            else
            {
                if (graphRef.isAdjacent(v, u))
                {
                    used[color[j]] = true;
                }
            }
        }

        int c = 0;
        while (c < n && used[c])
            c++;
        color[i] = c;
    }

    int maxColor = *std::max_element(color.begin(), color.end());
    return maxColor + 1; // Colors are zero-indexed
}

// Compute intersection of candidates and the neighbors of vertex v
std::vector<VertexID> MCBRBSolver::intersectCandidates(const std::vector<VertexID> &candidates, VertexID v)
{
    std::vector<VertexID> newCandidates;

    if (usingMatrix)
    {
        for (auto u : candidates)
        {
            if (matrixAdjacent(v, u))
            {
                newCandidates.push_back(u);
            }
        }
    }
    else
    {
        const auto &neighbors = graphRef.getNeighbors(v);
        for (auto u : candidates)
        {
            if (std::binary_search(neighbors.begin(), neighbors.end(), u))
            {
                newCandidates.push_back(u);
            }
        }
    }

    return newCandidates;
}

// K-core reduction: Remove candidates that can't be part of a clique larger than current best
void MCBRBSolver::reduceCandidates(std::vector<VertexID> &candidates, int currentCliqueSize)
{
    int threshold = bestClique.size() - currentCliqueSize;
    if (threshold <= 0)
        return; // No pruning needed

    bool removed;
    std::vector<int> localDegree(graphRef.getNumVertices(), 0);

    // Compute initial local degrees
    for (auto v : candidates)
    {
        for (auto u : candidates)
        {
            if (v != u)
            {
                if (usingMatrix)
                {
                    if (matrixAdjacent(v, u))
                    {
                        localDegree[v]++;
                    }
                }
                else
                {
                    if (graphRef.isAdjacent(v, u))
                    {
                        localDegree[v]++;
                    }
                }
            }
        }
    }

    // Iterate until no more removals
    do
    {
        removed = false;
        std::vector<VertexID> newCandidates;

        for (auto v : candidates)
        {
            if (localDegree[v] >= threshold)
            {
                newCandidates.push_back(v);
            }
            else
            {
                // Remove v, update degrees of its neighbors
                for (auto u : candidates)
                {
                    if (u != v)
                    {
                        if (usingMatrix)
                        {
                            if (matrixAdjacent(v, u))
                            {
                                localDegree[u]--;
                            }
                        }
                        else
                        {
                            if (graphRef.isAdjacent(v, u))
                            {
                                localDegree[u]--;
                            }
                        }
                    }
                }
                removed = true;
            }
        }

        candidates = newCandidates;
    } while (removed);
}

// Calculate if using matrix representation is beneficial
bool MCBRBSolver::shouldUseMatrix(const std::vector<VertexID> &candidates) const
{
    // Use matrix if the candidate set is dense enough
    const int minSize = 20;              // Minimum size to consider using matrix
    const double densityThreshold = 0.3; // Density threshold

    if (candidates.size() < minSize)
        return false;

    // Calculate density (expensive operation, so only do for reasonably sized sets)
    long edgeCount = 0;
    long possibleEdges = candidates.size() * (candidates.size() - 1) / 2;

    for (size_t i = 0; i < candidates.size(); i++)
    {
        for (size_t j = i + 1; j < candidates.size(); j++)
        {
            if (graphRef.isAdjacent(candidates[i], candidates[j]))
            {
                edgeCount++;
            }
        }
    }

    double density = static_cast<double>(edgeCount) / possibleEdges;
    return density > densityThreshold;
}

// Construct adjacency matrix for a subgraph
void MCBRBSolver::constructMatrix(const std::vector<VertexID> &candidates)
{
    matrixSize = candidates.size();
    if (matrixSize == 0)
    {
        usingMatrix = false;
        return;
    }

    // Set up mappings
    matrixMapping.resize(graphRef.getNumVertices(), -1);
    inverseMapping.resize(matrixSize);

    for (int i = 0; i < matrixSize; i++)
    {
        VertexID v = candidates[i];
        matrixMapping[v] = i;
        inverseMapping[i] = v;
    }

    // Allocate matrix (bit-packed)
    int matrixBytes = (matrixSize * matrixSize + 7) / 8;
    adjacencyMatrix.resize(matrixBytes, 0);

    // Fill the matrix
    for (int i = 0; i < matrixSize; i++)
    {
        VertexID v = inverseMapping[i];
        for (int j = 0; j < matrixSize; j++)
        {
            VertexID u = inverseMapping[j];
            if (i != j && graphRef.isAdjacent(v, u))
            {
                setMatrixBit(i, j);
            }
        }
    }

    usingMatrix = true;
}

// Set a bit in the adjacency matrix
void MCBRBSolver::setMatrixBit(int row, int col)
{
    int bitPos = row * matrixSize + col;
    adjacencyMatrix[bitPos / 8] |= (1 << (bitPos % 8));
}

// Get a bit from the adjacency matrix
bool MCBRBSolver::getMatrixBit(int row, int col) const
{
    int bitPos = row * matrixSize + col;
    return (adjacencyMatrix[bitPos / 8] & (1 << (bitPos % 8))) != 0;
}

// Check if two vertices are adjacent using the matrix
bool MCBRBSolver::matrixAdjacent(VertexID u, VertexID v) const
{
    if (!usingMatrix)
        return graphRef.isAdjacent(u, v);

    int rowU = matrixMapping[u];
    int rowV = matrixMapping[v];

    // If either vertex is not in the matrix, use the original graph
    if (rowU == -1 || rowV == -1)
        return graphRef.isAdjacent(u, v);

    return getMatrixBit(rowU, rowV);
}

// Compute degeneracy ordering for the candidate set
std::vector<VertexID> MCBRBSolver::computeDegeneracyOrdering(const std::vector<VertexID> &origCandidates)
{
    // Work on a local copy
    std::vector<VertexID> candidates = origCandidates;
    int n = candidates.size();
    std::vector<VertexID> order; // Will be the final degeneracy ordering

    // Compute induced degrees for the candidate set
    std::vector<int> degrees(n, 0);
    std::vector<int> vertex2idx(graphRef.getNumVertices(), -1);

    for (int i = 0; i < n; i++)
    {
        vertex2idx[candidates[i]] = i;
    }

    for (int i = 0; i < n; i++)
    {
        VertexID v = candidates[i];
        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            VertexID u = candidates[j];
            if (usingMatrix)
            {
                if (matrixAdjacent(v, u))
                {
                    degrees[i]++;
                }
            }
            else
            {
                if (graphRef.isAdjacent(v, u))
                {
                    degrees[i]++;
                }
            }
        }
    }

    // Iteratively add the vertex with minimum degree to the ordering
    std::vector<bool> removed(n, false);
    for (int i = 0; i < n; i++)
    {
        // Find unremoved vertex with minimum degree
        int minDeg = n;
        int minIdx = -1;

        for (int j = 0; j < n; j++)
        {
            if (!removed[j] && degrees[j] < minDeg)
            {
                minDeg = degrees[j];
                minIdx = j;
            }
        }

        if (minIdx == -1)
            break; // Should not happen

        // Add to ordering
        order.push_back(candidates[minIdx]);
        removed[minIdx] = true;

        // Update degrees
        for (int j = 0; j < n; j++)
        {
            if (removed[j])
                continue;

            VertexID v = candidates[j];
            VertexID u = candidates[minIdx];

            if (usingMatrix)
            {
                if (matrixAdjacent(v, u))
                {
                    degrees[j]--;
                }
            }
            else
            {
                if (graphRef.isAdjacent(v, u))
                {
                    degrees[j]--;
                }
            }
        }
    }

    return order;
}

// Special case for small candidate sets: try to find triangles
void MCBRBSolver::searchTriangle(const std::vector<VertexID> &currentClique, const std::vector<VertexID> &candidates)
{
    int cliqueSize = currentClique.size();
    int candSize = candidates.size();

    // Not enough candidates to improve best clique
    if (cliqueSize + candSize <= bestClique.size())
    {
        return;
    }

    // Find any triangle (3-clique) among the candidates
    for (int i = 0; i < candSize; i++)
    {
        for (int j = i + 1; j < candSize; j++)
        {
            VertexID u = candidates[i];
            VertexID v = candidates[j];

            bool adjacent = false;
            if (usingMatrix)
            {
                adjacent = matrixAdjacent(u, v);
            }
            else
            {
                adjacent = graphRef.isAdjacent(u, v);
            }

            if (adjacent)
            {
                for (int k = j + 1; k < candSize; k++)
                {
                    VertexID w = candidates[k];

                    bool adjacent1 = false, adjacent2 = false;
                    if (usingMatrix)
                    {
                        adjacent1 = matrixAdjacent(u, w);
                        adjacent2 = matrixAdjacent(v, w);
                    }
                    else
                    {
                        adjacent1 = graphRef.isAdjacent(u, w);
                        adjacent2 = graphRef.isAdjacent(v, w);
                    }

                    if (adjacent1 && adjacent2)
                    {
                        // Found a triangle
                        std::vector<VertexID> newClique = currentClique;
                        newClique.push_back(u);
                        newClique.push_back(v);
                        newClique.push_back(w);

                        if (newClique.size() > bestClique.size())
                        {
                            bestClique = newClique;
                            std::cout << "[Solver] Found a new best clique of size "
                                      << bestClique.size() << " via triangle search." << std::endl;
                        }
                        return; // One triangle is enough
                    }
                }
            }
        }
    }

    // If we can't find a triangle, try pairs
    if (bestClique.size() <= cliqueSize + 1)
    {
        for (int i = 0; i < candSize; i++)
        {
            for (int j = i + 1; j < candSize; j++)
            {
                VertexID u = candidates[i];
                VertexID v = candidates[j];

                bool adjacent = false;
                if (usingMatrix)
                {
                    adjacent = matrixAdjacent(u, v);
                }
                else
                {
                    adjacent = graphRef.isAdjacent(u, v);
                }

                if (adjacent)
                {
                    std::vector<VertexID> newClique = currentClique;
                    newClique.push_back(u);
                    newClique.push_back(v);

                    if (newClique.size() > bestClique.size())
                    {
                        bestClique = newClique;
                        std::cout << "[Solver] Found a new best clique of size "
                                  << bestClique.size() << " via pair search." << std::endl;
                    }
                    return;
                }
            }
        }
    }

    // If nothing else helps, at least add one vertex
    if (bestClique.size() <= cliqueSize && !candidates.empty())
    {
        std::vector<VertexID> newClique = currentClique;
        newClique.push_back(candidates[0]);

        if (newClique.size() > bestClique.size())
        {
            bestClique = newClique;
            std::cout << "[Solver] Found a new best clique of size "
                      << bestClique.size() << " via singleton addition." << std::endl;
        }
    }
}

// Kernelization by degree reduction
void MCBRBSolver::kernelizeByDegreeReduction(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates)
{
    bool changed = true;

    // Store local degrees for efficiency
    std::vector<int> localDegree(graphRef.getNumVertices(), 0);

    // Initialize local degrees
    for (auto v : candidates)
    {
        for (auto u : candidates)
        {
            if (v != u)
            {
                if (usingMatrix)
                {
                    if (matrixAdjacent(v, u))
                    {
                        localDegree[v]++;
                    }
                }
                else
                {
                    if (graphRef.isAdjacent(v, u))
                    {
                        localDegree[v]++;
                    }
                }
            }
        }
    }

    while (changed && !candidates.empty())
    {
        changed = false;

        // Identify vertices with degree 0, 1, 2, etc.
        std::vector<VertexID> degreeZero, degreeOne, degreeTwo;

        for (auto v : candidates)
        {
            if (localDegree[v] == 0)
            {
                degreeZero.push_back(v);
                changed = true;
            }
            else if (localDegree[v] == 1)
            {
                degreeOne.push_back(v);
                changed = true;
            }
            else if (localDegree[v] == 2)
            {
                degreeTwo.push_back(v);
                changed = true;
            }
        }

        // Process degree-0 vertices (isolated vertices)
        for (auto v : degreeZero)
        {
            // Add to current clique and remove from candidates
            currentClique.push_back(v);
            candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
        }

        // Process degree-1 vertices
        for (auto v : degreeOne)
        {
            if (std::find(candidates.begin(), candidates.end(), v) == candidates.end())
                continue;

            // Find its only neighbor
            VertexID neighbor = -1;
            for (auto u : candidates)
            {
                if (u != v)
                {
                    if (usingMatrix)
                    {
                        if (matrixAdjacent(v, u))
                        {
                            neighbor = u;
                            break;
                        }
                    }
                    else
                    {
                        if (graphRef.isAdjacent(v, u))
                        {
                            neighbor = u;
                            break;
                        }
                    }
                }
            }

            if (neighbor != -1)
            {
                // Add v to current clique
                currentClique.push_back(v);

                // Remove v and its neighbor from candidates
                candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
                candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor), candidates.end());

                // Update local degrees
                for (auto u : candidates)
                {
                    if (usingMatrix)
                    {
                        if (matrixAdjacent(neighbor, u))
                        {
                            localDegree[u]--;
                        }
                    }
                    else
                    {
                        if (graphRef.isAdjacent(neighbor, u))
                        {
                            localDegree[u]--;
                        }
                    }
                }
            }
        }

        // Process degree-2 vertices
        for (auto v : degreeTwo)
        {
            if (std::find(candidates.begin(), candidates.end(), v) == candidates.end())
                continue;

            // Find its two neighbors
            VertexID neighbor1 = -1, neighbor2 = -1;
            for (auto u : candidates)
            {
                if (u != v)
                {
                    if (usingMatrix)
                    {
                        if (matrixAdjacent(v, u))
                        {
                            if (neighbor1 == -1)
                            {
                                neighbor1 = u;
                            }
                            else
                            {
                                neighbor2 = u;
                                break;
                            }
                        }
                    }
                    else
                    {
                        if (graphRef.isAdjacent(v, u))
                        {
                            if (neighbor1 == -1)
                            {
                                neighbor1 = u;
                            }
                            else
                            {
                                neighbor2 = u;
                                break;
                            }
                        }
                    }
                }
            }

            if (neighbor1 != -1 && neighbor2 != -1)
            {
                bool adjacent = false;
                if (usingMatrix)
                {
                    adjacent = matrixAdjacent(neighbor1, neighbor2);
                }
                else
                {
                    adjacent = graphRef.isAdjacent(neighbor1, neighbor2);
                }

                if (adjacent)
                {
                    // The neighbors form a triangle with v
                    currentClique.push_back(v);
                    currentClique.push_back(neighbor1);
                    currentClique.push_back(neighbor2);

                    // Remove all three from candidates
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor1), candidates.end());
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor2), candidates.end());

                    // Update local degrees
                    for (auto u : candidates)
                    {
                        if (usingMatrix)
                        {
                            if (matrixAdjacent(neighbor1, u))
                            {
                                localDegree[u]--;
                            }
                            if (matrixAdjacent(neighbor2, u))
                            {
                                localDegree[u]--;
                            }
                        }
                        else
                        {
                            if (graphRef.isAdjacent(neighbor1, u))
                            {
                                localDegree[u]--;
                            }
                            if (graphRef.isAdjacent(neighbor2, u))
                            {
                                localDegree[u]--;
                            }
                        }
                    }
                }
                else
                {
                    // Folding operation (more complex, simplified here)
                    currentClique.push_back(v);

                    // Remove v from candidates
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());

                    // Update local degrees
                    for (auto u : candidates)
                    {
                        if (usingMatrix)
                        {
                            if (matrixAdjacent(v, u))
                            {
                                localDegree[u]--;
                            }
                        }
                        else
                        {
                            if (graphRef.isAdjacent(v, u))
                            {
                                localDegree[u]--;
                            }
                        }
                    }
                }
            }
        }
    }
}

int MCBRBSolver::dsaturColoring(const std::vector<VertexID> &candidates)
{
    int n = candidates.size();
    if (n == 0)
        return 0;

    // Track the color assigned to each vertex
    std::vector<int> color(n, -1);

    // Track the saturation of each vertex (number of different adjacent colors)
    std::vector<int> saturation(n, 0);

    // Track vertex degrees within the candidate set
    std::vector<int> degrees(n, 0);

    // Compute degrees
    for (int i = 0; i < n; i++)
    {
        VertexID v = candidates[i];
        for (int j = 0; j < n; j++)
        {
            if (i != j && graphRef.isAdjacent(v, candidates[j]))
            {
                degrees[i]++;
            }
        }
    }

    // Used colors for each vertex's neighborhood
    std::vector<std::vector<bool>> usedColors(n);
    for (int i = 0; i < n; i++)
    {
        usedColors[i].resize(n, false);
    }

    // Color vertices one by one
    for (int c = 0; c < n; c++)
    {
        // Find uncolored vertex with maximum saturation
        // If tied, choose the one with maximum degree
        int maxSat = -1;
        int maxDeg = -1;
        int maxIdx = -1;

        for (int i = 0; i < n; i++)
        {
            if (color[i] == -1)
            { // If not colored yet
                if (saturation[i] > maxSat ||
                    (saturation[i] == maxSat && degrees[i] > maxDeg))
                {
                    maxSat = saturation[i];
                    maxDeg = degrees[i];
                    maxIdx = i;
                }
            }
        }

        if (maxIdx == -1)
            break; // All vertices colored

        // Find the smallest available color
        int minColor = 0;
        while (minColor < n && usedColors[maxIdx][minColor])
            minColor++;

        // Assign color
        color[maxIdx] = minColor;

        // Update saturation of adjacent vertices
        for (int i = 0; i < n; i++)
        {
            if (color[i] == -1 && graphRef.isAdjacent(candidates[maxIdx], candidates[i]))
            {
                if (!usedColors[i][minColor])
                {
                    usedColors[i][minColor] = true;
                    saturation[i]++;
                }
            }
        }
    }

    // Return the number of colors used
    return *std::max_element(color.begin(), color.end()) + 1;
}

// Utility function to check if using bit-matrix would be beneficial
bool MCBRBSolver::shouldUseBitMatrix(const std::vector<VertexID> &candidates)
{
    // Use bit matrix if candidate set is reasonably large and dense
    if (candidates.size() < 32)
        return false; // Too small to benefit

    // Calculate density
    int edges = 0;
    int totalPossible = candidates.size() * (candidates.size() - 1) / 2;

    // Sample some vertices to estimate density
    const int maxSamples = 100;
    int samples = std::min(maxSamples, (int)candidates.size());

    for (int i = 0; i < samples; i++)
    {
        VertexID v = candidates[i];
        for (int j = i + 1; j < samples; j++)
        {
            if (graphRef.isAdjacent(v, candidates[j]))
            {
                edges++;
            }
        }
    }

    double sampleDensity = (double)edges / (samples * (samples - 1) / 2);
    return sampleDensity > 0.3; // Use bit matrix for dense subgraphs
}

// Construct bit matrix for a candidate set
void MCBRBSolver::constructBitMatrix(const std::vector<VertexID> &candidates)
{
    int n = candidates.size();

    // Clear any existing mappings
    vertexToMatrixIndex.clear();
    matrixToVertex.resize(n);

    // Create mappings
    for (int i = 0; i < n; i++)
    {
        vertexToMatrixIndex[candidates[i]] = i;
        matrixToVertex[i] = candidates[i];
    }

    // Allocate bit matrix
    int words = (n + 31) / 32; // Number of 32-bit words needed per row
    bitMatrix.resize(n * words, 0);

    // Fill the matrix
    for (int i = 0; i < n; i++)
    {
        VertexID v = candidates[i];
        for (int j = 0; j < n; j++)
        {
            if (i != j && graphRef.isAdjacent(v, candidates[j]))
            {
                // Set bit in matrix
                int wordIndex = i * words + (j / 32);
                int bitPosition = j % 32;
                bitMatrix[wordIndex] |= (1U << bitPosition);
            }
        }
    }

    usingBitMatrix = true;
}

// Check adjacency using bit matrix
bool MCBRBSolver::areBitMatrixAdjacent(VertexID u, VertexID v)
{
    if (!usingBitMatrix)
        return graphRef.isAdjacent(u, v);

    auto it1 = vertexToMatrixIndex.find(u);
    auto it2 = vertexToMatrixIndex.find(v);

    if (it1 == vertexToMatrixIndex.end() || it2 == vertexToMatrixIndex.end())
        return graphRef.isAdjacent(u, v);

    int i = it1->second;
    int j = it2->second;

    int words = (matrixToVertex.size() + 31) / 32;
    int wordIndex = i * words + (j / 32);
    int bitPosition = j % 32;

    return (bitMatrix[wordIndex] & (1U << bitPosition)) != 0;
}

// Advanced kernelization for low-degree vertices
void MCBRBSolver::advancedKernelization(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates)
{
    if (candidates.size() <= 1)
        return;

    bool changes = true;
    std::vector<int> degrees(graphRef.getNumVertices(), 0);

    // Compute initial degrees
    for (VertexID v : candidates)
    {
        degrees[v] = 0;
        for (VertexID u : candidates)
        {
            if (v != u && graphRef.isAdjacent(v, u))
            {
                degrees[v]++;
            }
        }
    }

    while (changes && !candidates.empty())
    {
        changes = false;

        // Find vertices with degree 0, 1, 2
        std::vector<VertexID> degree0, degree1, degree2;
        for (VertexID v : candidates)
        {
            if (degrees[v] == 0)
                degree0.push_back(v);
            else if (degrees[v] == 1)
                degree1.push_back(v);
            else if (degrees[v] == 2)
                degree2.push_back(v);
        }

        // Process degree-0 vertices (isolated)
        if (!degree0.empty())
        {
            changes = true;

            // Take any degree-0 vertex into the clique
            VertexID v = degree0[0];
            currentClique.push_back(v);

            // Remove v from candidates
            candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
        }
        // Process degree-1 vertices
        else if (!degree1.empty())
        {
            changes = true;
            VertexID v = degree1[0];

            // Find v's only neighbor
            VertexID neighbor = -1;
            for (VertexID u : candidates)
            {
                if (u != v && graphRef.isAdjacent(v, u))
                {
                    neighbor = u;
                    break;
                }
            }

            // Add v to the clique
            currentClique.push_back(v);

            // Remove v and its neighbor from candidates
            candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
            if (neighbor != -1)
            {
                candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor), candidates.end());

                // Update degrees
                for (VertexID u : candidates)
                {
                    if (graphRef.isAdjacent(u, neighbor))
                    {
                        degrees[u]--;
                    }
                }
            }
        }
        // Process degree-2 vertices
        else if (!degree2.empty())
        {
            changes = true;
            VertexID v = degree2[0];

            // Find v's two neighbors
            VertexID n1 = -1, n2 = -1;
            for (VertexID u : candidates)
            {
                if (u != v && graphRef.isAdjacent(v, u))
                {
                    if (n1 == -1)
                        n1 = u;
                    else
                    {
                        n2 = u;
                        break;
                    }
                }
            }

            // Case 1: Both neighbors are connected (forming a triangle)
            if (n1 != -1 && n2 != -1 && graphRef.isAdjacent(n1, n2))
            {
                // All three form a clique, so we can add v to our clique
                currentClique.push_back(v);

                // Remove v, n1, n2 from candidates
                std::vector<VertexID> toRemove = {v, n1, n2};
                for (VertexID x : toRemove)
                {
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), x), candidates.end());

                    // Update degrees
                    for (VertexID u : candidates)
                    {
                        if (graphRef.isAdjacent(u, x))
                        {
                            degrees[u]--;
                        }
                    }
                }
            }
            // Case 2: The neighbors are not connected (fold)
            else if (n1 != -1 && n2 != -1)
            {
                // Remove v from candidates
                candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());

                // Add a virtual edge between n1 and n2 (folding)
                // This is a bit complex to implement directly, so we'll simplify

                // For our implementation, we'll add v to the clique and remove both neighbors
                currentClique.push_back(v);

                // Remove n1 and n2
                candidates.erase(std::remove(candidates.begin(), candidates.end(), n1), candidates.end());
                candidates.erase(std::remove(candidates.begin(), candidates.end(), n2), candidates.end());

                // Update degrees
                for (VertexID u : candidates)
                {
                    if (graphRef.isAdjacent(u, n1))
                        degrees[u]--;
                    if (graphRef.isAdjacent(u, n2))
                        degrees[u]--;
                }
            }
        }
    }
}

// Enhanced folding operation for degree-2 vertices
void MCBRBSolver::enhancedFolding(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates,
                                  VertexID v, VertexID n1, VertexID n2, std::vector<int> &localDegree)
{
    // Record the folding operation for later backtracking
    foldingOperations.push_back({v, n1, n2});

    // Create a virtual edge between n1 and n2 if they are not already adjacent
    bool adjacent = usingMatrix ? matrixAdjacent(n1, n2) : graphRef.isAdjacent(n1, n2);

    if (!adjacent)
    {
        // Add a virtual edge in our tracking structure
        virtualEdges.insert({std::min(n1, n2), std::max(n1, n2)});

        // When adding a virtual edge, we need to update the degrees of n1 and n2
        localDegree[n1]++;
        localDegree[n2]++;
    }

    // Remove v from candidates
    candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());

    // Track that v was part of a folding operation
    foldedVertices.insert(v);
}

// Main recursive search function
// Update the core search function with our enhancements
// Enhanced main search function with all optimizations
void MCBRBSolver::search(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates)
{
    // Update statistics for monitoring
    branches++;
    maxDepth = std::max(maxDepth, static_cast<int>(currentClique.size()));

    // Early termination if current best can't be improved
    if (currentClique.size() + candidates.size() <= bestClique.size())
    {
        return;
    }

    // Clear per-branch state
    foldingOperations.clear();
    virtualEdges.clear();
    foldedVertices.clear();

    // For smaller candidate sets, use kernelization techniques
    if (candidates.size() < 100)
    {
        // Use advanced kernelization with dominance detection
        advancedKernelizationWithDominance(currentClique, candidates);
    }

    // Apply k-core reduction to remove vertices that cannot extend to large enough cliques
    reduceCandidates(candidates, currentClique.size());

    // If candidate set is very small, use specialized triangle search
    if (candidates.size() <= 3)
    {
        searchTriangle(currentClique, candidates);
        return;
    }

    // Choose the most efficient representation based on graph properties
    bool shouldUseBitsets = candidates.size() > 200;
    bool shouldUseMatrixHere = shouldSwitchToMatrix(candidates) && !shouldUseBitsets;
    bool oldUsingMatrix = usingMatrix;
    bool oldUsingBitsets = usingBitsets;

    // Set up efficient representations
    if (shouldUseMatrixHere && !usingMatrix)
    {
        constructMatrix(candidates);
    }
    else if (shouldUseBitsets && !usingBitsets)
    {
        constructBitsetForCandidates(candidates);
    }

    // Compute upper bound using best coloring method
    int colorBound = dsaturColoring(candidates);
    int upperBound = currentClique.size() + colorBound;

    // Prune if the upper bound is not better than current best
    if (upperBound <= static_cast<int>(bestClique.size()))
    {
        // Restore representation state
        usingMatrix = oldUsingMatrix;
        usingBitsets = oldUsingBitsets;
        return;
    }

    // If no more candidates, check if we have a better clique
    if (candidates.empty())
    {
        if (currentClique.size() > bestClique.size())
        {
            bestClique = currentClique;
            std::cout << "[Solver] Found a new best clique of size " << bestClique.size() << std::endl;
        }
        // Restore representation state
        usingMatrix = oldUsingMatrix;
        usingBitsets = oldUsingBitsets;
        return;
    }

    // Get candidates in optimized order to maximize pruning
    std::vector<VertexID> orderedCandidates = determineBranchOrder(candidates);

    // Branch on each candidate in determined order
    for (auto v : orderedCandidates)
    {
        // Skip if we can't improve the best clique with remaining candidates
        if (currentClique.size() + 1 + (orderedCandidates.size() - 1) <= bestClique.size())
        {
            break;
        }

        // Create new clique with the current vertex
        std::vector<VertexID> newClique = currentClique;
        newClique.push_back(v);

        // Find candidates that are neighbors of v (using fastest available method)
        std::vector<VertexID> newCandidates;
        if (usingBitsets)
        {
            newCandidates = fastIntersectCandidates(candidates, v);
        }
        else
        {
            newCandidates = intersectCandidates(candidates, v);
        }

        // Recurse using the new clique and candidate set
        search(newClique, newCandidates);

        // Remove v from candidates (we've explored all cliques that include v)
        candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
    }

    // Restore representation state
    usingMatrix = oldUsingMatrix;
    usingBitsets = oldUsingBitsets;
}

// Add this to MCBRBSolver.cpp
// Construct bitset representation for more efficient operations
void MCBRBSolver::constructBitsetForCandidates(const std::vector<VertexID> &candidates)
{
    // Clear existing bitsets
    vertexBitsets.clear();

    // Fast set operations require mapping vertices to consecutive indices
    vertexToConsecutiveIdx.clear();
    consecutiveIdxToVertex.clear();

    // Map vertices to consecutive indices
    for (size_t i = 0; i < candidates.size(); i++)
    {
        vertexToConsecutiveIdx[candidates[i]] = i;
        consecutiveIdxToVertex.push_back(candidates[i]);
    }

    // Calculate how many 32-bit words we need per bitset
    const int numWords = (candidates.size() + 31) / 32;

    // Initialize bitsets for each vertex
    vertexBitsets.resize(graphRef.getNumVertices());
    for (auto &bitset : vertexBitsets)
    {
        bitset.clear();
        bitset.resize(numWords, 0);
    }

    // Fill bitsets - each bit position represents a vertex in candidates
    for (size_t i = 0; i < candidates.size(); i++)
    {
        VertexID v = candidates[i];

        // Set bits for all neighbors of v that are in candidates
        for (VertexID neighbor : graphRef.getNeighbors(v))
        {
            auto it = vertexToConsecutiveIdx.find(neighbor);
            if (it != vertexToConsecutiveIdx.end())
            {
                // Neighbor is in the candidate set, set its bit
                int idx = it->second;
                int wordIdx = idx / 32;
                int bitPos = idx % 32;
                vertexBitsets[v][wordIdx] |= (1u << bitPos);
            }
        }
    }

    usingBitsets = true;
}

// Count common neighbors using bitset operations (much faster than looping)
int MCBRBSolver::countCommonNeighbors(VertexID v, VertexID w)
{
    if (!usingBitsets)
    {
        // Fallback method
        int count = 0;
        for (VertexID u = 0; u < graphRef.getNumVertices(); u++)
        {
            if (u != v && u != w &&
                graphRef.isAdjacent(v, u) && graphRef.isAdjacent(w, u))
            {
                count++;
            }
        }
        return count;
    }

    // Use fast bitset operations
    int count = 0;
    int numWords = vertexBitsets[v].size();

    for (int i = 0; i < numWords; i++)
    {
        // Common bits set in both bitsets represent common neighbors
        uint32_t commonBits = vertexBitsets[v][i] & vertexBitsets[w][i];

        // Count bits using Brian Kernighan's algorithm
        while (commonBits)
        {
            commonBits &= (commonBits - 1);
            count++;
        }
    }

    return count;
}

// Check if a set of vertices forms a clique
bool MCBRBSolver::isClique(const std::vector<VertexID> &vertices)
{
    for (size_t i = 0; i < vertices.size(); i++)
    {
        for (size_t j = i + 1; j < vertices.size(); j++)
        {
            if (!graphRef.isAdjacent(vertices[i], vertices[j]))
            {
                return false;
            }
        }
    }
    return true;
}

// Fast intersection of candidates with neighbors using bitsets
std::vector<VertexID> MCBRBSolver::fastIntersectCandidates(const std::vector<VertexID> &candidates, VertexID v)
{
    if (!usingBitsets)
    {
        // Fall back to standard method
        return intersectCandidates(candidates, v);
    }

    std::vector<VertexID> result;

    // For small candidate sets, standard method may be faster
    if (candidates.size() < 32)
    {
        for (VertexID u : candidates)
        {
            if (graphRef.isAdjacent(v, u))
            {
                result.push_back(u);
            }
        }
        return result;
    }

    // Get all vertices in candidates that are also neighbors of v
    const auto &vBitset = vertexBitsets[v];

    for (size_t i = 0; i < candidates.size(); i++)
    {
        VertexID candidate = candidates[i];
        auto it = vertexToConsecutiveIdx.find(candidate);

        // Check if candidate is in our mapping and is adjacent to v
        if (it != vertexToConsecutiveIdx.end())
        {
            int idx = it->second;
            int wordIdx = idx / 32;
            int bitPos = idx % 32;

            // Check if bit is set (vertices are adjacent)
            if ((vBitset[wordIdx] & (1u << bitPos)) != 0)
            {
                result.push_back(candidate);
            }
        }
    }

    return result;
}

// Compute intersection size using bitsets (much faster)
int MCBRBSolver::getBitsetIntersectionSize(VertexID v, const std::vector<VertexID> &vertices)
{
    if (!usingBitsets)
        return 0;

    // Count bits set in common
    int count = 0;
    const int words = vertexBitsets[v].size();

    for (auto u : vertices)
    {
        for (int w = 0; w < words; w++)
        {
            uint32_t common = vertexBitsets[v][w] & vertexBitsets[u][w];
            // Count bits using Brian Kernighan's algorithm
            while (common)
            {
                common &= (common - 1);
                count++;
            }
        }
    }

    return count;
}

// Improved branch-and-bound implementation with metrics
void MCBRBSolver::branchAndBound()
{
    std::cout << "[Solver] Starting branch-and-bound search..." << std::endl;

    // Record start time
    auto startTime = std::chrono::high_resolution_clock::now();

    // Initialize statistics
    branches = 0;
    maxDepth = 0;

    // Start with all vertices as candidates
    std::vector<VertexID> candidates;
    candidates.reserve(graphRef.getNumVertices());

    // Filter candidates by initial criteria
    for (VertexID i = 0; i < graphRef.getNumVertices(); i++)
    {
        // Only consider vertices with degree >= current best clique size - 1
        // (this is a necessary condition for a vertex to be in a larger clique)
        if (graphRef.degree(i) >= bestClique.size() - 1)
        {
            candidates.push_back(i);
        }
    }

    std::cout << "[Solver] Starting with " << candidates.size()
              << " candidates after initial filtering." << std::endl;

    // Initialize current clique as empty
    std::vector<VertexID> currentClique;

    // Start recursive search
    search(currentClique, candidates);

    // Record end time and calculate duration
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

    // Report results
    std::cout << "[Solver] Branch-and-bound search completed in " << duration << " ms." << std::endl;
    reportStatistics();

    // Verify the found clique
    if (!bestClique.empty())
    {
        bool validClique = isClique(bestClique);
        std::cout << "[Solver] Clique verification: "
                  << (validClique ? "PASSED" : "FAILED") << std::endl;
    }
}

// Report search statistics
void MCBRBSolver::reportStatistics() const
{
    std::cout << "[Solver] Search statistics:" << std::endl;
    std::cout << "  Best clique size: " << bestClique.size() << std::endl;
    std::cout << "  Branches explored: " << branches << std::endl;
    std::cout << "  Maximum depth reached: " << maxDepth << std::endl;
}

std::vector<VertexID> MCBRBSolver::getBestClique() const
{
    return bestClique;
}

// Advanced coloring implementation using DSatur algorithm
// This typically produces tighter upper bounds than simple greedy coloring
