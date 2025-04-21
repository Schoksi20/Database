#include "MCBRBSolver.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <cstring>
#include <chrono>

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

void MCBRBSolver::initialize()
{
    auto cliqueDD = HelperAlgorithms::runMCDD(graphRef);
    auto cliqueEGO = HelperAlgorithms::runImprovedMCEGO(graphRef);

    if (cliqueDD.size() >= cliqueEGO.size())
        bestClique = cliqueDD;
    else
        bestClique = cliqueEGO;

    std::cout << "[Solver] Initial best clique size: " << bestClique.size() << std::endl;
}

std::vector<VertexID> MCBRBSolver::determineBranchOrder(const std::vector<VertexID> &candidates)
{

    if (candidates.size() > 100)
    {
        return computeDegeneracyOrdering(candidates);
    }

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

        int degPosition = std::find(graphRef.getDegeneracyOrder().begin(),
                                    graphRef.getDegeneracyOrder().end(), v) -
                          graphRef.getDegeneracyOrder().begin();

        double score = 0.7 * internalDegree + 0.3 * (graphRef.getNumVertices() - degPosition);

        scoredCandidates.push_back(std::make_pair(score, v));
    }

    std::sort(scoredCandidates.begin(), scoredCandidates.end(),
              [](const std::pair<int, VertexID> &a, const std::pair<int, VertexID> &b)
              { return a.first > b.first; });

    std::vector<VertexID> orderedCandidates;
    for (const auto &pair : scoredCandidates)
    {
        orderedCandidates.push_back(pair.second);
    }

    return orderedCandidates;
}

void MCBRBSolver::advancedKernelizationWithDominance(std::vector<VertexID> &currentClique,
                                                     std::vector<VertexID> &candidates)
{

    advancedKernelization(currentClique, candidates);

    if (candidates.size() <= 3)
        return;

    std::vector<bool> dominated(graphRef.getNumVertices(), false);

    for (size_t i = 0; i < candidates.size(); i++)
    {
        VertexID v = candidates[i];
        if (dominated[v])
            continue;

        std::vector<VertexID> vNeighbors;
        for (size_t j = 0; j < candidates.size(); j++)
        {
            VertexID u = candidates[j];
            if (u != v && (usingMatrix ? matrixAdjacent(v, u) : graphRef.isAdjacent(v, u)))
            {
                vNeighbors.push_back(u);
            }
        }

        for (size_t j = i + 1; j < candidates.size(); j++)
        {
            VertexID w = candidates[j];
            if (dominated[w])
                continue;

            bool adjacent = (usingMatrix ? matrixAdjacent(v, w) : graphRef.isAdjacent(v, w));
            if (!adjacent)
                continue;

            std::vector<VertexID> wNeighbors;
            for (size_t k = 0; k < candidates.size(); k++)
            {
                VertexID u = candidates[k];
                if (u != w && (usingMatrix ? matrixAdjacent(w, u) : graphRef.isAdjacent(w, u)))
                {
                    wNeighbors.push_back(u);
                }
            }

            bool vDominatesW = true;
            for (VertexID u : wNeighbors)
            {
                if (u == v)
                    continue;

                if (std::find(vNeighbors.begin(), vNeighbors.end(), u) == vNeighbors.end())
                {
                    vDominatesW = false;
                    break;
                }
            }

            bool wDominatesV = true;
            for (VertexID u : vNeighbors)
            {
                if (u == w)
                    continue;

                if (std::find(wNeighbors.begin(), wNeighbors.end(), u) == wNeighbors.end())
                {
                    wDominatesV = false;
                    break;
                }
            }

            if (vDominatesW)
            {
                dominated[w] = true;
            }
            else if (wDominatesV)
            {
                dominated[v] = true;
                break;
            }
        }
    }

    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(),
                       [&dominated](VertexID v)
                       { return dominated[v]; }),
        candidates.end());
}

bool MCBRBSolver::shouldSwitchToMatrix(const std::vector<VertexID> &candidates)
{
    const int sizeThreshold = 50;
    const double densityThreshold = 0.4;

    if (candidates.size() < sizeThreshold)
    {
        return false;
    }

    int sampleSize = std::min(100, static_cast<int>(candidates.size()));
    int edges = 0;
    int possibleEdges = sampleSize * (sampleSize - 1) / 2;

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

int MCBRBSolver::advancedColoring(const std::vector<VertexID> &candidates)
{
    int n = candidates.size();
    if (n == 0)
        return 0;

    std::vector<int> color(n, -1);

    std::vector<int> saturation(n, 0);

    std::vector<std::vector<bool>> usedColors(n);
    for (int i = 0; i < n; i++)
    {
        usedColors[i].resize(n, false);
    }

    for (int i = 0; i < n; i++)
    {

        int maxSat = -1;
        int maxDeg = -1;
        int maxIdx = -1;

        for (int j = 0; j < n; j++)
        {
            if (color[j] == -1)
            {
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

        int c = 0;
        while (c < n && usedColors[v][c])
            c++;

        color[v] = c;

        for (int j = 0; j < n; j++)
        {
            if (color[j] == -1)
            {
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

    int maxColor = *std::max_element(color.begin(), color.end());
    return maxColor + 1;
}

int MCBRBSolver::greedyColoring(const std::vector<VertexID> &origCandidates)
{
    std::vector<VertexID> candidates = origCandidates;
    int n = candidates.size();
    if (n == 0)
        return 0;

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

    auto cmpLocalDegree = [&localDegree](int a, int b)
    {
        return localDegree[a] > localDegree[b];
    };

    std::vector<int> ordering(n);
    for (int i = 0; i < n; i++)
        ordering[i] = i;
    std::sort(ordering.begin(), ordering.end(), cmpLocalDegree);

    std::vector<VertexID> orderedCandidates(n);
    for (int i = 0; i < n; i++)
    {
        orderedCandidates[i] = candidates[ordering[i]];
    }

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
    return maxColor + 1;
}

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

void MCBRBSolver::reduceCandidates(std::vector<VertexID> &candidates, int currentCliqueSize)
{
    int threshold = bestClique.size() - currentCliqueSize;
    if (threshold <= 0)
        return;

    bool removed;
    std::vector<int> localDegree(graphRef.getNumVertices(), 0);

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

bool MCBRBSolver::shouldUseMatrix(const std::vector<VertexID> &candidates) const
{

    const int minSize = 20;
    const double densityThreshold = 0.3;

    if (candidates.size() < minSize)
        return false;

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

void MCBRBSolver::constructMatrix(const std::vector<VertexID> &candidates)
{
    matrixSize = candidates.size();
    if (matrixSize == 0)
    {
        usingMatrix = false;
        return;
    }

    matrixMapping.resize(graphRef.getNumVertices(), -1);
    inverseMapping.resize(matrixSize);

    for (int i = 0; i < matrixSize; i++)
    {
        VertexID v = candidates[i];
        matrixMapping[v] = i;
        inverseMapping[i] = v;
    }

    int matrixBytes = (matrixSize * matrixSize + 7) / 8;
    adjacencyMatrix.resize(matrixBytes, 0);

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

void MCBRBSolver::setMatrixBit(int row, int col)
{
    int bitPos = row * matrixSize + col;
    adjacencyMatrix[bitPos / 8] |= (1 << (bitPos % 8));
}

bool MCBRBSolver::getMatrixBit(int row, int col) const
{
    int bitPos = row * matrixSize + col;
    return (adjacencyMatrix[bitPos / 8] & (1 << (bitPos % 8))) != 0;
}

bool MCBRBSolver::matrixAdjacent(VertexID u, VertexID v) const
{
    if (!usingMatrix)
        return graphRef.isAdjacent(u, v);

    int rowU = matrixMapping[u];
    int rowV = matrixMapping[v];

    if (rowU == -1 || rowV == -1)
        return graphRef.isAdjacent(u, v);

    return getMatrixBit(rowU, rowV);
}

std::vector<VertexID> MCBRBSolver::computeDegeneracyOrdering(const std::vector<VertexID> &origCandidates)
{

    std::vector<VertexID> candidates = origCandidates;
    int n = candidates.size();
    std::vector<VertexID> order;

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

    std::vector<bool> removed(n, false);
    for (int i = 0; i < n; i++)
    {

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
            break;

        order.push_back(candidates[minIdx]);
        removed[minIdx] = true;

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

void MCBRBSolver::searchTriangle(const std::vector<VertexID> &currentClique, const std::vector<VertexID> &candidates)
{
    int cliqueSize = currentClique.size();
    int candSize = candidates.size();

    if (cliqueSize + candSize <= bestClique.size())
    {
        return;
    }

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
                        return;
                    }
                }
            }
        }
    }

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

void MCBRBSolver::kernelizeByDegreeReduction(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates)
{
    bool changed = true;

    std::vector<int> localDegree(graphRef.getNumVertices(), 0);

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

        for (auto v : degreeZero)
        {

            currentClique.push_back(v);
            candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
        }

        for (auto v : degreeOne)
        {
            if (std::find(candidates.begin(), candidates.end(), v) == candidates.end())
                continue;

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

                currentClique.push_back(v);

                candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
                candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor), candidates.end());

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

        for (auto v : degreeTwo)
        {
            if (std::find(candidates.begin(), candidates.end(), v) == candidates.end())
                continue;

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

                    currentClique.push_back(v);
                    currentClique.push_back(neighbor1);
                    currentClique.push_back(neighbor2);

                    candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor1), candidates.end());
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor2), candidates.end());

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

                    currentClique.push_back(v);

                    candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());

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

    std::vector<int> color(n, -1);

    std::vector<int> saturation(n, 0);

    std::vector<int> degrees(n, 0);

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

    std::vector<std::vector<bool>> usedColors(n);
    for (int i = 0; i < n; i++)
    {
        usedColors[i].resize(n, false);
    }

    for (int c = 0; c < n; c++)
    {

        int maxSat = -1;
        int maxDeg = -1;
        int maxIdx = -1;

        for (int i = 0; i < n; i++)
        {
            if (color[i] == -1)
            {
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
            break;

        int minColor = 0;
        while (minColor < n && usedColors[maxIdx][minColor])
            minColor++;

        color[maxIdx] = minColor;

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

    return *std::max_element(color.begin(), color.end()) + 1;
}

bool MCBRBSolver::shouldUseBitMatrix(const std::vector<VertexID> &candidates)
{

    if (candidates.size() < 32)
        return false;

    int edges = 0;
    int totalPossible = candidates.size() * (candidates.size() - 1) / 2;

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
    return sampleDensity > 0.3;
}

void MCBRBSolver::constructBitMatrix(const std::vector<VertexID> &candidates)
{
    int n = candidates.size();

    vertexToMatrixIndex.clear();
    matrixToVertex.resize(n);

    for (int i = 0; i < n; i++)
    {
        vertexToMatrixIndex[candidates[i]] = i;
        matrixToVertex[i] = candidates[i];
    }

    int words = (n + 31) / 32;
    bitMatrix.resize(n * words, 0);

    for (int i = 0; i < n; i++)
    {
        VertexID v = candidates[i];
        for (int j = 0; j < n; j++)
        {
            if (i != j && graphRef.isAdjacent(v, candidates[j]))
            {

                int wordIndex = i * words + (j / 32);
                int bitPosition = j % 32;
                bitMatrix[wordIndex] |= (1U << bitPosition);
            }
        }
    }

    usingBitMatrix = true;
}

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

void MCBRBSolver::advancedKernelization(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates)
{
    if (candidates.size() <= 1)
        return;

    bool changes = true;
    std::vector<int> degrees(graphRef.getNumVertices(), 0);

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

        if (!degree0.empty())
        {
            changes = true;

            VertexID v = degree0[0];
            currentClique.push_back(v);

            candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
        }

        else if (!degree1.empty())
        {
            changes = true;
            VertexID v = degree1[0];

            VertexID neighbor = -1;
            for (VertexID u : candidates)
            {
                if (u != v && graphRef.isAdjacent(v, u))
                {
                    neighbor = u;
                    break;
                }
            }

            currentClique.push_back(v);

            candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
            if (neighbor != -1)
            {
                candidates.erase(std::remove(candidates.begin(), candidates.end(), neighbor), candidates.end());

                for (VertexID u : candidates)
                {
                    if (graphRef.isAdjacent(u, neighbor))
                    {
                        degrees[u]--;
                    }
                }
            }
        }

        else if (!degree2.empty())
        {
            changes = true;
            VertexID v = degree2[0];

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

            if (n1 != -1 && n2 != -1 && graphRef.isAdjacent(n1, n2))
            {

                currentClique.push_back(v);

                std::vector<VertexID> toRemove = {v, n1, n2};
                for (VertexID x : toRemove)
                {
                    candidates.erase(std::remove(candidates.begin(), candidates.end(), x), candidates.end());

                    for (VertexID u : candidates)
                    {
                        if (graphRef.isAdjacent(u, x))
                        {
                            degrees[u]--;
                        }
                    }
                }
            }

            else if (n1 != -1 && n2 != -1)
            {

                candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());

                currentClique.push_back(v);

                candidates.erase(std::remove(candidates.begin(), candidates.end(), n1), candidates.end());
                candidates.erase(std::remove(candidates.begin(), candidates.end(), n2), candidates.end());

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

void MCBRBSolver::enhancedFolding(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates,
                                  VertexID v, VertexID n1, VertexID n2, std::vector<int> &localDegree)
{

    foldingOperations.push_back({v, n1, n2});

    bool adjacent = usingMatrix ? matrixAdjacent(n1, n2) : graphRef.isAdjacent(n1, n2);

    if (!adjacent)
    {

        virtualEdges.insert({std::min(n1, n2), std::max(n1, n2)});

        localDegree[n1]++;
        localDegree[n2]++;
    }

    candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());

    foldedVertices.insert(v);
}

void MCBRBSolver::search(std::vector<VertexID> &currentClique, std::vector<VertexID> &candidates)
{

    branches++;
    maxDepth = std::max(maxDepth, static_cast<int>(currentClique.size()));

    if (currentClique.size() + candidates.size() <= bestClique.size())
    {
        return;
    }

    foldingOperations.clear();
    virtualEdges.clear();
    foldedVertices.clear();

    if (candidates.size() < 100)
    {

        advancedKernelizationWithDominance(currentClique, candidates);
    }

    reduceCandidates(candidates, currentClique.size());

    if (candidates.size() <= 3)
    {
        searchTriangle(currentClique, candidates);
        return;
    }

    bool shouldUseBitsets = candidates.size() > 200;
    bool shouldUseMatrixHere = shouldSwitchToMatrix(candidates) && !shouldUseBitsets;
    bool oldUsingMatrix = usingMatrix;
    bool oldUsingBitsets = usingBitsets;

    if (shouldUseMatrixHere && !usingMatrix)
    {
        constructMatrix(candidates);
    }
    else if (shouldUseBitsets && !usingBitsets)
    {
        constructBitsetForCandidates(candidates);
    }

    int colorBound = dsaturColoring(candidates);
    int upperBound = currentClique.size() + colorBound;

    if (upperBound <= static_cast<int>(bestClique.size()))
    {

        usingMatrix = oldUsingMatrix;
        usingBitsets = oldUsingBitsets;
        return;
    }

    if (candidates.empty())
    {
        if (currentClique.size() > bestClique.size())
        {
            bestClique = currentClique;
            std::cout << "[Solver] Found a new best clique of size " << bestClique.size() << std::endl;
        }

        usingMatrix = oldUsingMatrix;
        usingBitsets = oldUsingBitsets;
        return;
    }

    std::vector<VertexID> orderedCandidates = determineBranchOrder(candidates);

    for (auto v : orderedCandidates)
    {

        if (currentClique.size() + 1 + (orderedCandidates.size() - 1) <= bestClique.size())
        {
            break;
        }

        std::vector<VertexID> newClique = currentClique;
        newClique.push_back(v);

        std::vector<VertexID> newCandidates;
        if (usingBitsets)
        {
            newCandidates = fastIntersectCandidates(candidates, v);
        }
        else
        {
            newCandidates = intersectCandidates(candidates, v);
        }

        search(newClique, newCandidates);

        candidates.erase(std::remove(candidates.begin(), candidates.end(), v), candidates.end());
    }

    usingMatrix = oldUsingMatrix;
    usingBitsets = oldUsingBitsets;
}

void MCBRBSolver::constructBitsetForCandidates(const std::vector<VertexID> &candidates)
{

    vertexBitsets.clear();

    vertexToConsecutiveIdx.clear();
    consecutiveIdxToVertex.clear();

    for (size_t i = 0; i < candidates.size(); i++)
    {
        vertexToConsecutiveIdx[candidates[i]] = i;
        consecutiveIdxToVertex.push_back(candidates[i]);
    }

    const int numWords = (candidates.size() + 31) / 32;

    vertexBitsets.resize(graphRef.getNumVertices());
    for (auto &bitset : vertexBitsets)
    {
        bitset.clear();
        bitset.resize(numWords, 0);
    }

    for (size_t i = 0; i < candidates.size(); i++)
    {
        VertexID v = candidates[i];

        for (VertexID neighbor : graphRef.getNeighbors(v))
        {
            auto it = vertexToConsecutiveIdx.find(neighbor);
            if (it != vertexToConsecutiveIdx.end())
            {

                int idx = it->second;
                int wordIdx = idx / 32;
                int bitPos = idx % 32;
                vertexBitsets[v][wordIdx] |= (1u << bitPos);
            }
        }
    }

    usingBitsets = true;
}

int MCBRBSolver::countCommonNeighbors(VertexID v, VertexID w)
{
    if (!usingBitsets)
    {

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

    int count = 0;
    int numWords = vertexBitsets[v].size();

    for (int i = 0; i < numWords; i++)
    {

        uint32_t commonBits = vertexBitsets[v][i] & vertexBitsets[w][i];

        while (commonBits)
        {
            commonBits &= (commonBits - 1);
            count++;
        }
    }

    return count;
}

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

std::vector<VertexID> MCBRBSolver::fastIntersectCandidates(const std::vector<VertexID> &candidates, VertexID v)
{
    if (!usingBitsets)
    {

        return intersectCandidates(candidates, v);
    }

    std::vector<VertexID> result;

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

    const auto &vBitset = vertexBitsets[v];

    for (size_t i = 0; i < candidates.size(); i++)
    {
        VertexID candidate = candidates[i];
        auto it = vertexToConsecutiveIdx.find(candidate);

        if (it != vertexToConsecutiveIdx.end())
        {
            int idx = it->second;
            int wordIdx = idx / 32;
            int bitPos = idx % 32;

            if ((vBitset[wordIdx] & (1u << bitPos)) != 0)
            {
                result.push_back(candidate);
            }
        }
    }

    return result;
}

int MCBRBSolver::getBitsetIntersectionSize(VertexID v, const std::vector<VertexID> &vertices)
{
    if (!usingBitsets)
        return 0;

    int count = 0;
    const int words = vertexBitsets[v].size();

    for (auto u : vertices)
    {
        for (int w = 0; w < words; w++)
        {
            uint32_t common = vertexBitsets[v][w] & vertexBitsets[u][w];

            while (common)
            {
                common &= (common - 1);
                count++;
            }
        }
    }

    return count;
}

void MCBRBSolver::branchAndBound()
{
    std::cout << "[Solver] Starting branch-and-bound search..." << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    branches = 0;
    maxDepth = 0;

    std::vector<VertexID> candidates;
    candidates.reserve(graphRef.getNumVertices());

    for (VertexID i = 0; i < graphRef.getNumVertices(); i++)
    {

        if (graphRef.degree(i) >= bestClique.size() - 1)
        {
            candidates.push_back(i);
        }
    }

    std::cout << "[Solver] Starting with " << candidates.size()
              << " candidates after initial filtering." << std::endl;

    std::vector<VertexID> currentClique;

    search(currentClique, candidates);

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

    std::cout << "[Solver] Branch-and-bound search completed in " << duration << " ms." << std::endl;
    reportStatistics();

    if (!bestClique.empty())
    {
        bool validClique = isClique(bestClique);
        std::cout << "[Solver] Clique verification: "
                  << (validClique ? "PASSED" : "FAILED") << std::endl;
    }
}

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
