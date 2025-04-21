#include "HelperAlgorithms.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <unordered_set>

Clique HelperAlgorithms::runMCDD(const Graph &graph)
{
    std::cout << "[Helper] Running MC-DD algorithm..." << std::endl;

    Clique maximalClique = findMaximalCliqueByDegeneracy(graph);

    std::cout << "[Helper] MC-DD returns a clique of size " << maximalClique.size() << std::endl;
    return maximalClique;
}

Clique HelperAlgorithms::runMCEGO(const Graph &graph)
{
    std::cout << "[Helper] Running MC-EGO algorithm..." << std::endl;

    Clique bestClique;

    for (VertexID v = 0; v < graph.getNumVertices(); v++)
    {

        if (graph.degree(v) < bestClique.size())
        {
            continue;
        }

        Clique egoClique = exploreEgoNetwork(graph, v);

        if (egoClique.size() > bestClique.size())
        {
            bestClique = egoClique;
        }
    }

    std::cout << "[Helper] MC-EGO returns a clique of size " << bestClique.size() << std::endl;
    return bestClique;
}

Clique HelperAlgorithms::findMaximalCliqueByDegeneracy(const Graph &graph)
{

    std::vector<VertexID> degeneracyOrder;
    std::vector<int> degrees(graph.getNumVertices());
    std::vector<bool> removed(graph.getNumVertices(), false);

    for (VertexID v = 0; v < graph.getNumVertices(); v++)
    {
        degrees[v] = graph.degree(v);
    }

    for (VertexID i = 0; i < graph.getNumVertices(); i++)
    {

        VertexID minDegVertex = 0;
        int minDeg = graph.getNumVertices() + 1;

        for (VertexID v = 0; v < graph.getNumVertices(); v++)
        {
            if (!removed[v] && degrees[v] < minDeg)
            {
                minDeg = degrees[v];
                minDegVertex = v;
            }
        }

        degeneracyOrder.push_back(minDegVertex);

        removed[minDegVertex] = true;

        for (VertexID neighbor : graph.getNeighbors(minDegVertex))
        {
            if (!removed[neighbor])
            {
                degrees[neighbor]--;
            }
        }
    }

    Clique maximalClique;
    std::vector<bool> inClique(graph.getNumVertices(), false);

    for (auto it = degeneracyOrder.rbegin(); it != degeneracyOrder.rend(); ++it)
    {
        VertexID v = *it;

        bool canAdd = true;
        for (VertexID u : maximalClique)
        {
            if (!graph.isAdjacent(v, u))
            {
                canAdd = false;
                break;
            }
        }

        if (canAdd)
        {
            maximalClique.push_back(v);
            inClique[v] = true;
        }
    }

    return maximalClique;
}

Clique HelperAlgorithms::runImprovedMCEGO(const Graph &graph)
{
    std::cout << "[Helper] Running Improved MC-EGO algorithm..." << std::endl;

    Clique bestClique;

    std::vector<VertexID> degOrder;
    std::vector<int> degrees(graph.getNumVertices());
    std::vector<bool> removed(graph.getNumVertices(), false);

    for (VertexID v = 0; v < graph.getNumVertices(); v++)
    {
        degrees[v] = graph.degree(v);
    }

    for (VertexID i = 0; i < graph.getNumVertices(); i++)
    {

        VertexID minDegVertex = 0;
        int minDeg = graph.getNumVertices() + 1;

        for (VertexID v = 0; v < graph.getNumVertices(); v++)
        {
            if (!removed[v] && degrees[v] < minDeg)
            {
                minDeg = degrees[v];
                minDegVertex = v;
            }
        }

        degOrder.push_back(minDegVertex);

        removed[minDegVertex] = true;

        for (VertexID neighbor : graph.getNeighbors(minDegVertex))
        {
            if (!removed[neighbor])
            {
                degrees[neighbor]--;
            }
        }
    }

    std::reverse(degOrder.begin(), degOrder.end());

    for (VertexID v : degOrder)
    {

        if (graph.degree(v) < bestClique.size())
        {
            continue;
        }

        Clique egoClique = improvedExploreEgoNetwork(graph, v);

        if (egoClique.size() > bestClique.size())
        {
            bestClique = egoClique;
            std::cout << "[Helper] MC-EGO found improved clique of size "
                      << bestClique.size() << std::endl;
        }
    }

    std::cout << "[Helper] Improved MC-EGO returns a clique of size " << bestClique.size() << std::endl;
    return bestClique;
}

Clique HelperAlgorithms::improvedExploreEgoNetwork(const Graph &graph, VertexID center)
{
    const auto &neighbors = graph.getNeighbors(center);

    if (neighbors.empty())
    {
        Clique singletonClique;
        singletonClique.push_back(center);
        return singletonClique;
    }

    std::vector<VertexID> candidates(neighbors.begin(), neighbors.end());

    double density = 0.0;
    int edgeCount = 0;
    int possibleEdges = candidates.size() * (candidates.size() - 1) / 2;

    if (possibleEdges > 0)
    {
        for (size_t i = 0; i < candidates.size(); i++)
        {
            for (size_t j = i + 1; j < candidates.size(); j++)
            {
                if (graph.isAdjacent(candidates[i], candidates[j]))
                {
                    edgeCount++;
                }
            }
        }
        density = static_cast<double>(edgeCount) / possibleEdges;
    }

    if (density > 0.5 && candidates.size() > 10)
    {
        return findMaximalCliqueInDenseNeighborhood(graph, center, candidates);
    }

    std::vector<int> localDegree(graph.getNumVertices(), 0);

    for (auto v : candidates)
    {
        for (auto u : candidates)
        {
            if (v != u && graph.isAdjacent(v, u))
            {
                localDegree[v]++;
            }
        }
    }

    int maxK = 0;
    for (auto v : candidates)
    {
        maxK = std::max(maxK, localDegree[v]);
    }

    int k = std::max(1, maxK / 2);
    bool reduced = true;

    while (reduced)
    {
        reduced = false;
        std::vector<VertexID> newCandidates;

        for (auto v : candidates)
        {
            if (localDegree[v] >= k)
            {
                newCandidates.push_back(v);
            }
            else
            {
                reduced = true;

                for (auto u : candidates)
                {
                    if (u != v && graph.isAdjacent(u, v))
                    {
                        localDegree[u]--;
                    }
                }
            }
        }

        if (newCandidates.empty() && k > 1)
        {
            k = k / 2;
            reduced = true;
            newCandidates = candidates;
        }
        else
        {
            candidates = newCandidates;

            if (!reduced && !candidates.empty() && k < maxK)
            {
                k++;
                reduced = true;
            }
        }
    }

    std::vector<VertexID> degeneracyOrder;
    std::vector<bool> vertexRemoved(graph.getNumVertices(), false);

    while (!candidates.empty())
    {
        VertexID minDegVertex = candidates[0];
        VertexID minDeg = localDegree[minDegVertex];
        int minIdx = 0;

        for (size_t i = 1; i < candidates.size(); i++)
        {
            if (localDegree[candidates[i]] < minDeg)
            {
                minDeg = localDegree[candidates[i]];
                minDegVertex = candidates[i];
                minIdx = i;
            }
        }

        degeneracyOrder.push_back(minDegVertex);
        vertexRemoved[minDegVertex] = true;

        candidates[minIdx] = candidates.back();
        candidates.pop_back();

        for (auto v : candidates)
        {
            if (graph.isAdjacent(v, minDegVertex))
            {
                localDegree[v]--;
            }
        }
    }

    Clique maximalClique;

    for (auto it = degeneracyOrder.rbegin(); it != degeneracyOrder.rend(); ++it)
    {
        VertexID v = *it;

        bool canAdd = true;
        for (auto u : maximalClique)
        {
            if (!graph.isAdjacent(v, u))
            {
                canAdd = false;
                break;
            }
        }

        if (canAdd)
        {
            maximalClique.push_back(v);
        }
    }

    maximalClique.push_back(center);

    bool extended;
    do
    {
        extended = false;
        for (auto n : neighbors)
        {

            if (std::find(maximalClique.begin(), maximalClique.end(), n) != maximalClique.end())
            {
                continue;
            }

            bool canAdd = true;
            for (auto u : maximalClique)
            {
                if (!graph.isAdjacent(n, u))
                {
                    canAdd = false;
                    break;
                }
            }

            if (canAdd)
            {
                maximalClique.push_back(n);
                extended = true;
                break;
            }
        }
    } while (extended);

    return maximalClique;
}

Clique HelperAlgorithms::findMaximalCliqueInDenseNeighborhood(const Graph &graph, VertexID center, const std::vector<VertexID> &neighbors)
{

    Clique bestClique;
    bestClique.push_back(center);

    std::vector<VertexID> startVertices;

    if (neighbors.size() > 20)
    {

        std::vector<std::pair<int, VertexID>> sortedNeighbors;
        for (auto n : neighbors)
        {
            sortedNeighbors.push_back({graph.degree(n), n});
        }
        std::sort(sortedNeighbors.begin(), sortedNeighbors.end(),
                  [](const std::pair<int, VertexID> &a, const std::pair<int, VertexID> &b)
                  { return a.first > b.first; });

        int topCount = std::min(10, static_cast<int>(sortedNeighbors.size()));
        for (int i = 0; i < topCount; i++)
        {
            startVertices.push_back(sortedNeighbors[i].second);
        }

        if (neighbors.size() > topCount)
        {
            for (int i = 0; i < 5 && i + topCount < neighbors.size(); i++)
            {
                int randPos = topCount + (i * 7) % (neighbors.size() - topCount);
                startVertices.push_back(sortedNeighbors[randPos].second);
            }
        }
    }
    else
    {
        startVertices = neighbors;
    }

    for (auto start : startVertices)
    {
        Clique currentClique;
        currentClique.push_back(center);
        currentClique.push_back(start);

        bool extended;
        do
        {
            extended = false;
            VertexID bestVertex = graph.getNumVertices();
            int maxConnections = -1;

            for (auto n : neighbors)
            {

                if (std::find(currentClique.begin(), currentClique.end(), n) != currentClique.end())
                {
                    continue;
                }

                int connections = 0;
                bool canAdd = true;
                for (auto u : currentClique)
                {
                    if (graph.isAdjacent(n, u))
                    {
                        connections++;
                    }
                    else
                    {
                        canAdd = false;
                        break;
                    }
                }

                if (canAdd && connections > maxConnections)
                {
                    maxConnections = connections;
                    bestVertex = n;
                }
            }

            if (bestVertex < graph.getNumVertices())
            {
                currentClique.push_back(bestVertex);
                extended = true;
            }
        } while (extended);

        if (currentClique.size() > bestClique.size())
        {
            bestClique = currentClique;
        }
    }

    return bestClique;
}

int HelperAlgorithms::colorByDegeneracy(const Graph &graph, const std::vector<VertexID> &vertices)
{
    int n = vertices.size();
    if (n == 0)
        return 0;

    std::vector<VertexID> degeneracyOrder;
    std::vector<int> degrees(n);
    std::vector<bool> removed(n, false);

    for (int i = 0; i < n; i++)
    {
        VertexID v = vertices[i];
        int deg = 0;
        for (VertexID u : vertices)
        {
            if (v != u && graph.isAdjacent(v, u))
            {
                deg++;
            }
        }
        degrees[i] = deg;
    }

    for (int i = 0; i < n; i++)
    {

        int minDegIdx = -1;
        int minDeg = n;

        for (int j = 0; j < n; j++)
        {
            if (!removed[j] && degrees[j] < minDeg)
            {
                minDeg = degrees[j];
                minDegIdx = j;
            }
        }

        degeneracyOrder.push_back(vertices[minDegIdx]);

        removed[minDegIdx] = true;

        for (int j = 0; j < n; j++)
        {
            if (removed[j])
                continue;

            VertexID v = vertices[minDegIdx];
            VertexID u = vertices[j];

            if (graph.isAdjacent(v, u))
            {
                degrees[j]--;
            }
        }
    }

    std::vector<int> color(graph.getNumVertices(), -1);
    int maxColor = -1;

    for (VertexID v : degeneracyOrder)
    {
        std::vector<bool> used(n, false);

        for (VertexID u : graph.getNeighbors(v))
        {
            if (color[u] != -1)
            {
                used[color[u]] = true;
            }
        }

        int c = 0;
        while (c < n && used[c])
            c++;

        color[v] = c;
        maxColor = std::max(maxColor, c);
    }

    return maxColor + 1;
}

Clique HelperAlgorithms::exploreEgoNetwork(const Graph &graph, VertexID center)
{

    const auto &neighbors = graph.getNeighbors(center);

    if (neighbors.size() == 0)
    {
        Clique singletonClique;
        singletonClique.push_back(center);
        return singletonClique;
    }

    std::vector<VertexID> candidates(neighbors.begin(), neighbors.end());

    std::vector<int> localDegree(graph.getNumVertices(), 0);

    for (auto v : candidates)
    {
        for (auto u : candidates)
        {
            if (v != u && graph.isAdjacent(v, u))
            {
                localDegree[v]++;
            }
        }
    }

    bool reduced;
    int k = 1;

    do
    {
        reduced = false;
        std::vector<VertexID> newCandidates;

        for (auto v : candidates)
        {
            if (localDegree[v] >= k)
            {
                newCandidates.push_back(v);
            }
            else
            {

                reduced = true;

                for (auto u : candidates)
                {
                    if (u != v && graph.isAdjacent(u, v))
                    {
                        localDegree[u]--;
                    }
                }
            }
        }

        candidates = newCandidates;

        if (!reduced && !candidates.empty())
        {

            k++;
            reduced = true;
        }
    } while (reduced);

    Clique maximalClique;
    std::vector<bool> inClique(graph.getNumVertices(), false);

    std::vector<VertexID> degeneracyOrder;
    std::vector<bool> removed(graph.getNumVertices(), false);

    while (!candidates.empty())
    {

        VertexID minDegVertex = candidates[0];
        for (auto v : candidates)
        {
            if (localDegree[v] < localDegree[minDegVertex])
            {
                minDegVertex = v;
            }
        }

        degeneracyOrder.push_back(minDegVertex);

        candidates.erase(std::remove(candidates.begin(), candidates.end(), minDegVertex), candidates.end());
        removed[minDegVertex] = true;

        for (auto v : candidates)
        {
            if (graph.isAdjacent(v, minDegVertex))
            {
                localDegree[v]--;
            }
        }
    }

    for (auto it = degeneracyOrder.rbegin(); it != degeneracyOrder.rend(); ++it)
    {
        VertexID v = *it;

        bool canAdd = true;
        for (auto u : maximalClique)
        {
            if (!graph.isAdjacent(v, u))
            {
                canAdd = false;
                break;
            }
        }

        if (canAdd)
        {
            maximalClique.push_back(v);
        }
    }

    maximalClique.push_back(center);

    return maximalClique;
}

int HelperAlgorithms::computeUpperBound(const Graph &graph, const std::vector<VertexID> &candidates)
{

    return colorByDegeneracy(graph, candidates);
}