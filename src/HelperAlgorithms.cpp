#include "HelperAlgorithms.h"
#include <iostream>
#include <algorithm>
#include <queue>
#include <unordered_set>

// MC-DD: Degeneracy-based maximum clique algorithm
Clique HelperAlgorithms::runMCDD(const Graph &graph)
{
    std::cout << "[Helper] Running MC-DD algorithm..." << std::endl;

    // Implementation of the degeneracy-based approach
    Clique maximalClique = findMaximalCliqueByDegeneracy(graph);

    std::cout << "[Helper] MC-DD returns a clique of size " << maximalClique.size() << std::endl;
    return maximalClique;
}

// MC-EGO: Ego-centric maximum clique algorithm
Clique HelperAlgorithms::runMCEGO(const Graph &graph)
{
    std::cout << "[Helper] Running MC-EGO algorithm..." << std::endl;

    // Initialize best clique
    Clique bestClique;

    // Try each vertex as the center of an ego network
    for (VertexID v = 0; v < graph.getNumVertices(); v++)
    {
        // Only process vertices with enough neighbors
        if (graph.degree(v) < bestClique.size())
        {
            continue;
        }

        // Find maximal clique in ego network of v
        Clique egoClique = exploreEgoNetwork(graph, v);

        // Update best clique if better
        if (egoClique.size() > bestClique.size())
        {
            bestClique = egoClique;
        }
    }

    std::cout << "[Helper] MC-EGO returns a clique of size " << bestClique.size() << std::endl;
    return bestClique;
}

// Find a maximal clique by degeneracy ordering
Clique HelperAlgorithms::findMaximalCliqueByDegeneracy(const Graph &graph)
{
    // Compute degeneracy ordering
    std::vector<VertexID> degeneracyOrder;
    std::vector<int> degrees(graph.getNumVertices());
    std::vector<bool> removed(graph.getNumVertices(), false);

    // Initialize degrees
    for (VertexID v = 0; v < graph.getNumVertices(); v++)
    {
        degrees[v] = graph.degree(v);
    }

    // Process vertices in non-decreasing order of degree
    for (VertexID i = 0; i < graph.getNumVertices(); i++)
    {
        // Find unremoved vertex with minimum degree
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

        // Add to degeneracy ordering
        degeneracyOrder.push_back(minDegVertex);

        // Mark as removed
        removed[minDegVertex] = true;

        // Update degrees of adjacent vertices
        for (VertexID neighbor : graph.getNeighbors(minDegVertex))
        {
            if (!removed[neighbor])
            {
                degrees[neighbor]--;
            }
        }
    }

    // Find a maximal clique by adding vertices in reverse degeneracy order
    Clique maximalClique;
    std::vector<bool> inClique(graph.getNumVertices(), false);

    for (auto it = degeneracyOrder.rbegin(); it != degeneracyOrder.rend(); ++it)
    {
        VertexID v = *it;

        // Check if v is adjacent to all vertices in the current clique
        bool canAdd = true;
        for (VertexID u : maximalClique)
        {
            if (!graph.isAdjacent(v, u))
            {
                canAdd = false;
                break;
            }
        }

        // Add v to clique if possible
        if (canAdd)
        {
            maximalClique.push_back(v);
            inClique[v] = true;
        }
    }

    return maximalClique;
}

// Enhanced MC-EGO implementation
Clique HelperAlgorithms::runImprovedMCEGO(const Graph &graph)
{
    std::cout << "[Helper] Running Improved MC-EGO algorithm..." << std::endl;

    // Initialize best clique
    Clique bestClique;

    // Get the degeneracy ordering - we'll explore vertices in this order for better results
    std::vector<VertexID> degOrder;
    std::vector<int> degrees(graph.getNumVertices());
    std::vector<bool> removed(graph.getNumVertices(), false);

    // Initialize degrees
    for (VertexID v = 0; v < graph.getNumVertices(); v++)
    {
        degrees[v] = graph.degree(v);
    }

    // Create degeneracy ordering
    for (VertexID i = 0; i < graph.getNumVertices(); i++)
    {
        // Find unremoved vertex with minimum degree
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

        // Add to degeneracy ordering
        degOrder.push_back(minDegVertex);

        // Mark as removed
        removed[minDegVertex] = true;

        // Update degrees of adjacent vertices
        for (VertexID neighbor : graph.getNeighbors(minDegVertex))
        {
            if (!removed[neighbor])
            {
                degrees[neighbor]--;
            }
        }
    }

    // Reverse the order to put high-degree vertices first
    std::reverse(degOrder.begin(), degOrder.end());

    // Process vertices in the ordering
    for (VertexID v : degOrder)
    {
        // Only process vertices with enough neighbors
        if (graph.degree(v) < bestClique.size())
        {
            continue;
        }

        // Find maximal clique in ego network of v
        Clique egoClique = improvedExploreEgoNetwork(graph, v);

        // Update best clique if better
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

// Add this function to HelperAlgorithms.cpp to improve ego network processing
Clique HelperAlgorithms::improvedExploreEgoNetwork(const Graph &graph, VertexID center)
{
    // Get neighbors of center
    const auto &neighbors = graph.getNeighbors(center);

    // If not enough neighbors, return singleton clique
    if (neighbors.empty())
    {
        Clique singletonClique;
        singletonClique.push_back(center);
        return singletonClique;
    }

    // Create candidate set (neighbors of center)
    std::vector<VertexID> candidates(neighbors.begin(), neighbors.end());

    // Compute density of the neighborhood
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

    // For dense neighborhoods, use a more sophisticated approach
    if (density > 0.5 && candidates.size() > 10)
    {
        return findMaximalCliqueInDenseNeighborhood(graph, center, candidates);
    }

    // Filter candidates to find max(k)-core in the neighborhood
    std::vector<int> localDegree(graph.getNumVertices(), 0);

    // Compute local degrees
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

    // Improved k-core reduction - determine max k dynamically
    int maxK = 0;
    for (auto v : candidates)
    {
        maxK = std::max(maxK, localDegree[v]);
    }

    // Start with a higher k-value to quickly reduce the subgraph
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
                // Update degrees of neighbors
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
            // If all vertices removed, try a lower k
            k = k / 2;
            reduced = true;
            newCandidates = candidates;
        }
        else
        {
            candidates = newCandidates;

            // If stable at this k, try increasing it
            if (!reduced && !candidates.empty() && k < maxK)
            {
                k++;
                reduced = true;
            }
        }
    }

    // Use degeneracy ordering for better clique finding
    std::vector<VertexID> degeneracyOrder;
    std::vector<bool> vertexRemoved(graph.getNumVertices(), false);

    // Create degeneracy ordering
    while (!candidates.empty())
    {
        // Find vertex with minimum local degree
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

        // Add to ordering
        degeneracyOrder.push_back(minDegVertex);
        vertexRemoved[minDegVertex] = true;

        // Remove from candidates
        candidates[minIdx] = candidates.back();
        candidates.pop_back();

        // Update local degrees
        for (auto v : candidates)
        {
            if (graph.isAdjacent(v, minDegVertex))
            {
                localDegree[v]--;
            }
        }
    }

    // Find maximal clique in reverse degeneracy order
    Clique maximalClique;

    // Process vertices in reverse degeneracy order
    for (auto it = degeneracyOrder.rbegin(); it != degeneracyOrder.rend(); ++it)
    {
        VertexID v = *it;

        // Check if v is adjacent to all vertices in the current clique
        bool canAdd = true;
        for (auto u : maximalClique)
        {
            if (!graph.isAdjacent(v, u))
            {
                canAdd = false;
                break;
            }
        }

        // Add v to clique if possible
        if (canAdd)
        {
            maximalClique.push_back(v);
        }
    }

    // Add the center vertex
    maximalClique.push_back(center);

    // Try to extend the clique greedily
    bool extended;
    do
    {
        extended = false;
        for (auto n : neighbors)
        {
            // Skip vertices already in the clique
            if (std::find(maximalClique.begin(), maximalClique.end(), n) != maximalClique.end())
            {
                continue;
            }

            // Check if n is adjacent to all vertices in the clique
            bool canAdd = true;
            for (auto u : maximalClique)
            {
                if (!graph.isAdjacent(n, u))
                {
                    canAdd = false;
                    break;
                }
            }

            // Add n to clique if possible
            if (canAdd)
            {
                maximalClique.push_back(n);
                extended = true;
                break; // Start over with the new clique
            }
        }
    } while (extended);

    return maximalClique;
}

// Helper for dense neighborhoods
Clique HelperAlgorithms::findMaximalCliqueInDenseNeighborhood(const Graph &graph, VertexID center, const std::vector<VertexID> &neighbors)
{
    // For dense neighborhoods, try several starting points
    Clique bestClique;
    bestClique.push_back(center); // Start with just the center

    // Try multiple starting vertices for diversity
    std::vector<VertexID> startVertices;

    // Select a diverse set of high-degree vertices
    if (neighbors.size() > 20)
    {
        // Sort neighbors by degree
        std::vector<std::pair<int, VertexID>> sortedNeighbors;
        for (auto n : neighbors)
        {
            sortedNeighbors.push_back({graph.degree(n), n});
        }
        std::sort(sortedNeighbors.begin(), sortedNeighbors.end(),
                  [](const std::pair<int, VertexID> &a, const std::pair<int, VertexID> &b)
                  { return a.first > b.first; });

        // Take top 10 and some random ones for diversity
        int topCount = std::min(10, static_cast<int>(sortedNeighbors.size()));
        for (int i = 0; i < topCount; i++)
        {
            startVertices.push_back(sortedNeighbors[i].second);
        }

        // Add some random vertices for diversity
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

    // Try each starting vertex
    for (auto start : startVertices)
    {
        Clique currentClique;
        currentClique.push_back(center);
        currentClique.push_back(start);

        // Add vertices greedily
        bool extended;
        do
        {
            extended = false;
            VertexID bestVertex = graph.getNumVertices();
            int maxConnections = -1;

            for (auto n : neighbors)
            {
                // Skip vertices already in the clique
                if (std::find(currentClique.begin(), currentClique.end(), n) != currentClique.end())
                {
                    continue;
                }

                // Count connections to the current clique
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

                // If fully connected and has more overall connections, prefer this vertex
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

        // Keep the best clique found
        if (currentClique.size() > bestClique.size())
        {
            bestClique = currentClique;
        }
    }

    return bestClique;
}

// Color vertices by degeneracy ordering
int HelperAlgorithms::colorByDegeneracy(const Graph &graph, const std::vector<VertexID> &vertices)
{
    int n = vertices.size();
    if (n == 0)
        return 0;

    // Compute degeneracy ordering
    std::vector<VertexID> degeneracyOrder;
    std::vector<int> degrees(n);
    std::vector<bool> removed(n, false);

    // Initialize degrees
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

    // Process vertices in non-decreasing order of degree
    for (int i = 0; i < n; i++)
    {
        // Find unremoved vertex with minimum degree
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

        // Add to degeneracy ordering
        degeneracyOrder.push_back(vertices[minDegIdx]);

        // Mark as removed
        removed[minDegIdx] = true;

        // Update degrees of adjacent vertices
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

    // Color vertices in the computed order
    std::vector<int> color(graph.getNumVertices(), -1);
    int maxColor = -1;

    for (VertexID v : degeneracyOrder)
    {
        std::vector<bool> used(n, false);

        // Mark colors used by neighbors
        for (VertexID u : graph.getNeighbors(v))
        {
            if (color[u] != -1)
            {
                used[color[u]] = true;
            }
        }

        // Find smallest available color
        int c = 0;
        while (c < n && used[c])
            c++;

        color[v] = c;
        maxColor = std::max(maxColor, c);
    }

    return maxColor + 1; // Number of colors used
}

// Explore ego network of a vertex to find a maximal clique
Clique HelperAlgorithms::exploreEgoNetwork(const Graph &graph, VertexID center)
{
    // Get neighbors of center
    const auto &neighbors = graph.getNeighbors(center);

    // If not enough neighbors, return empty clique
    if (neighbors.size() == 0)
    {
        Clique singletonClique;
        singletonClique.push_back(center);
        return singletonClique;
    }

    // Create candidate set (neighbors of center)
    std::vector<VertexID> candidates(neighbors.begin(), neighbors.end());

    // Filter candidates to find max(k)-core in the neighborhood
    std::vector<int> localDegree(graph.getNumVertices(), 0);

    // Compute local degrees
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

    // K-core reduction
    bool reduced;
    int k = 1; // Start with 1-core

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
                // Remove vertices with degree < k
                reduced = true;

                // Update degrees of neighbors
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
            // Try next k-value
            k++;
            reduced = true;
        }
    } while (reduced);

    // Use degeneracy-based approach to find maximal clique in the filtered neighborhood
    Clique maximalClique;
    std::vector<bool> inClique(graph.getNumVertices(), false);

    // Compute degeneracy ordering on candidates
    std::vector<VertexID> degeneracyOrder;
    std::vector<bool> removed(graph.getNumVertices(), false);

    while (!candidates.empty())
    {
        // Find vertex with minimum local degree
        VertexID minDegVertex = candidates[0];
        for (auto v : candidates)
        {
            if (localDegree[v] < localDegree[minDegVertex])
            {
                minDegVertex = v;
            }
        }

        // Add to degeneracy ordering
        degeneracyOrder.push_back(minDegVertex);

        // Remove from candidates
        candidates.erase(std::remove(candidates.begin(), candidates.end(), minDegVertex), candidates.end());
        removed[minDegVertex] = true;

        // Update local degrees
        for (auto v : candidates)
        {
            if (graph.isAdjacent(v, minDegVertex))
            {
                localDegree[v]--;
            }
        }
    }

    // Find maximal clique in reverse degeneracy order
    for (auto it = degeneracyOrder.rbegin(); it != degeneracyOrder.rend(); ++it)
    {
        VertexID v = *it;

        // Check if v is adjacent to all vertices in current clique
        bool canAdd = true;
        for (auto u : maximalClique)
        {
            if (!graph.isAdjacent(v, u))
            {
                canAdd = false;
                break;
            }
        }

        // Add to clique if possible
        if (canAdd)
        {
            maximalClique.push_back(v);
        }
    }

    // Add the center vertex to the clique
    maximalClique.push_back(center);

    return maximalClique;
}

// Compute an upper bound for the maximum clique size
int HelperAlgorithms::computeUpperBound(const Graph &graph, const std::vector<VertexID> &candidates)
{
    // Use coloring-based bound
    return colorByDegeneracy(graph, candidates);
}