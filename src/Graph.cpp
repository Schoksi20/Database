#include "Graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cstring>

// Constructor: initialize members
Graph::Graph() : numVertices(0), numEdges(0) {}

// Destructor
Graph::~Graph() {}

// Main loading function that tries binary first, then text
bool Graph::loadFromFiles(const std::string &directory)
{
    if (loadBinary(directory))
    {
        std::cout << "Loaded graph from binary files." << std::endl;
        return true;
    }
    else
    {
        std::cout << "Binary files not found, trying text format..." << std::endl;
        return loadText(directory);
    }
}

// Load binary format (similar to original implementation)
bool Graph::loadBinary(const std::string &directory)
{
    // Try to open the binary degree file
    std::string degreePath = directory + "/b_degree.bin";
    std::ifstream degreeFile(degreePath, std::ios::binary);
    if (!degreeFile.is_open())
    {
        return false; // Binary file not found
    }

    // Read file header and validate format
    int sizeCheck;
    degreeFile.read(reinterpret_cast<char *>(&sizeCheck), sizeof(int));
    if (sizeCheck != sizeof(int))
    {
        std::cerr << "Binary file format mismatch: int size differs." << std::endl;
        return false;
    }

    // Read vertex and edge counts
    degreeFile.read(reinterpret_cast<char *>(&numVertices), sizeof(VertexID));
    degreeFile.read(reinterpret_cast<char *>(&numEdges), sizeof(EdgePtr));

    // Read vertex degrees
    std::vector<VertexID> degrees(numVertices);
    degreeFile.read(reinterpret_cast<char *>(degrees.data()), numVertices * sizeof(VertexID));
    degreeFile.close();

    // Open adjacency list file
    std::string adjPath = directory + "/b_adj.bin";
    std::ifstream adjFile(adjPath, std::ios::binary);
    if (!adjFile.is_open())
    {
        std::cerr << "Failed to open adjacency file: " << adjPath << std::endl;
        return false;
    }

    // Resize adjacency list
    adjList.resize(numVertices);

    // For each vertex, read its neighbors
    for (VertexID v = 0; v < numVertices; v++)
    {
        // Allocate space for neighbors
        VertexID degree = degrees[v];
        adjList[v].resize(degree);

        // Read neighbor IDs if degree > 0
        if (degree > 0)
        {
            adjFile.read(reinterpret_cast<char *>(adjList[v].data()), degree * sizeof(VertexID));

            // Remove self loops and parallel edges
            std::sort(adjList[v].begin(), adjList[v].end());
            auto newEnd = std::unique(adjList[v].begin(), adjList[v].end());
            adjList[v].erase(newEnd, adjList[v].end());

            // Remove self loops (vertex v itself)
            adjList[v].erase(
                std::remove(adjList[v].begin(), adjList[v].end(), v),
                adjList[v].end());
        }
    }

    adjFile.close();

    // Recalculate edge count (undirected)
    EdgePtr totalDegree = 0;
    for (const auto &neighbors : adjList)
    {
        totalDegree += neighbors.size();
    }
    numEdges = totalDegree / 2; // Each edge counted twice in undirected graph

    return true;
}

// Load from text format as a fallback
bool Graph::loadText(const std::string &directory)
{
    // Try to open a text format file
    std::string filePath = directory + "/graph.txt";
    std::ifstream infile(filePath);
    if (!infile.is_open())
    {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return false;
    }

    // Assume the first line contains number of vertices and edges
    std::string line;
    if (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (!(iss >> numVertices >> numEdges))
        {
            std::cerr << "Error reading graph summary from file." << std::endl;
            return false;
        }
    }

    // Resize the adjacency list
    adjList.resize(numVertices);

    // Read edges. Each subsequent line contains two integers: u and v
    VertexID u, v;
    while (infile >> u >> v)
    {
        if (u >= numVertices || v >= numVertices)
        {
            std::cerr << "Invalid vertex IDs in edge: " << u << " - " << v << std::endl;
            continue;
        }

        // Add edge (undirected)
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }

    // Sort adjacency lists and remove duplicates
    for (VertexID v = 0; v < numVertices; v++)
    {
        std::sort(adjList[v].begin(), adjList[v].end());
        auto newEnd = std::unique(adjList[v].begin(), adjList[v].end());
        adjList[v].erase(newEnd, adjList[v].end());

        // Remove self loops
        adjList[v].erase(
            std::remove(adjList[v].begin(), adjList[v].end(), v),
            adjList[v].end());
    }

    infile.close();
    return true;
}

// Check if two vertices are adjacent
bool Graph::isAdjacent(VertexID u, VertexID v) const
{
    if (u >= numVertices || v >= numVertices)
        return false;

    // Use binary search on sorted adjacency list
    return std::binary_search(adjList[u].begin(), adjList[u].end(), v);
}

// Preprocessing: This can include computing degeneracy ordering and coreness
void Graph::preprocess()
{
    std::cout << "Preprocessing the graph..." << std::endl;

    // Initialize degeneracy ordering and coreness
    degeneracyOrder.resize(numVertices);
    coreness.resize(numVertices, 0);

    // Compute degeneracy ordering and coreness
    std::vector<VertexID> vertexDegrees(numVertices);
    std::vector<bool> removed(numVertices, false);

    // Initialize degrees
    for (VertexID v = 0; v < numVertices; v++)
    {
        vertexDegrees[v] = adjList[v].size();
    }

    // Process vertices in non-decreasing order of degree
    for (VertexID i = 0; i < numVertices; i++)
    {
        // Find unremoved vertex with minimum degree
        VertexID minDegVertex = 0;
        VertexID minDeg = numVertices + 1;

        for (VertexID v = 0; v < numVertices; v++)
        {
            if (!removed[v] && vertexDegrees[v] < minDeg)
            {
                minDeg = vertexDegrees[v];
                minDegVertex = v;
            }
        }

        // Add to degeneracy ordering and record coreness
        degeneracyOrder[i] = minDegVertex;
        coreness[minDegVertex] = minDeg;

        // Mark as removed
        removed[minDegVertex] = true;

        // Update degrees of adjacent vertices
        for (VertexID neighbor : adjList[minDegVertex])
        {
            if (!removed[neighbor])
            {
                vertexDegrees[neighbor]--;
            }
        }
    }

    std::cout << "Preprocessing completed. Maximum coreness: "
              << *std::max_element(coreness.begin(), coreness.end()) << std::endl;
}

// Print a summary of the graph
void Graph::printSummary() const
{
    std::cout << "Graph Summary:" << std::endl;
    std::cout << "  Vertices: " << numVertices << std::endl;
    std::cout << "  Edges: " << numEdges << " (undirected)" << std::endl;
}