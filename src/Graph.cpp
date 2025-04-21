#include "Graph.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <cstring>

Graph::Graph() : numVertices(0), numEdges(0) {}

Graph::~Graph() {}

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

bool Graph::loadBinary(const std::string &directory)
{

    std::string degreePath = directory + "/b_degree.bin";
    std::ifstream degreeFile(degreePath, std::ios::binary);
    if (!degreeFile.is_open())
    {
        return false;
    }

    int sizeCheck;
    degreeFile.read(reinterpret_cast<char *>(&sizeCheck), sizeof(int));
    if (sizeCheck != sizeof(int))
    {
        std::cerr << "Binary file format mismatch: int size differs." << std::endl;
        return false;
    }

    degreeFile.read(reinterpret_cast<char *>(&numVertices), sizeof(VertexID));
    degreeFile.read(reinterpret_cast<char *>(&numEdges), sizeof(EdgePtr));

    std::vector<VertexID> degrees(numVertices);
    degreeFile.read(reinterpret_cast<char *>(degrees.data()), numVertices * sizeof(VertexID));
    degreeFile.close();

    std::string adjPath = directory + "/b_adj.bin";
    std::ifstream adjFile(adjPath, std::ios::binary);
    if (!adjFile.is_open())
    {
        std::cerr << "Failed to open adjacency file: " << adjPath << std::endl;
        return false;
    }

    adjList.resize(numVertices);

    for (VertexID v = 0; v < numVertices; v++)
    {

        VertexID degree = degrees[v];
        adjList[v].resize(degree);

        if (degree > 0)
        {
            adjFile.read(reinterpret_cast<char *>(adjList[v].data()), degree * sizeof(VertexID));

            std::sort(adjList[v].begin(), adjList[v].end());
            auto newEnd = std::unique(adjList[v].begin(), adjList[v].end());
            adjList[v].erase(newEnd, adjList[v].end());

            adjList[v].erase(
                std::remove(adjList[v].begin(), adjList[v].end(), v),
                adjList[v].end());
        }
    }

    adjFile.close();

    EdgePtr totalDegree = 0;
    for (const auto &neighbors : adjList)
    {
        totalDegree += neighbors.size();
    }
    numEdges = totalDegree / 2;

    return true;
}

bool Graph::loadText(const std::string &directory)
{

    std::string filePath = directory + "/graph.txt";
    std::ifstream infile(filePath);
    if (!infile.is_open())
    {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        return false;
    }

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

    adjList.resize(numVertices);

    VertexID u, v;
    while (infile >> u >> v)
    {
        if (u >= numVertices || v >= numVertices)
        {
            std::cerr << "Invalid vertex IDs in edge: " << u << " - " << v << std::endl;
            continue;
        }

        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }

    for (VertexID v = 0; v < numVertices; v++)
    {
        std::sort(adjList[v].begin(), adjList[v].end());
        auto newEnd = std::unique(adjList[v].begin(), adjList[v].end());
        adjList[v].erase(newEnd, adjList[v].end());

        adjList[v].erase(
            std::remove(adjList[v].begin(), adjList[v].end(), v),
            adjList[v].end());
    }

    infile.close();
    return true;
}

bool Graph::isAdjacent(VertexID u, VertexID v) const
{
    if (u >= numVertices || v >= numVertices)
        return false;

    return std::binary_search(adjList[u].begin(), adjList[u].end(), v);
}

void Graph::preprocess()
{
    std::cout << "Preprocessing the graph..." << std::endl;

    degeneracyOrder.resize(numVertices);
    coreness.resize(numVertices, 0);

    std::vector<VertexID> vertexDegrees(numVertices);
    std::vector<bool> removed(numVertices, false);

    for (VertexID v = 0; v < numVertices; v++)
    {
        vertexDegrees[v] = adjList[v].size();
    }

    for (VertexID i = 0; i < numVertices; i++)
    {

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

        degeneracyOrder[i] = minDegVertex;
        coreness[minDegVertex] = minDeg;

        removed[minDegVertex] = true;

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

void Graph::printSummary() const
{
    std::cout << "Graph Summary:" << std::endl;
    std::cout << "  Vertices: " << numVertices << std::endl;
    std::cout << "  Edges: " << numEdges << " (undirected)" << std::endl;
}