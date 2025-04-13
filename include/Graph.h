#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <vector>
#include <cstdint>

// Define types to match the original implementation's efficiency
typedef unsigned int VertexID;
typedef unsigned long EdgePtr;

class Graph
{
public:
	Graph();
	~Graph();

	// Add to public section of Graph.h:
	const std::vector<VertexID> &getDegeneracyOrder() const { return degeneracyOrder; }
	// Load graph from files (either binary or text format)
	bool loadFromFiles(const std::string &directory);

	// Load specifically from binary format
	bool loadBinary(const std::string &directory);

	// Load from text format as fallback
	bool loadText(const std::string &directory);

	// Preprocess the graph (compute degeneracy ordering, etc.)
	void preprocess();

	// Print summary statistics
	void printSummary() const;

	// Graph properties
	VertexID getNumVertices() const { return numVertices; }
	EdgePtr getNumEdges() const { return numEdges; }

	// Access vertex neighbors (for algorithms)
	const std::vector<VertexID> &getNeighbors(VertexID v) const { return adjList[v]; }

	// Check if two vertices are adjacent
	bool isAdjacent(VertexID u, VertexID v) const;

	// Get vertex degree
	VertexID degree(VertexID v) const { return adjList[v].size(); }

	// Data needed by algorithms
	std::vector<std::vector<VertexID>> adjList; // Adjacency list

private:
	// Graph properties
	VertexID numVertices; // Number of vertices
	EdgePtr numEdges;	  // Number of edges (undirected)

	// Computed properties for algorithms
	std::vector<VertexID> degeneracyOrder; // Vertices in degeneracy order
	std::vector<VertexID> coreness;		   // Coreness value for each vertex
};

#endif // GRAPH_H