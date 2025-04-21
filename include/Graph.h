#ifndef GRAPH_H
#define GRAPH_H

#include <string>
#include <vector>
#include <cstdint>

typedef unsigned int VertexID;
typedef unsigned long EdgePtr;

class Graph
{
public:
	Graph();
	~Graph();

	const std::vector<VertexID> &getDegeneracyOrder() const { return degeneracyOrder; }

	bool loadFromFiles(const std::string &directory);

	bool loadBinary(const std::string &directory);

	bool loadText(const std::string &directory);

	void preprocess();

	void printSummary() const;

	VertexID getNumVertices() const { return numVertices; }
	EdgePtr getNumEdges() const { return numEdges; }

	const std::vector<VertexID> &getNeighbors(VertexID v) const { return adjList[v]; }

	bool isAdjacent(VertexID u, VertexID v) const;

	VertexID degree(VertexID v) const { return adjList[v].size(); }

	std::vector<std::vector<VertexID>> adjList;

private:
	VertexID numVertices;
	EdgePtr numEdges;

	std::vector<VertexID> degeneracyOrder;
	std::vector<VertexID> coreness;
};

#endif