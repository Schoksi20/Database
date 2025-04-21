#include "Graph.h"
#include "MCBRBSolver.h"
#include <iostream>
#include <string>
#include <cstring>
#include <chrono>
#include <fstream>
void printUsage()
{
    std::cout << "Usage: MC-BRB [algorithm] [graph-dir] [output]" << std::endl;
    std::cout << "  algorithm: MC-DD, MC-EGO, MC-BRB, verify" << std::endl;
    std::cout << "  graph-dir: Directory containing graph files" << std::endl;
    std::cout << "  output: Optional. If 'output' is specified, write clique to file" << std::endl;
}

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        printUsage();
        return 1;
    }

#ifndef NDEBUG
    std::cout << "**** MC-BRB (Debug) build at " << __TIME__ << " " << __DATE__ << " ****" << std::endl;
#else
    std::cout << "**** MC-BRB (Release) build at " << __TIME__ << " " << __DATE__ << " ****" << std::endl;
#endif

    std::cout << "**** Alg: " << argv[1] << ", Graph: " << argv[2] << " ****" << std::endl;

    auto startTime = std::chrono::high_resolution_clock::now();

    Graph graph;
    if (!graph.loadFromFiles(argv[2]))
    {
        std::cerr << "Failed to load graph from " << argv[2] << std::endl;
        return 1;
    }

    graph.printSummary();
    graph.preprocess();

    std::string algorithm = argv[1];
    std::vector<VertexID> bestClique;

    if (algorithm == "MC-DD")
    {
        bestClique = HelperAlgorithms::runMCDD(graph);
    }
    else if (algorithm == "MC-EGO")
    {
        bestClique = HelperAlgorithms::runMCEGO(graph);
    }
    else if (algorithm == "MC-BRB")
    {

        MCBRBSolver solver(graph);
        solver.initialize();

        solver.branchAndBound();

        bestClique = solver.getBestClique();
    }
    else if (algorithm == "verify")
    {

        std::cout << "Verification not yet implemented." << std::endl;
        return 0;
    }
    else
    {
        printUsage();
        return 1;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);

    std::cout << "Final best clique size: " << bestClique.size() << std::endl;
    std::cout << "Total execution time: " << duration.count() << " microseconds" << std::endl;

    if (bestClique.size() <= 100)
    {
        std::cout << "Clique vertices: ";
        for (auto v : bestClique)
        {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }

    if (argc >= 4 && std::string(argv[3]) == "output")
    {
        std::cout << "Writing clique to file: clique.txt" << std::endl;
        std::ofstream outFile("clique.txt");
        if (outFile.is_open())
        {
            outFile << bestClique.size() << std::endl;
            for (auto v : bestClique)
            {
                outFile << v << std::endl;
            }
            outFile.close();
        }
        else
        {
            std::cerr << "Failed to open output file." << std::endl;
        }
    }

    return 0;
}