# MC-BRB: Maximum Clique Solver

MC-BRB is a high-performance implementation of a branch-and-reduce-bound algorithm for finding the maximum clique in undirected graphs. This project provides efficient solutions for the Maximum Clique Problem, which has applications in social network analysis, computational biology, pattern recognition, and various optimization problems.

## Building the Project

### Prerequisites

- C++11 compatible compiler (GCC, Clang, or MSVC)
- Make

### Compilation

1. Clone the repository and navigate to the project directory
2. Run `make` to compile the project

```bash
git clone https://github.com/Schoksi20/Database
cd MC-BRB
make
```

This will create the executable `MC-BRB` in the current directory.

## Usage

The basic syntax for running the program is:

```bash
./MC-BRB MC-BRB .\datasets\Amazon\
```
Here you can replace the .\datasets\{your dataset}

## Input Format

The program supports two input formats:

### Binary Format

This is the preferred format for large graphs. The directory should contain:
- `b_degree.bin`: Contains vertex count, edge count, and degree of each vertex
- `b_adj.bin`: Contains adjacency lists for all vertices

### Text Format

If binary files are not found, the program will look for a text file:
- `graph.txt`: The first line contains the number of vertices and edges
  Following lines contain the edges (one per line, space-separated vertex IDs)

Example:
```
5 7
0 1
0 2
1 2
1 3
2 3
2 4
3 4
```

## Output

The program outputs:
- Maximum clique size
- Execution time
- The vertices in the maximum clique (if the clique size is ≤ 100)

## Performance Notes

- **Large Graphs**: For very large graphs (millions of vertices), the MC-DD or MC-EGO algorithms may be more practical than MC-BRB
- **Dense Graphs**: The algorithm automatically uses specialized data structures for dense graphs
- **Memory Usage**: Memory consumption scales with the size of the graph and the maximum degree of vertices

## Algorithm Details

- **MC-DD**: Degeneracy-based maximal clique finder. Fast but may not find the optimal solution.
- **MC-EGO**: Explores ego networks of high-degree vertices to find large cliques. Often provides good solutions quickly.
- **MC-BRB**: Exact algorithm that is guaranteed to find a maximum clique. Uses a combination of branch-and-bound with advanced reduction techniques.
- **verify**: Verifies if the clique in `clique.txt` is valid.

## Example Output

```
**** MC-BRB (Release) build at 09:45:32 Apr 10 2023 ****
**** Alg: MC-BRB, Graph: /data/socfb-Penn94 ****
Loaded graph from binary files.
Graph Summary:
  Vertices: 41536
  Edges: 1362220 (undirected)
Preprocessing the graph...
Preprocessing completed. Maximum coreness: 24
[Helper] Running MC-DD algorithm...
[Helper] MC-DD returns a clique of size 19
[Helper] Running Improved MC-EGO algorithm...
[Helper] MC-EGO found improved clique of size 21
[Helper] Improved MC-EGO returns a clique of size 21
[Solver] Initial best clique size: 21
[Solver] Starting branch-and-bound search...
[Solver] Found a new best clique of size 22
[Solver] Branch-and-bound search completed in 5823 ms.
[Solver] Search statistics:
  Best clique size: 22
  Branches explored: 1048576
  Maximum depth reached: 31
[Solver] Clique verification: PASSED
Final best clique size: 22
Total execution time: 9726815 microseconds
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
