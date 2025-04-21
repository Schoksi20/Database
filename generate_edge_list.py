import random
import argparse

def generate_edge_list(num_vertices, num_edges, output_file):
    """
    Generate a random undirected graph edge list.
    
    Parameters:
        num_vertices: Number of vertices in the graph
        num_edges: Approximate number of undirected edges
        output_file: File to save the edge list
    """
    print(f"Generating random graph with {num_vertices} vertices and ~{num_edges} edges...")
    
    # Check if requested edges exceeds maximum possible
    possible_edges = num_vertices * (num_vertices - 1) // 2
    if num_edges > possible_edges:
        print(f"Warning: Requested edges ({num_edges}) exceeds maximum possible ({possible_edges})")
        num_edges = possible_edges
    
    # Keep track of edges to avoid duplicates
    edges = set()
    
    with open(output_file, 'w') as f:
        # Write header
        f.write(f"# Undirected graph: {output_file}\n")
        f.write(f"# Vertices: {num_vertices} Edges: {num_edges}\n")
        
        # Generate random edges
        while len(edges) < num_edges:
            u = random.randint(0, num_vertices - 1)
            v = random.randint(0, num_vertices - 1)
            
            # Skip self-loops and ensure u < v to avoid duplicates
            if u != v:
                edge = (min(u, v), max(u, v))
                if edge not in edges:
                    edges.add(edge)
                    
                    if len(edges) % 100000 == 0:
                        print(f"Generated {len(edges)} edges...")
        
        # Write edges to file
        for u, v in edges:
            f.write(f"{u} {v}\n")
    
    print(f"Generated {len(edges)} edges and saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a random graph edge list')
    parser.add_argument('--vertices', type=int, default=1000, 
                        help='Number of vertices in the graph')
    parser.add_argument('--edges', type=int, default=5000, 
                        help='Number of edges in the graph')
    parser.add_argument('--output', type=str, default='edges.txt', 
                        help='Output edge list file')
    
    args = parser.parse_args()
    
    generate_edge_list(args.vertices, args.edges, args.output)