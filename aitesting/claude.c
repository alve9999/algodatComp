#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>

// Constants
#define MAX_VERTICES 10000000  // Adjust as needed
#define NIL 0
#define INF INT_MAX

// Graph representation
typedef struct {
    int* edges;
    int* last;
    int* next;
    int edge_count;
    int max_edges;
} Graph;

// Function prototypes
void init_graph(Graph* graph, int max_edges);
void add_edge(Graph* graph, int u, int v);
void free_graph(Graph* graph);
bool bfs(Graph* graph, int u_size, int v_size, int* pairU, int* pairV, int* dist);
bool dfs(Graph* graph, int u, int* pairU, int* pairV, int* dist);
int hopcroft_karp(Graph* graph, int u_size, int v_size);

// Initialize graph
void init_graph(Graph* graph, int max_edges) {
    graph->max_edges = max_edges;
    graph->edge_count = 0;
    
    // +1 to accommodate 1-indexed vertices
    graph->edges = (int*)malloc((max_edges + 1) * sizeof(int));
    graph->next = (int*)malloc((max_edges + 1) * sizeof(int));
    graph->last = (int*)calloc(MAX_VERTICES + 1, sizeof(int));
    
    if (!graph->edges || !graph->next || !graph->last) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
}

// Add edge from u to v
void add_edge(Graph* graph, int u, int v) {
    if (graph->edge_count >= graph->max_edges) {
        fprintf(stderr, "Graph edge capacity exceeded\n");
        exit(1);
    }
    
    graph->edge_count++;
    graph->edges[graph->edge_count] = v;
    graph->next[graph->edge_count] = graph->last[u];
    graph->last[u] = graph->edge_count;
}

// Free graph resources
void free_graph(Graph* graph) {
    free(graph->edges);
    free(graph->next);
    free(graph->last);
}

// BFS to find augmenting paths
bool bfs(Graph* graph, int u_size, int v_size, int* pairU, int* pairV, int* dist) {
    int *queue = calloc(MAX_VERTICES + 1, sizeof(int));
    int front = 0, rear = 0;
    
    // Set distance to NIL for all vertices in U
    for (int u = 1; u <= u_size; u++) {
        if (pairU[u] == NIL) {
            dist[u] = 0;
            queue[rear++] = u;
        } else {
            dist[u] = INF;
        }
    }
    
    dist[NIL] = INF;
    
    // BFS
    while (front < rear) {
        int u = queue[front++];
        
        if (dist[u] < dist[NIL]) {
            // Iterate through all adjacent vertices of u
            for (int e = graph->last[u]; e; e = graph->next[e]) {
                int v = graph->edges[e];
                
                // If v is not matched or the distance to the matched vertex is not set
                if (dist[pairV[v]] == INF) {
                    dist[pairV[v]] = dist[u] + 1;
                    queue[rear++] = pairV[v];
                }
            }
        }
    }
    
    // Return true if we found an augmenting path to NIL
    return dist[NIL] != INF;
}

// DFS to find augmenting paths
bool dfs(Graph* graph, int u, int* pairU, int* pairV, int* dist) {
    if (u == NIL) {
        return true;
    }
    
    // Try all adjacent vertices
    for (int e = graph->last[u]; e; e = graph->next[e]) {
        int v = graph->edges[e];
        
        // Follow the distances in the shortest path
        if (dist[pairV[v]] == dist[u] + 1) {
            if (dfs(graph, pairV[v], pairU, pairV, dist)) {
                pairV[v] = u;
                pairU[u] = v;
                return true;
            }
        }
    }
    
    // No augmenting path found through u
    dist[u] = INF;
    return false;
}

// Hopcroft-Karp algorithm
int hopcroft_karp(Graph* graph, int u_size, int v_size) {
    // Allocate memory for pair arrays and distance array
    int* pairU = (int*)calloc(u_size + 1, sizeof(int));
    int* pairV = (int*)calloc(v_size + 1, sizeof(int));
    int* dist = (int*)malloc((u_size + 1) * sizeof(int));
    
    if (!pairU || !pairV || !dist) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    
    int max_matching = 0;
    
    // Main algorithm loop
    while (bfs(graph, u_size, v_size, pairU, pairV, dist)) {
        // Try to find augmenting paths for all free vertices in U
        for (int u = 1; u <= u_size; u++) {
            if (pairU[u] == NIL && dfs(graph, u, pairU, pairV, dist)) {
                max_matching++;
            }
        }
    }
    
    // Free allocated memory
    free(pairU);
    free(pairV);
    free(dist);
    
    return max_matching;
}

int main() {
    // Buffer for reading lines
    char line[100];
    int max_edges = 1000000000;  // Initial estimate
    
    Graph graph;
    init_graph(&graph, max_edges);
    
    int max_u = 0, max_v = 0;
    int u, v;
    
    // Read edges from stdin
    while (fgets(line, sizeof(line), stdin) != NULL) {
        if (sscanf(line, "%d %d", &u, &v) == 2) {
            add_edge(&graph, u, v);
            
            // Track the max vertex indices for U and V
            if (u > max_u) max_u = u;
            if (v > max_v) max_v = v;
        }
    }
    
    // Find maximum matching
    int result = hopcroft_karp(&graph, max_u, max_v);
    
    // Output the result
    printf("%d\n", result);
    
    // Clean up
    free_graph(&graph);
    
    return 0;
}
