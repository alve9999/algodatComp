#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>

// Constants
#define NIL 0
#define INF INT_MAX

// Graph representation
typedef struct {
    int n;
    int* edges;
    int* last;
    int* next;
    int edge_count;
    int max_edges;
} Graph;

// Function prototypes
void add_edge(Graph* graph, int u, int v);
void free_graph(Graph* graph);
bool bfs(Graph* graph, int u_size, int v_size, int* pairU, int* pairV, int* dist);
bool dfs(Graph* graph, int u, int* pairU, int* pairV, int* dist);
int hopcroft_karp(Graph* graph, int u_size, int v_size);

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
    int *queue = calloc(graph->n + 1, sizeof(int));
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


typedef struct {
    size_t        u;    /* odd.        */
    size_t        v;    /* even.    */
} xedge_t;


size_t matching(size_t n, size_t m, xedge_t e[]) {
    Graph graph1;
    graph1.n = n;
    {
        Graph *graph = &graph1;
        graph->max_edges = m;
        graph->edge_count = 0;

        // +1 to accommodate 1-indexed vertices
        graph->edges = (int*)malloc((m + 1) * sizeof(int));
        graph->next = (int*)malloc((m + 1) * sizeof(int));
        graph->last = (int*)calloc(n + 1, sizeof(int));

        if (!graph->edges || !graph->next || !graph->last) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
    }
    for (int i = 0; i < m; ++i) {
        xedge_t *ed = &e[i];
        add_edge(&graph1, ed->u, ed->v);
    }
    int result = hopcroft_karp(&graph1, n, n);
    return result;
}
