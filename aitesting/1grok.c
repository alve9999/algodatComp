#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

// Edge structure to represent an edge between nodes u and v
typedef struct {
    size_t u;  // Node in set U (must be odd)
    size_t v;  // Node in set V (must be even)
} xedge_t;

// Structure to hold matching data and graph representation
typedef struct {
    size_t *pairU;      // Matching for nodes in U (odd nodes)
    size_t *pairV;      // Matching for nodes in V (even nodes)
    size_t *dist;       // Distance labels for BFS
    size_t *visited;    // Visited flags for DFS
    size_t *start;      // Start indices for adjacency list in CSR format
    size_t *adj;        // Adjacency list (neighbors in V for each U node)
    size_t n;           // Total number of nodes
} matching_data_t;

// Depth-First Search (DFS) to find an augmenting path
static bool DFS(size_t u, matching_data_t *data) {
    for (size_t i = data->start[u]; i < data->start[u + 1]; i++) {
        size_t v = data->adj[i];
        if (data->dist[v] == data->dist[u] + 1 && !data->visited[v]) {
            data->visited[v] = true;
            if (data->pairV[v] == 0 || DFS(data->pairV[v], data)) {
                data->pairU[u] = v;
                data->pairV[v] = u;
                return true;
            }
        }
    }
    return false;
}

// Public function to compute the maximum cardinality matching
size_t matching(size_t n, size_t m, xedge_t e[]) {
    // Allocate array to count degrees for odd nodes (U)
    size_t *deg = calloc(n + 1, sizeof(size_t));
    for (size_t i = 0; i < m; i++) {
        size_t u = e[i].u;
        size_t v = e[i].v;
        if (u % 2 == 1 && v % 2 == 0) { // Ensure u is odd, v is even
            deg[u]++;
        }
        // Invalid edges are silently ignored
    }

    // Build the start array for Compressed Sparse Row (CSR) format
    size_t *start = malloc((n + 2) * sizeof(size_t));
    start[1] = 0;
    for (size_t u = 1; u <= n; u++) {
        start[u + 1] = start[u] + deg[u];
    }

    // Build the adjacency list for odd nodes
    size_t *adj = malloc(m * sizeof(size_t));
    size_t *pos = malloc((n + 1) * sizeof(size_t));
    for (size_t u = 1; u <= n; u++) {
        if (u % 2 == 1) { // Only odd nodes have neighbors
            pos[u] = start[u];
        }
    }
    for (size_t i = 0; i < m; i++) {
        size_t u = e[i].u;
        size_t v = e[i].v;
        if (u % 2 == 1 && v % 2 == 0) {
            adj[pos[u]++] = v;
        }
    }

    // Allocate arrays for the matching algorithm
    size_t *pairU = calloc(n + 1, sizeof(size_t));  // Matches for U nodes
    size_t *pairV = calloc(n + 1, sizeof(size_t));  // Matches for V nodes
    size_t *dist = malloc((n + 1) * sizeof(size_t)); // Distance labels
    size_t *visited = malloc((n + 1) * sizeof(size_t)); // Visited flags
    size_t *queue = malloc((n + 1) * sizeof(size_t));   // Queue for BFS

    // Initialize matching data structure
    matching_data_t data = {pairU, pairV, dist, visited, start, adj, n};

    size_t max_matching = 0;
    while (true) {
        // BFS Phase: Build the layered graph
        for (size_t i = 1; i <= n; i++) {
            dist[i] = (size_t)-1; // Sentinel value for unvisited nodes
        }
        size_t head = 0, tail = 0;
        for (size_t u = 1; u <= n; u++) {
            if (u % 2 == 1 && pairU[u] == 0) { // Free odd nodes
                dist[u] = 0;
                queue[tail++] = u;
            }
        }
        bool can_reach_free_v = false;
        while (head < tail) {
            size_t current = queue[head++];
            if (current % 2 == 1) { // Node in U (odd)
                for (size_t i = start[current]; i < start[current + 1]; i++) {
                    size_t v = adj[i];
                    if (pairU[current] != v && dist[v] == (size_t)-1) {
                        dist[v] = dist[current] + 1;
                        queue[tail++] = v;
                        if (pairV[v] == 0) can_reach_free_v = true;
                    }
                }
            } else { // Node in V (even)
                size_t v = current;
                if (pairV[v] != 0 && dist[pairV[v]] == (size_t)-1) {
                    dist[pairV[v]] = dist[v] + 1;
                    queue[tail++] = pairV[v];
                }
            }
        }

        // If no free vertex in V is reachable, the maximum matching is found
        if (!can_reach_free_v) break;

        // DFS Phase: Find augmenting paths
        for (size_t i = 1; i <= n; i++) visited[i] = 0;
        for (size_t u = 1; u <= n; u++) {
            if (u % 2 == 1 && pairU[u] == 0 && DFS(u, &data)) {
                max_matching++;
            }
        }
    }

    // Clean up allocated memory
    free(deg);
    free(start);
    free(adj);
    free(pos);
    free(pairU);
    free(pairV);
    free(dist);
    free(visited);
    free(queue);

    return max_matching;
}
