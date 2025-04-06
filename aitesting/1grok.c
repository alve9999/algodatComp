#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

// Edge structure
typedef struct {
    size_t u;
    size_t v;
} xedge_t;

// Structure for matching data
typedef struct {
    size_t *pairU;      // Matching for nodes in U
    size_t *pairV;      // Matching for nodes in V
    size_t *dist;       // Distance labels from BFS
    size_t *visited;    // Visited flags for DFS
    size_t *start;      // Start indices for adjacency list
    size_t *adj;        // Adjacency list (neighbors in V)
    size_t left_size;   // Number of nodes in U (n/2)
    size_t n;           // Total number of nodes
} matching_data_t;

// DFS to find an augmenting path
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

// Compute maximum matching
size_t matching(size_t n, size_t m, xedge_t e[]) {
    size_t left_size = n / 2;

    // Build adjacency list in CSR format
    size_t *deg = calloc(left_size + 1, sizeof(size_t));
    for (size_t i = 0; i < m; i++) {
        deg[e[i].u]++;
    }
    size_t *start = malloc((left_size + 2) * sizeof(size_t));
    start[1] = 0;
    for (size_t u = 1; u <= left_size; u++) {
        start[u + 1] = start[u] + deg[u];
    }
    size_t *adj = malloc(m * sizeof(size_t));
    size_t *pos = malloc((left_size + 1) * sizeof(size_t));
    for (size_t u = 1; u <= left_size; u++) {
        pos[u] = start[u];
    }
    for (size_t i = 0; i < m; i++) {
        size_t u = e[i].u;
        size_t v = e[i].v;
        adj[pos[u]++] = v;
    }

    // Allocate arrays for matching
    size_t *pairU = calloc(left_size + 1, sizeof(size_t));
    size_t *pairV = calloc(n + 1, sizeof(size_t));
    size_t *dist = malloc((n + 1) * sizeof(size_t));
    size_t *visited = malloc((n + 1) * sizeof(size_t));
    size_t *queue = malloc((n + 1) * sizeof(size_t));

    // Initialize matching data
    matching_data_t data = {pairU, pairV, dist, visited, start, adj, left_size, n};

    size_t max_matching = 0;
    while (true) {
        // BFS Phase: Build layered graph
        for (size_t i = 1; i <= n; i++) {
            dist[i] = (size_t)-1; // Use (size_t)-1 as sentinel
        }
        size_t head = 0, tail = 0;
        for (size_t u = 1; u <= left_size; u++) {
            if (pairU[u] == 0) {
                dist[u] = 0;
                queue[tail++] = u;
            }
        }
        bool can_reach_free_v = false;
        while (head < tail) {
            size_t current = queue[head++];
            if (current <= left_size) { // Node in U
                for (size_t i = start[current]; i < start[current + 1]; i++) {
                    size_t v = adj[i];
                    if (pairU[current] != v && dist[v] == (size_t)-1) {
                        dist[v] = dist[current] + 1;
                        queue[tail++] = v;
                        if (pairV[v] == 0) can_reach_free_v = true;
                    }
                }
            } else { // Node in V
                size_t v = current;
                if (pairV[v] != 0 && dist[pairV[v]] == (size_t)-1) {
                    dist[pairV[v]] = dist[v] + 1;
                    queue[tail++] = pairV[v];
                }
            }
        }

        // If no free vertex in V is reachable, we're done
        if (!can_reach_free_v) break;

        // DFS Phase: Find augmenting paths
        for (size_t i = 1; i <= n; i++) visited[i] = 0;
        for (size_t u = 1; u <= left_size; u++) {
            if (pairU[u] == 0 && DFS(u, &data)) {
                max_matching++;
            }
        }
    }

    // Clean up
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
