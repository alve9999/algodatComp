#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

typedef struct {
    int u;
    int v;
} Edge;

int compare_ints(const void *a, const void *b) {
    int arg1 = *(const int *)a;
    int arg2 = *(const int *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

int dfs(int u, int **adj, int *adj_sizes, int *pair_u, int *pair_v, int *dist) {
    for (int i = 0; i < adj_sizes[u]; i++) {
        int v = adj[u][i];
        if (pair_v[v] == -1 || (dist[pair_v[v]] == dist[u] + 1 && dfs(pair_v[v], adj, adj_sizes, pair_u, pair_v, dist))) {
            pair_u[u] = v;
            pair_v[v] = u;
            dist[u] = INT_MAX;
            return 1;
        }
    }
    dist[u] = INT_MAX;
    return 0;
}

int main() {
    Edge *edges = NULL;
    int num_edges = 0;
    int capacity_edges = 0;
    int u, v;

    // Read all edges
    while (scanf("%d %d", &u, &v) == 2) {
        if (num_edges >= capacity_edges) {
            capacity_edges = capacity_edges == 0 ? 1 : capacity_edges * 2;
            edges = realloc(edges, capacity_edges * sizeof(Edge));
            if (!edges) {
                fprintf(stderr, "Memory error\n");
                return 1;
            }
        }
        edges[num_edges].u = u;
        edges[num_edges].v = v;
        num_edges++;
    }

    if (num_edges == 0) {
        printf("0\n");
        return 0;
    }

    // Collect left and right nodes
    int *left_temp = malloc(num_edges * sizeof(int));
    int *right_temp = malloc(num_edges * sizeof(int));
    for (int i = 0; i < num_edges; i++) {
        left_temp[i] = edges[i].u;
        right_temp[i] = edges[i].v;
    }

    // Process left nodes
    qsort(left_temp, num_edges, sizeof(int), compare_ints);
    int left_unique_size = 0;
    for (int i = 0; i < num_edges; i++) {
        if (i == 0 || left_temp[i] != left_temp[i - 1]) {
            left_temp[left_unique_size++] = left_temp[i];
        }
    }
    int *left_nodes = malloc(left_unique_size * sizeof(int));
    memcpy(left_nodes, left_temp, left_unique_size * sizeof(int));
    free(left_temp);

    // Process right nodes
    qsort(right_temp, num_edges, sizeof(int), compare_ints);
    int right_unique_size = 0;
    for (int i = 0; i < num_edges; i++) {
        if (i == 0 || right_temp[i] != right_temp[i - 1]) {
            right_temp[right_unique_size++] = right_temp[i];
        }
    }
    int *right_nodes = malloc(right_unique_size * sizeof(int));
    memcpy(right_nodes, right_temp, right_unique_size * sizeof(int));
    free(right_temp);

    // Build adjacency lists
    int U = left_unique_size;
    int V = right_unique_size;
    int **adj = malloc(U * sizeof(int *));
    int *adj_sizes = calloc(U, sizeof(int));

    int *counts = calloc(U, sizeof(int));
    for (int i = 0; i < num_edges; i++) {
        int u_val = edges[i].u;
        int *found = bsearch(&u_val, left_nodes, left_unique_size, sizeof(int), compare_ints);
        int u_index = found - left_nodes;
        counts[u_index]++;
    }

    for (int u = 0; u < U; u++) {
        adj[u] = malloc(counts[u] * sizeof(int));
    }

    memset(counts, 0, U * sizeof(int));
    for (int i = 0; i < num_edges; i++) {
        int u_val = edges[i].u;
        int *found_u = bsearch(&u_val, left_nodes, left_unique_size, sizeof(int), compare_ints);
        int u_index = found_u - left_nodes;
        int v_val = edges[i].v;
        int *found_v = bsearch(&v_val, right_nodes, right_unique_size, sizeof(int), compare_ints);
        int v_index = found_v - right_nodes;
        adj[u_index][counts[u_index]++] = v_index;
    }

    for (int u = 0; u < U; u++) {
        int size = counts[u];
        if (size == 0) {
            adj_sizes[u] = 0;
            continue;
        }
        qsort(adj[u], size, sizeof(int), compare_ints);
        int new_size = 0;
        for (int i = 0; i < size; i++) {
            if (i == 0 || adj[u][i] != adj[u][i - 1]) {
                adj[u][new_size++] = adj[u][i];
            }
        }
        adj[u] = realloc(adj[u], new_size * sizeof(int));
        adj_sizes[u] = new_size;
    }
    free(counts);

    // Initialize pair arrays
    int *pair_u = malloc(U * sizeof(int));
    int *pair_v = malloc(V * sizeof(int));
    for (int i = 0; i < U; i++) pair_u[i] = -1;
    for (int j = 0; j < V; j++) pair_v[j] = -1;

    int *dist = malloc(U * sizeof(int));
    int result = 0;

    while (1) {
        // BFS to find layered paths
        int *queue = malloc(U * sizeof(int));
        int front = 0, rear = 0;
        int found = 0;

        for (int u = 0; u < U; u++) {
            if (pair_u[u] == -1) {
                dist[u] = 0;
                queue[rear++] = u;
            } else {
                dist[u] = INT_MAX;
            }
        }

        while (front < rear) {
            int u = queue[front++];
            for (int i = 0; i < adj_sizes[u]; i++) {
                int v = adj[u][i];
                if (pair_v[v] == -1) {
                    found = 1;
                } else if (dist[pair_v[v]] == INT_MAX) {
                    dist[pair_v[v]] = dist[u] + 1;
                    queue[rear++] = pair_v[v];
                }
            }
        }

        free(queue);
        if (!found) break;

        // Perform DFS for each unmatched left node
        for (int u = 0; u < U; u++) {
            if (pair_u[u] == -1) {
                if (dfs(u, adj, adj_sizes, pair_u, pair_v, dist)) {
                    result++;
                }
            }
        }
    }

    printf("%d\n", result);

    // Cleanup
    free(pair_u);
    free(pair_v);
    free(dist);
    for (int u = 0; u < U; u++) {
        free(adj[u]);
    }
    free(adj);
    free(adj_sizes);
    free(left_nodes);
    free(right_nodes);
    free(edges);

    return 0;
}
