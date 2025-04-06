#include <stdlib.h>
#include <stddef.h>

typedef struct {
    size_t u;    /* odd.        */
    size_t v;    /* even.    */
} xedge_t;

static int dfs(size_t u, size_t left_size, size_t **adj, size_t *degree, size_t *Pair_U, size_t *Pair_V, size_t *dist) {
    if (u != left_size) { // u is in U
        for (size_t i = 0; i < degree[u]; i++) {
            size_t v = adj[u][i];
            size_t v_idx = v - left_size;
            size_t u_prime = Pair_V[v_idx];
            if (u_prime == left_size) { // v is free
                if (dist[left_size] == dist[u] + 1) {
                    if (dfs(left_size, left_size, adj, degree, Pair_U, Pair_V, dist)) {
                        Pair_U[u] = v;
                        Pair_V[v_idx] = u;
                        return 1;
                    }
                }
            } else {
                if (dist[u_prime] == dist[u] + 1) {
                    if (dfs(u_prime, left_size, adj, degree, Pair_U, Pair_V, dist)) {
                        Pair_U[u] = v;
                        Pair_V[v_idx] = u;
                        return 1;
                    }
                }
            }
        }
        dist[u] = (size_t)-1; // mark as unreachable in next BFS
        return 0;
    }
    return 1; // reached NIL, augmenting path found
}

size_t matching(size_t n, size_t m, xedge_t e[]) {
    size_t left_size = n / 2;
    size_t right_size = n / 2;
    if (n == 0 || left_size + right_size != n) {
        return 0; // invalid input
    }

    // Build adjacency lists for U (left partition).
    size_t *degree = calloc(left_size, sizeof(size_t));
    for (size_t i = 0; i < m; i++) {
        size_t u = e[i].u;
        if (u < left_size) {
            degree[u]++;
        }
    }

    size_t **adj = malloc(left_size * sizeof(size_t *));
    for (size_t u = 0; u < left_size; u++) {
        adj[u] = malloc(degree[u] * sizeof(size_t));
    }

    size_t *index = calloc(left_size, sizeof(size_t));
    for (size_t i = 0; i < m; i++) {
        size_t u = e[i].u;
        size_t v = e[i].v;
        if (u < left_size) {
            adj[u][index[u]++] = v;
        }
    }
    free(index);

    // Initialize Pair_U and Pair_V.
    size_t NIL = left_size;
    size_t *Pair_U = malloc(left_size * sizeof(size_t));
    size_t *Pair_V = malloc(right_size * sizeof(size_t));
    for (size_t u = 0; u < left_size; u++) {
        Pair_U[u] = NIL;
    }
    for (size_t v_idx = 0; v_idx < right_size; v_idx++) {
        Pair_V[v_idx] = NIL;
    }

    size_t *dist = malloc((left_size + 1) * sizeof(size_t));
    size_t *queue = malloc(left_size * sizeof(size_t));
    size_t matching_size = 0;

    while (1) {
        // BFS to find layers.
        size_t q_head = 0, q_tail = 0;
        for (size_t u = 0; u < left_size; u++) {
            if (Pair_U[u] == NIL) {
                dist[u] = 0;
                queue[q_tail++] = u;
            } else {
                dist[u] = (size_t)-1;
            }
        }
        dist[NIL] = (size_t)-1;

        while (q_head < q_tail) {
            size_t u = queue[q_head++];
            if (dist[u] < dist[NIL]) {
                for (size_t i = 0; i < degree[u]; i++) {
                    size_t v = adj[u][i];
                    size_t v_idx = v - left_size;
                    size_t u_prime = Pair_V[v_idx];
                    if (u_prime == NIL) {
                        if (dist[NIL] == (size_t)-1) {
                            dist[NIL] = dist[u] + 1;
                        }
                    } else {
                        if (dist[u_prime] == (size_t)-1) {
                            dist[u_prime] = dist[u] + 1;
                            queue[q_tail++] = u_prime;
                        }
                    }
                }
            }
        }

        if (dist[NIL] == (size_t)-1) {
            break; // no more augmenting paths
        }

        // Run DFS to find augmenting paths.
        size_t augmenting_paths = 0;
        for (size_t u = 0; u < left_size; u++) {
            if (Pair_U[u] == NIL) {
                if (dfs(u, left_size, adj, degree, Pair_U, Pair_V, dist)) {
                    augmenting_paths++;
                }
            }
        }

        if (augmenting_paths == 0) {
            break;
        }

        matching_size += augmenting_paths;
    }

    // Free allocated memory.
    free(dist);
    free(queue);
    free(Pair_U);
    free(Pair_V);
    for (size_t u = 0; u < left_size; u++) {
        free(adj[u]);
    }
    free(adj);
    free(degree);

    return matching_size;
}
