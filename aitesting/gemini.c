#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For memset
#include <stdbool.h>
#include <limits.h> // For INT_MAX

// Define NIL (sentinel for unmatched nodes)
#define NIL 0
// Define INF (infinity for distance calculation)
#define INF INT_MAX
// Define UNASSIGNED for partition checking
#define UNASSIGNED -1

// Structure for edges provided by the user (comments ignored for partition logic)
typedef struct {
    size_t u;
    size_t v;
} xedge_t;

// --- Internal Data Structures ---

// Node for adjacency list
typedef struct edge_node {
    size_t to;          // Neighbor vertex
    size_t next_edge;   // Index in adj_edges of next edge from same source
} edge_node_t;

// Structure to hold graph and algorithm state
typedef struct {
    size_t n;               // Number of vertices (1 to n)
    size_t m;               // Number of edges input
    edge_node_t *adj_edges; // Array storing all edge destinations and next pointers
    size_t *adj_head;       // head[i] = index in adj_edges of first edge from node i
    size_t edge_capacity;   // Current allocated capacity for edges
    size_t edge_count;      // Current number of edge entries used (2*actual edges)

    int *partition;         // partition[i] = 0 or 1 if node i is in that partition, UNASSIGNED otherwise
    size_t *match;          // match[i] = j if i is matched with j, else NIL
    int *dist;              // dist[i] = distance from a free node in partition 0
    size_t *queue;          // Queue for BFS traversals
} graph_hk_t;

// --- Helper Functions ---

// Initialize the graph structure
static bool graph_init(graph_hk_t *g, size_t n, size_t m) {
    g->n = n;
    g->m = m; // Store original edge count if needed, though not directly used later
    g->edge_capacity = m * 2; // Start with capacity for all edges (undirected)
    g->edge_count = 0;

    // Allocate memory (+1 for 1-based indexing, NIL = 0)
    g->adj_head = (size_t *)malloc((n + 1) * sizeof(size_t));
    g->adj_edges = (edge_node_t *)malloc(g->edge_capacity * sizeof(edge_node_t));
    g->partition = (int *)malloc((n + 1) * sizeof(int));
    g->match = (size_t *)calloc(n + 1, sizeof(size_t)); // Init to NIL (0)
    g->dist = (int *)malloc((n + 1) * sizeof(int));
    g->queue = (size_t *)malloc((n + 1) * sizeof(size_t));

    // Check for allocation failures
    if (!g->adj_head || !g->adj_edges || !g->partition || !g->match || !g->dist || !g->queue) {
        // Free any successfully allocated memory before returning failure
        free(g->adj_head);
        free(g->adj_edges);
        free(g->partition);
        free(g->match);
        free(g->dist);
        free(g->queue);
        return false; // Indicate initialization failure
    }

    // Initialize adjacency list heads to invalid index (using SIZE_MAX as sentinel)
    memset(g->adj_head, -1, (n + 1) * sizeof(size_t)); // Use -1 bit pattern (SIZE_MAX)
    // Initialize partition array
    for (size_t i = 0; i <= n; ++i) {
        g->partition[i] = UNASSIGNED;
    }
    // Match array already initialized to 0 (NIL) by calloc

    return true;
}

// Free allocated graph resources
static void graph_destroy(graph_hk_t *g) {
    free(g->adj_head);
    free(g->adj_edges);
    free(g->partition);
    free(g->match);
    free(g->dist);
    free(g->queue);
    // Reset pointers to prevent double free if called again
    g->adj_head = NULL;
    g->adj_edges = NULL;
    g->partition = NULL;
    g->match = NULL;
    g->dist = NULL;
    g->queue = NULL;
}


// Add an undirected edge to the adjacency list
static void add_edge_undirected(graph_hk_t *g, size_t u, size_t v) {
    // Add edge u -> v
    // Simple dynamic resizing if needed (though unlikely with initial sizing)
    if (g->edge_count >= g->edge_capacity) {
         // In a real-world scenario, realloc would be needed, but we pre-allocate enough.
         // For simplicity, we assume initial allocation is sufficient. If not, add realloc logic.
         fprintf(stderr, "Warning: Edge capacity exceeded (should not happen with pre-allocation).\n");
         return;
    }
    g->adj_edges[g->edge_count].to = v;
    g->adj_edges[g->edge_count].next_edge = g->adj_head[u];
    g->adj_head[u] = g->edge_count;
    g->edge_count++;

    // Add edge v -> u
     if (g->edge_count >= g->edge_capacity) {
         fprintf(stderr, "Warning: Edge capacity exceeded (should not happen with pre-allocation).\n");
         return;
     }
    g->adj_edges[g->edge_count].to = u;
    g->adj_edges[g->edge_count].next_edge = g->adj_head[v];
    g->adj_head[v] = g->edge_count;
    g->edge_count++;
}

// Check for bipartiteness and assign partitions using BFS
// Returns true if bipartite, false otherwise. Assigns g->partition.
static bool check_bipartite_and_assign_partition(graph_hk_t *g) {
    size_t q_head = 0, q_tail = 0;

    for (size_t start_node = 1; start_node <= g->n; ++start_node) {
        if (g->partition[start_node] == UNASSIGNED) {
            // Start BFS for this connected component
            g->partition[start_node] = 0; // Assign to partition 0
            g->queue[q_tail++] = start_node;

            while (q_head < q_tail) {
                size_t u = g->queue[q_head++];
                int current_partition = g->partition[u];
                int neighbor_partition = 1 - current_partition;

                for (size_t idx = g->adj_head[u]; idx != (size_t)-1; idx = g->adj_edges[idx].next_edge) {
                    size_t v = g->adj_edges[idx].to;

                    if (g->partition[v] == UNASSIGNED) {
                        g->partition[v] = neighbor_partition;
                        g->queue[q_tail++] = v;
                    } else if (g->partition[v] == current_partition) {
                        // Edge connects two nodes in the same partition - not bipartite!
                        return false;
                    }
                    // If partition[v] == neighbor_partition, it's already correctly assigned, do nothing.
                }
            }
        }
    }
    // If we processed all components without conflict, it's bipartite.
    return true;
}


// Breadth-First Search for Hopcroft-Karp (adapted for arbitrary partitions)
// Finds layers starting from free nodes in partition 0.
static bool hk_bfs(graph_hk_t *g) {
    size_t q_head = 0, q_tail = 0;

    // Initialize distances
    for (size_t i = 0; i <= g->n; ++i) {
        g->dist[i] = INF;
    }

    // Enqueue all free nodes from partition 0
    for (size_t u = 1; u <= g->n; ++u) {
        // Ensure the node was assigned a partition (relevant for disconnected graphs if not fully processed)
        if (g->partition[u] == 0 && g->match[u] == NIL) {
            g->dist[u] = 0;
            g->queue[q_tail++] = u;
        }
    }

    g->dist[NIL] = INF; // Sentinel distance for reaching *any* free node in partition 1

    while (q_head < q_tail) {
        size_t u = g->queue[q_head++]; // u is always from partition 0 or matched from 1

        // If we explore beyond the shortest path length found so far, stop for this node
        if (g->dist[u] >= g->dist[NIL]) continue;

        // Explore neighbors v of u (v must be in partition 1)
        for (size_t idx = g->adj_head[u]; idx != (size_t)-1; idx = g->adj_edges[idx].next_edge) {
            size_t v = g->adj_edges[idx].to;
            size_t next_node = g->match[v]; // Node matched with v (must be partition 0 or NIL)

            // Check if the matched node 'next_node' (partner of v) is reachable
            // and hasn't been visited/layered yet (dist == INF).
            // This also handles the case where v is free (next_node == NIL).
            if (g->dist[next_node] == INF) {
                 g->dist[next_node] = g->dist[u] + 1;
                 if (next_node == NIL) {
                     // Found an augmenting path ending at free node v (partition 1)
                     g->dist[NIL] = g->dist[u] + 1; // Update shortest path length
                 } else {
                     // Enqueue the matched node 'next_node' (which is in partition 0)
                     g->queue[q_tail++] = next_node;
                 }
             }
        }
    }

    // Returns true if at least one augmenting path length was found (dist[NIL] updated from INF)
    return g->dist[NIL] != INF;
}


// Depth-First Search for Hopcroft-Karp (adapted for arbitrary partitions)
// Finds an augmenting path starting from node u (assumed partition 0).
static bool hk_dfs(graph_hk_t *g, size_t u) {
    // Base case: If u is NIL, it means we followed a matched edge from a node v
    // where v was previously free. So we have found the end of an augmenting path.
    if (u == NIL) return true;

    // Explore neighbors v of u
    for (size_t idx = g->adj_head[u]; idx != (size_t)-1; idx = g->adj_edges[idx].next_edge) {
        size_t v = g->adj_edges[idx].to;
        size_t next_node = g->match[v]; // Node matched with v

        // Check if edge (u, v) leads to a node 'next_node' on the next layer
        // according to distances computed by BFS.
        // Includes the case where v is free (next_node == NIL), because
        // dist[NIL] holds the shortest augmenting path length.
        if (g->dist[next_node] == g->dist[u] + 1) {
            // Recursively search from the partner node 'next_node'
            if (hk_dfs(g, next_node)) {
                // Augmenting path found - flip matching status for edge (u, v)
                g->match[v] = u;
                g->match[u] = v;
                return true; // Path found and augmented
            }
        }
    }

    // No augmenting path found starting from u in this DFS traversal for this phase
    // Mark u as visited in this phase by setting its distance to INF
    g->dist[u] = INF;
    return false;
}

// --- Public Interface Function ---

/**
 * @brief Computes the maximum cardinality matching in a bipartite graph.
 *
 * Determines the bipartition automatically. If the graph is not bipartite,
 * returns 0. Uses the Hopcroft-Karp algorithm.
 *
 * @param n Total number of vertices (nodes are numbered 1 to n).
 * @param m Number of edges in the input array `e`.
 * @param e Array of edges connecting vertices.
 * @return The size of the maximum cardinality matching, or 0 if the graph
 * is not bipartite or an error occurs.
 */
size_t matching(size_t n, size_t m, xedge_t e[]) {
    if (n == 0) {
        return 0;
    }

    graph_hk_t g;
    size_t matching_size = 0;

    // --- Setup ---
    if (!graph_init(&g, n, m)) {
        fprintf(stderr, "Error: Failed to initialize graph structures.\n");
        return 0; // Indicate error
    }

    // --- Build Adjacency List (Undirected) ---
    for (size_t i = 0; i < m; ++i) {
        // Basic validation
        if (e[i].u > n || e[i].v > n || e[i].u == 0 || e[i].v == 0 || e[i].u == e[i].v) {
             fprintf(stderr, "Warning: Skipping invalid edge (%zu, %zu)\n", e[i].u, e[i].v);
            continue;
        }
        add_edge_undirected(&g, e[i].u, e[i].v);
    }

    // --- Determine Bipartition ---
    if (!check_bipartite_and_assign_partition(&g)) {
        fprintf(stderr, "Error: Input graph is not bipartite.\n");
        graph_destroy(&g);
        return 0; // Cannot perform bipartite matching
    }

    // --- Hopcroft-Karp Algorithm ---
    while (hk_bfs(&g)) { // While shortest augmenting paths exist
        // Try to find augmenting paths starting from each free node in partition 0
        for (size_t u = 1; u <= n; ++u) {
            if (g.partition[u] == 0 && g.match[u] == NIL) {
                if (hk_dfs(&g, u)) {
                    matching_size++; // Increment count if DFS found a path
                }
            }
        }
    }

    // --- Cleanup ---
    graph_destroy(&g);

    return matching_size;
}
