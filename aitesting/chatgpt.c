#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define INF INT_MAX

// Structure for dynamic adjacency list (for vertices on the left side)
typedef struct {
    int *neighbors;   // Array of neighbors (right side vertices)
    int size;         // Number of neighbors currently stored
    int capacity;     // Allocated capacity of the array
} AdjList;

// Global arrays for matching and distances.
// pairU[u] is the right vertex matched with left vertex u (or 0 if none)
// pairV[v] is the left vertex matched with right vertex v (or 0 if none)
// dist[u] is used during BFS.
int *pairU, *pairV, *dist;

// The number of vertices in the left and right sets.
int nLeft, nRight;

// Array of adjacency lists for the left vertices (indexed 1..nLeft)
AdjList *adj;

// Queue structure used in BFS
typedef struct {
    int *data;
    int front, back, capacity;
} Queue;

void initQueue(Queue *q, int cap) {
    q->data = (int *)malloc(cap * sizeof(int));
    q->capacity = cap;
    q->front = 0;
    q->back = 0;
}

void freeQueue(Queue *q) {
    free(q->data);
}

int isEmpty(Queue *q) {
    return q->front == q->back;
}

void enqueue(Queue *q, int val) {
    if(q->back == q->capacity) {
        q->capacity *= 2;
        q->data = (int *)realloc(q->data, q->capacity * sizeof(int));
    }
    q->data[q->back++] = val;
}

int dequeue(Queue *q) {
    return q->data[q->front++];
}

// BFS routine of Hopcroft-Karp. Returns 1 if there is an augmenting path.
int bfs() {
    Queue q;
    initQueue(&q, nLeft);
    
    // Initialize distances for all left vertices
    for (int u = 1; u <= nLeft; u++) {
        if (pairU[u] == 0) {
            dist[u] = 0;
            enqueue(&q, u);
        } else {
            dist[u] = INF;
        }
    }
    
    // Set the distance for the dummy vertex 0
    dist[0] = INF;
    
    // Process the queue
    while (!isEmpty(&q)) {
        int u = dequeue(&q);
        if (dist[u] < dist[0]) {
            // Traverse all neighbors v of u
            for (int i = 0; i < adj[u].size; i++) {
                int v = adj[u].neighbors[i];
                if (dist[pairV[v]] == INF) {
                    dist[pairV[v]] = dist[u] + 1;
                    enqueue(&q, pairV[v]);
                }
            }
        }
    }
    
    freeQueue(&q);
    return dist[0] != INF;
}

// DFS routine to search for augmenting paths starting from left vertex u.
int dfs(int u) {
    if (u != 0) {
        for (int i = 0; i < adj[u].size; i++) {
            int v = adj[u].neighbors[i];
            if (dist[pairV[v]] == dist[u] + 1 && dfs(pairV[v])) {
                pairV[v] = u;
                pairU[u] = v;
                return 1;
            }
        }
        dist[u] = INF;
        return 0;
    }
    return 1;
}

// Hopcroft-Karp algorithm. Returns the maximum matching.
int hopcroftKarp() {
    int matching = 0;
    while (bfs()) {
        for (int u = 1; u <= nLeft; u++) {
            if (pairU[u] == 0 && dfs(u))
                matching++;
        }
    }
    return matching;
}

// Function to add an edge from left vertex u to right vertex v.
void addEdge(int u, int v) {
    // Expand the adjacency list if needed.
    if (adj[u].size == adj[u].capacity) {
        int newCapacity = (adj[u].capacity == 0 ? 2 : adj[u].capacity * 2);
        adj[u].neighbors = (int *)realloc(adj[u].neighbors, newCapacity * sizeof(int));
        adj[u].capacity = newCapacity;
    }
    adj[u].neighbors[adj[u].size++] = v;
}

int main() {
    int u, v;
    int capacityEdges = 100;
    int countEdges = 0;
    
    // Temporary storage for edges
    int (*edges)[2] = malloc(capacityEdges * sizeof(*edges));
    
    // Determine the maximum vertex id on each side.
    nLeft = 0;
    nRight = 0;
    
    // Read input edges from stdin
    while (scanf("%d %d", &u, &v) == 2) {
        if (countEdges == capacityEdges) {
            capacityEdges *= 2;
            edges = realloc(edges, capacityEdges * sizeof(*edges));
        }
        edges[countEdges][0] = u;
        edges[countEdges][1] = v;
        countEdges++;
        if (u > nLeft) nLeft = u;
        if (v > nRight) nRight = v;
    }
    
    // Allocate memory for the left adjacency lists (1-indexed)
    adj = (AdjList *)malloc((nLeft + 1) * sizeof(AdjList));
    for (int i = 0; i <= nLeft; i++) {
        adj[i].neighbors = NULL;
        adj[i].size = 0;
        adj[i].capacity = 0;
    }
    
    // Populate the graph with edges
    for (int i = 0; i < countEdges; i++) {
        u = edges[i][0];
        v = edges[i][1];
        addEdge(u, v);
    }
    free(edges);
    
    // Allocate arrays for pairU, pairV, and dist.
    pairU = (int *)calloc(nLeft + 1, sizeof(int)); // indices 1..nLeft; 0 used as NIL.
    pairV = (int *)calloc(nRight + 1, sizeof(int)); // indices 1..nRight.
    dist  = (int *)malloc((nLeft + 1) * sizeof(int));
    
    // Compute the maximum matching.
    int matching = hopcroftKarp();
    
    // Output the number of matches.
    printf("%d\n", matching);
    
    // Free allocated memory.
    for (int i = 0; i <= nLeft; i++) {
        free(adj[i].neighbors);
    }
    free(adj);
    free(pairU);
    free(pairV);
    free(dist);
    
    return matching;
}

