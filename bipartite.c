#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <limits.h>

bool bfs(int** graph, int graphSize, int* graphColSize, int* pairU, int* pairV, int* dist) {
    int queue[graphSize];
    int front = 0, rear = 0;
    
    for (int u = 0; u < graphSize; u++) {
        if (pairU[u] == -1) {
            dist[u] = 0;
            queue[rear++] = u; 
        } else {
            dist[u] = INT_MAX;
        }
    }
    
    dist[graphSize] = INT_MAX;
    

    while (front < rear) {
        int u = queue[front++];
        
        if (u != graphSize) {
            for (int j = 0; j < graphColSize[u]; j++) {
                int v = graph[u][j];
                
                if (pairV[v] == -1) {
                    dist[graphSize] = dist[u] + 1;
                } else if (dist[pairV[v]] == INT_MAX) {
                    dist[pairV[v]] = dist[u] + 1;
                    queue[rear++] = pairV[v];
                }
            }
        }
    }
    
    return dist[graphSize] != INT_MAX;
}

bool dfs(int** graph, int graphSize, int* graphColSize, int* pairU, int* pairV, int* dist, int u) {
    if (u == graphSize) return true;
    for (int j = 0; j < graphColSize[u]; j++) {
        int v = graph[u][j];
        
        if (pairV[v] == -1 || (dist[pairV[v]] == dist[u] + 1 && dfs(graph, graphSize, graphColSize, pairU, pairV, dist, pairV[v]))) {
            pairV[v] = u;
            pairU[u] = v;
            return true;
        }
    }
    
    dist[u] = INT_MAX;
    return false;
}


int hopcroftKarp(int** graph, int graphSize, int* graphColSize, int* matching) {

    int* pairA = (int*)malloc(graphSize * sizeof(int));
    int* pairB = (int*)malloc(graphSize * sizeof(int));
    int* dist = (int*)malloc((graphSize + 1) * sizeof(int));

    memset(pairA, -1, graphSize * sizeof(int));
    memset(pairB, -1, graphSize * sizeof(int));
    
    int matchingSize = 0;
    while (bfs(graph, graphSize, graphColSize, pairA, pairB, dist)) {
        for (int u = 0; u < graphSize; u++) {
            if (pairA[u] == -1 && dfs(graph, graphSize, graphColSize, pairA, pairB, dist, u)) {
                matchingSize++;
            }
        }
    }
    #pragma omp parallel for
    for (int u = 0; u < graphSize; u++) {
        matching[u] = pairA[u];
    }

    
    return matchingSize;
}



int maximumBipartiteMatching(int** graph, int graphSize, int* graphColSize, int* matching) {
    return hopcroftKarp(graph, graphSize, graphColSize, matching);
}



int** createGraphFromArrays(int setASize, int** edgeLists, int* edgeCounts) {
    int** graph = (int**)malloc(setASize * sizeof(int*));
    #pragma omp parallel for
    for (int i = 0; i < setASize; i++) {
        graph[i] = (int*)malloc(edgeCounts[i] * sizeof(int));
        for (int j = 0; j < edgeCounts[i]; j++) {
            graph[i][j] = edgeLists[i][j];
        }
    }
    return graph;
}

#include <time.h>

void runTestCaseFromFile(const char* filename, int setASize, int setBSize) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    int** edgeLists = (int**)malloc(setASize * sizeof(int*));
    int* edgeCounts = (int*)calloc(setASize, sizeof(int));

    int u, v;
    while (fscanf(file, "%d %d", &u, &v) == 2) {
        if (u >= setASize || v >= setBSize) continue; 

        edgeCounts[u]++;
        edgeLists[u] = (int*)realloc(edgeLists[u], edgeCounts[u] * sizeof(int));
        edgeLists[u][edgeCounts[u] - 1] = v;
    }
    fclose(file);

    int** graph = createGraphFromArrays(setASize, edgeLists, edgeCounts);


    int* matching = (int*)malloc(setASize * sizeof(int));
    memset(matching, -1, setASize * sizeof(int));
    double start = omp_get_wtime();
    int maxMatching = maximumBipartiteMatching(graph, setASize, edgeCounts, matching);
    double end = omp_get_wtime();
    double time_taken = end - start;
    printf("Function took %f seconds to execute\n", time_taken);
    printf("Matching size: %d\n", maxMatching);


}



int main() {
    runTestCaseFromFile("lastfm.txt", 2000000, 200000);
    return 0;
}
