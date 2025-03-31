#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#define BIG 1410065408



bool dfs(int u, int setASize, int* pairA, int* pairB, int* dist, int* visited, int iteration, int** graph, int* graphColSize) {
    if (u == setASize) return true;
    
    visited[u] = iteration; 

    for (int i = 0; i < graphColSize[u]; i++) {
        int v = graph[u][i];
        int nextU = (pairB[v] == -1) ? setASize : pairB[v];

        if (visited[nextU] == iteration) continue;
        if (dist[nextU] == dist[u] + 1) {
            if (dfs(nextU, setASize, pairA, pairB, dist, visited, iteration, graph, graphColSize)) {
                    pairB[v] = u;
                    pairA[u] = v;
                return true;
            }
        }
    }
    return false;
}

int maximumBipartiteMatching(int** graph, int graphSize, int* graphColSize, int* matching) {
    int setASize = graphSize;
    
    int* pairA = (int*)calloc(setASize, sizeof(int));
    int* pairB = (int*)calloc(graphSize, sizeof(int)); 
    int* dist = (int*)malloc((setASize + 1) * sizeof(int));
    int* visited = (int*)calloc(setASize, sizeof(int));
    
    memset(pairA, -1, setASize * sizeof(int));
    memset(pairB, -1, graphSize * sizeof(int));
    
    int matchingSize = 0;
    int iteration = 0;


    int* level = (int*)malloc((setASize + 1) * sizeof(int));

    while (1) {
        memset(dist,BIG,(setASize + 1) * sizeof(int));
        for (int i = 0; i <= setASize; i++) {
            dist[i] = BIG;
        }
        
        int levelStart = 0;
        int levelEnd = 0;

        for (int u = 0; u < setASize; u++) {
            if (pairA[u] == -1) {
                dist[u] = 0;
                level[levelEnd++] = u;
            }
        }

        bool pathFound = false;
        while (levelStart < levelEnd && !pathFound) {
            int levelSize = levelEnd - levelStart;
            int newLevelEnd = levelEnd;
            for (int i = 0; i < levelSize; i++) {
                int u = level[levelStart + i];

                if (u == setASize) {
                    pathFound = true;
                    continue;
                }

                for (int j = 0; j < graphColSize[u]; j++) {
                    int v = graph[u][j];
                    int nextU = (pairB[v] == -1) ? setASize : pairB[v];

                    if (dist[nextU] == BIG) {
                        dist[nextU] = dist[u] + 1;
                        level[newLevelEnd++] = nextU;
                    }
                }
            }

            levelStart += levelSize;
            levelEnd = newLevelEnd;

        }
        
        if (dist[setASize] == BIG) break;
        iteration++;
        for (int u = 0; u < setASize; u++) {
            if (pairA[u] == -1) {
                if (dfs(u, setASize, pairA, pairB, dist, visited, iteration, graph, graphColSize)) {
                    matchingSize++;
                }
            }
        }
    }
    
    for (int i = 0; i < setASize; i++) {
        matching[i] = pairA[i];
    }
    
    free(pairA);
    free(pairB);
    free(dist);
    free(visited);
    
    return matchingSize;
}

#define MAX_NODES 200000

int** createGraphFromArrays(int setASize, int** edgeLists, int* edgeCounts) {
    int** graph = (int**)malloc(setASize * sizeof(int*));
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
