#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <limits.h>
#include <time.h>
#include <stdatomic.h>

typedef struct {
    size_t        u;
    size_t        v;
} xedge_t;

int dfs_la_ts(int** graph, int* graphColSize,int graphSize,int* lookahead,int u,atomic_int* visited,int* pairA,int* pairB){
    int j = 0;
    //int j = lookahead[u];
    for (int i = 0; i < graphColSize[u]; i++, j = (j+1)%graphColSize[u]) {
        int v = graph[u][j];
        if(pairB[v] == -1){
            if (atomic_fetch_add(&visited[v], 1) == 0) {
                //lookahead[u] = (j+1)%graphColSize[u];
                return v;
            }
            if(atomic_fetch_add(&visited[v],1) == 0){
                int index = dfs_la_ts(graph,graphColSize,graphSize,lookahead,pairB[v],visited,pairA,pairB);
                if(index != -1){
                    return index;
                }
            }
        }
    }

    
    for (int i = 0; i < graphColSize[u]; i++) {
        int v = graph[u][i];
        if(atomic_fetch_add(&visited[v],1) == 0){
            int index = dfs_la_ts(graph,graphColSize,graphSize,lookahead,pairB[v],visited,pairA,pairB);
            if(index != -1){
                return index;
            }
        }
    }
    return -1;
}


int parallel_pothen_fan(int** graph, int graphSize, int* graphColSize, int* pairA, int* pairB){
    int* lookahead = (int*)malloc(graphSize * sizeof(int));
    for (int i = 0; i < graphSize; i++) {
        lookahead[i] = 0;
    }
    int path_found = 1;
    atomic_int matchings = 0;
    atomic_int* visited = (atomic_int*)calloc(graphSize, sizeof(atomic_int));
    while(path_found){
        for(int i = 0; i < graphSize; i+=2) {
            atomic_store(&visited[i], 0);
        }
        path_found = 0;
        for(int i = 1; i < graphSize; i+=2){
            int u = dfs_la_ts(graph,graphColSize,graphSize,lookahead,i,visited,pairA,pairB);
            if(u != -1){
                path_found = 1;
                pairB[u] = i;
                pairA[i] = u;
                atomic_fetch_add(&matchings,1);
            }
        }
    }
    return matchings;
}


void matchAndUpdate(int** graph, int* graphColSize, int* deg, int* pairA, int* pairB, bool* visited, int u, int* matchingSize) {
    if (visited[u]) return;
    visited[u] = true;
    for (int j = 0; j < graphColSize[u]; j++) {
        int v = graph[u][j];
        if (!visited[v]) {
            visited[v] = true;
            *matchingSize += 1;
            pairA[u] = v;
            pairB[v] = u;
            for (int i = 0; i < graphColSize[v]; i++) {
                int w = graph[v][i];
                deg[w]--;
                if (deg[w] == 1) {
                    matchAndUpdate(graph, graphColSize, deg, pairA, pairB, visited, w, matchingSize);
                }
            }
            break;
        }
    }
}


int karpSipser(int** graph, int graphSize, int* graphColSize, int* pairA, int* pairB){
    int* deg = (int*)malloc(graphSize * sizeof(int));
    memcpy(deg, graphColSize, graphSize * sizeof(int));
    bool* visited = (bool*)calloc(graphSize, sizeof(bool));
    int* deg1 = (int*)malloc(graphSize * sizeof(int));
    int* deg2 = (int*)malloc(graphSize * sizeof(int));
    int deg1Size = 0;
    int deg2Size = 0;
    for (int i = 1; i < graphSize; i+=2) {
        if (deg[i] == 1) {
            deg1[deg1Size++] = i;
        } else{
            deg2[deg2Size++] = i;
        }
    }
    int matchingSize = 0;
    for(int i = 0; i < deg1Size; i++) {
        matchAndUpdate(graph, graphColSize, deg, pairA, pairB, visited, deg1[i], &matchingSize);
    }
    for(int i = 0; i < deg2Size; i++) {
        matchAndUpdate(graph, graphColSize, deg, pairA, pairB, visited, deg2[i], &matchingSize);
    }
    return matchingSize;
}

bool bfs(int** graph, int graphSize, int* graphColSize, int* pairU, int* pairV, int* dist) {
    int queue[graphSize + 1];
    int front = 0, rear = 0;

    for (int u = 1; u < graphSize; u+=2) {
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
        
        if (dist[u] < dist[graphSize]) { 
            for (int j = 0; j < graphColSize[u]; j++) {
                int v = graph[u][j];

                if (pairV[v] == -1) {
                    dist[graphSize] = dist[u] + 1;
                } 

                else if (dist[pairV[v]] == INT_MAX) {
                    dist[pairV[v]] = dist[u] + 1;
                    queue[rear++] = pairV[v];
                }
            }
        }
    }
    

    return dist[graphSize] != INT_MAX;
}

bool dfs(int** graph, int graphSize, int* graphColSize, int* pairU, int* pairV, int* dist, int u) {
    if (u == -1) return true; 
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

int hopcroftKarp(int** graph, int graphSize, int* graphColSize, int* pairA, int* pairB) {

    int* dist = (int*)malloc((graphSize + 1) * sizeof(int));

    int matchingSize = 0;
    while (bfs(graph, graphSize, graphColSize, pairA, pairB, dist)) {
        for (int u = 1; u < graphSize; u+=2) {
            if (pairA[u] == -1 && dfs(graph, graphSize, graphColSize, pairA, pairB, dist, u)) {
                matchingSize++;
            }
        }
    }
    
    return matchingSize;
}

int maximumBipartiteMatching(int** graph, int graphSize, int* graphColSize, int* matching) {
    int* pairA = (int*)malloc((graphSize+1) * sizeof(int));
    int* pairB = (int*)malloc((graphSize+1) * sizeof(int));
    memset(pairA, -1, (graphSize+1) * sizeof(int));
    memset(pairB, -1, (graphSize+1) * sizeof(int));

    
    clock_t start_time = clock();

    //int init_matchings = karpSipser(graph, graphSize, graphColSize, pairA, pairB);
    int init_matchings = 0;
    // Stop the timer after karpSipser and before hopcroftKarp
    clock_t end_time_1 = clock();
    double time_karpSipser = (double)(end_time_1 - start_time) / CLOCKS_PER_SEC;
    //printf("Time taken by karpSipser: %f seconds\n", time_karpSipser);

    // Now, call hopcroftKarp
    start_time = clock();
    int main_matchings = parallel_pothen_fan(graph, graphSize, graphColSize,pairA, pairB);
    //int main_matchings = hopcroftKarp(graph, graphSize, graphColSize, pairA, pairB);

    // Stop the timer after hopcroftKarp
    clock_t end_time_2 = clock();
    double time_hopcroftKarp = (double)(end_time_2 - start_time) / CLOCKS_PER_SEC;
    //printf("Time taken by hopcroftKarp: %f seconds\n", time_hopcroftKarp);
    //printf("%d %d\n", init_matchings, main_matchings);
    return init_matchings + main_matchings;
}

int** createGraphFromArrays(size_t n, size_t m, xedge_t e[]) {
    int* edgeCounts = (int*)calloc(n, sizeof(int));
    

    for (size_t i = 0; i < m; i++) {
        if (e[i].u < n) edgeCounts[e[i].u]++;
        if (e[i].v < n) edgeCounts[e[i].v]++;
    }

    
    int** graph = (int**)malloc(n * sizeof(int*));
    for (size_t i = 0; i < n; i++) {
        graph[i] = (int*)malloc(edgeCounts[i] * sizeof(int));
    }
    
    int* currentIndex = (int*)calloc(n, sizeof(int));
    

    for (size_t i = 0; i < m; i++) {
        if (e[i].u < n)
            graph[e[i].u][currentIndex[e[i].u]++] = e[i].v;
        if (e[i].v < n)
            graph[e[i].v][currentIndex[e[i].v]++] = e[i].u;
    }

    
    free(edgeCounts);
    free(currentIndex);
    
    return graph;
}

size_t matching(size_t n, size_t m, xedge_t e[]) {
    int* edgeCounts = (int*)calloc(n, sizeof(int));
    for (size_t i = 0; i < m; i++) {
        edgeCounts[e[i].u]++;
        edgeCounts[e[i].v]++;
    }
    int** graph = (int**)malloc((n+1) * sizeof(int*));
    
    for (size_t i = 0; i < (n+1); i++) {
        graph[i] = (int*)malloc(edgeCounts[i] * sizeof(int));

    }
    int* currentIndices = (int*)calloc((n+1), sizeof(int));
    
    for (size_t i = 0; i < m; i++) {
        if (e[i].u < n) {
            graph[e[i].u][currentIndices[e[i].u]++] = e[i].v;
        }
        if (e[i].v < n) {
            graph[e[i].v][currentIndices[e[i].v]++] = e[i].u;
        }
    }
    free(currentIndices);
    int* matchResult = (int*)malloc((n+1) * sizeof(int));

    
    memset(matchResult, -1, (n+1) * sizeof(int));
   
    size_t maxMatching = maximumBipartiteMatching(graph, n, edgeCounts, matchResult);

    
    return maxMatching;
}

/*
void runTestCaseFromFile(const char* filename, size_t setASize, size_t setBSize) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }
    
    size_t edgeCount = 0;
    int u, v;
    while (fscanf(file, "%d %d", &u, &v) == 2) {
        if (u < setASize && v < setBSize) {
            edgeCount++;
        }
    }
    
    rewind(file);
    
    xedge_t* edges = (xedge_t*)malloc(edgeCount * sizeof(xedge_t));
    
    size_t edgeIndex = 0;
    while (fscanf(file, "%d %d", &u, &v) == 2) {
        if (u < setASize && v < setBSize) {
            edges[edgeIndex].u = u;
            edges[edgeIndex].v = v;
            edgeIndex++;
        }
    }
    fclose(file);
    int matchingSize = matching(setASize, edgeCount, edges);
    printf("Matching size: %zu\n", matchingSize);

}

int main() {
    runTestCaseFromFile("mag_smaller.txt", 2000000, 200000);
    return 0;
}*/
