#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <limits.h>
#include <time.h>
#include <stdatomic.h>
#include <assert.h>

typedef struct {
    size_t        u;
    size_t        v;
} xedge_t;

typedef struct {
    int u;
    int i;
} u_i_t;

// uses two stacks, node u, and it's i'eth edge
static inline bool dfs_la_ts(int** graph, int* graphColSize,int graphSize,int* lookahead,int start, atomic_int* visited,int* pairA,int* pairB, u_i_t *stack, bool order){
    int stack_at = 0;
    int u = start;
    WORK_ON_NEW_U:;
    int i = order ? 0 : (graphColSize[u] - 1);
    // take from stack
    assert(u != 0); // should not be 0

    // TODO: do lookahead thing if we use stack-mechanism instead of recursive.
    for (int j = lookahead[u]; j < graphColSize[u]; j++) {
        int v = graph[u][j];
        if(pairB[v] == 0){ // unmatched v. Can match them together
            if (atomic_fetch_add(&visited[v], 1) == 0) {
                lookahead[u] = j+1;
                // found unmatched v, that wasn't visited. Let's match with this one.
                pairA[u] = v;
                pairB[v] = u;

                // go backwards through stack and match things
                while (stack_at > 0) {
                    stack_at--;
                    u = stack[stack_at].u;
                    assert(u != 0);
                    i = stack[stack_at].i;
                    const int v = graph[u][i];
                    pairA[u] = v;
                    pairB[v] = u;
                }
                return 1;
            }
        }
    }
    lookahead[u] = graphColSize[u]; // stop looking up things in future from this node
    CONTINUE_ON_SECOND_LOOP:;
    for (; order ? (i < graphColSize[u]) : (i >= 0); order ? i++ : i--) { // TODO: fariness
        int v = graph[u][i];
        if(atomic_fetch_add(&visited[v],1) == 0){
            // don't do recursion, instead, push current state on stack, and continue on pairB[v] ()
            int new_u = pairB[v];
            assert(new_u != 0); // v should be matched, because we already did lookahead before this.
            // push (u, v), work then on new_u
            assert(u != 0);
            stack[stack_at].u = u; stack[stack_at].i = i; stack_at++;
            u = new_u;
            goto WORK_ON_NEW_U;
        }
    }
    // we found no continuation from here. Continue from the previous u, at the i it was at in the second loop
    if (stack_at == 0) {
        return 0;
    } // nothing left
    stack_at--;
    u = stack[stack_at].u;
    i = stack[stack_at].i + (order ? (1) : -1); // after ourselves
    goto CONTINUE_ON_SECOND_LOOP;
}


static int parallel_pothen_fan(int** graph, int graphSize, int* graphColSize, int* pairA, int* pairB){
    int* lookahead = (int*)calloc(graphSize, sizeof(int));

    // allocate stacks
    int max_threads = 32;
    omp_set_num_threads(max_threads);
    u_i_t *stacks = calloc(graphSize * max_threads, sizeof(*stacks));

    int path_found = 1;
    long matchings = 0;
    bool order = false;
    // FOR V NODES
    atomic_int* visited = calloc(graphSize, sizeof(atomic_int));
    while(path_found){
        memset(visited, 0, sizeof(*visited)*graphSize);
        path_found = 0;
        order = !order;
        #pragma omp parallel for reduction(+:matchings)
        for(int i = 1; i < graphSize; i+=2){ // UNMATCHED u VERTICES
            // SAFETY: this read is okay. Because no unmatched u's will be visited from searches starting on other unmatches u's. only matched u's will be found.
            if (pairA[i] != 0) continue;
            int t_id = omp_get_thread_num();
            // printf("%d\n", t_id);
            u_i_t *stack = &stacks[graphSize * t_id];
            // these are inlined, so get two monomorphised versions
            int res;
            if (order) {
                res = dfs_la_ts(graph,graphColSize,graphSize,lookahead,i,visited,pairA,pairB, stack, true);
            }
            else {
                res = dfs_la_ts(graph,graphColSize,graphSize,lookahead,i,visited,pairA,pairB, stack, false);
            }
            if(res != 0){
                path_found = 1;
                matchings++;
            }
        }
    }
    free(lookahead);
    free(visited);
    free(stacks);
    return matchings;
}


static void matchAndUpdate(int** graph, int* graphColSize, int* deg, int* pairA, int* pairB, bool* visited, int u, int* matchingSize) {
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


static int karpSipser(int** graph, int graphSize, int* graphColSize, int* pairA, int* pairB){
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

static int maximumBipartiteMatching(int** graph, int graphSize, int* graphColSize) {
    graphSize++;
    int* pairA = (int*)malloc((graphSize) * sizeof(int));
    int* pairB = (int*)malloc((graphSize) * sizeof(int));
    memset(pairA, 0, (graphSize) * sizeof(int));
    memset(pairB, 0, (graphSize) * sizeof(int));


    clock_t start_time = clock();

    // int init_matchings = karpSipser(graph, graphSize, graphColSize, pairA, pairB);
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

static int** createGraphFromArrays(size_t n, size_t m, xedge_t e[]) {
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

    size_t maxMatching = maximumBipartiteMatching(graph, n, edgeCounts);

    return maxMatching;
}
