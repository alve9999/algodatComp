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

typedef struct {
    union {
        atomic_bool b;
    };
} visit_t;

// uses two stacks, node u, and it's i'eth edge
static inline bool dfs_la_ts(int* graph, int* graphOffset,int graphSize,int* lookahead,int start, visit_t* visited,int* v_pair, u_i_t *stack, bool order){
    int stack_at = 0;
    int u = start;
    WORK_ON_NEW_U:;
    int i = order ? graphOffset[u] : (graphOffset[u+1] - 1);
    // take from stack
    assert(u != 0); // should not be 0
    // TODO: do lookahead thing if we use stack-mechanism instead of recursive.
    for (int j = lookahead[u]; j < graphOffset[u+1]; j++) {
        int v = graph[j];
        if(v_pair[v] == 0){ // unmatched v. Can match them together
            if (!atomic_load_explicit(&visited[v].b, memory_order_relaxed)) {
                if (atomic_exchange_explicit(&visited[v].b,true, memory_order_relaxed) == 0){
                    lookahead[u] = j+1;
                    // found unmatched v, that wasn't visited. Let's match with this one.
                    // pairA[u] = v; // Not needed!
                    v_pair[v] = u;

                    // go backwards through stack and match things
                    while (stack_at > 0) {
                        stack_at--;
                        u = stack[stack_at].u;
                        i = stack[stack_at].i;
                        int v = graph[i];
                        assert(u != 0);
                        // pairA[u] = v; // not necessary, already matched to something, won't be checked
                        v_pair[v] = u;
                    }
                    return 1;
                }
            }
        }
    }
    lookahead[u] = graphOffset[u+1]; // stop looking up things in future from this node
    CONTINUE_ON_SECOND_LOOP:;
    for (; order ? (i < graphOffset[u+1]) : (i >= graphOffset[u]); order ? i++ : i--) { // TODO: fariness
        int v = graph[i];
        if (!atomic_load_explicit(&visited[v].b, memory_order_relaxed)) {
            if (atomic_exchange_explicit(&visited[v].b,true, memory_order_relaxed) == 0){
                // don't do recursion, instead, push current state on stack, and continue on pairB[v] ()
                int new_u = v_pair[v];
                assert(new_u != 0); // v should be matched, because we already did lookahead before this.
                // push (u, v), work then on new_u
                assert(u != 0);
                stack[stack_at].u = u; stack[stack_at].i = i; stack_at++;

                u = new_u;
                goto WORK_ON_NEW_U;
            }
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

static int parallel_pothen_fan(int* graph, int graphSize, int* graphOffset, int* v_pair, bool *u_matched){
    static bool inited = false;
    static int* lookahead;
    static u_i_t *stacks;
    static visit_t* visited;
    const int max_threads = 32;
    if (!inited) {
        inited = true;
        lookahead = (int*)calloc(graphSize, sizeof(int));
        stacks = calloc(graphSize * max_threads, sizeof(*stacks));
        visited = calloc(graphSize, sizeof(*visited));
    }
    memcpy(lookahead, graphOffset, (size_t)graphSize * sizeof(*lookahead)); // Copy graphOffset to lookahead

    // allocate stacks
    omp_set_num_threads(max_threads);

    int path_found = 1;
    long matchings = 0;
    bool order = false;
    // FOR V NODES
    while(path_found){
        memset(visited, 0, sizeof(*visited)*graphSize);
        path_found = 0;
        order = !order;
        int worklist_size = 0;
        static int* worklist;
        static int* local_matchings;;
        if (!worklist) worklist = malloc(graphSize * sizeof(int));
        for (int i = 0; i < (graphSize >> 1); ++i) {
            if (!u_matched[i]) {
                worklist[worklist_size++] = 1 + (i << 1);
            }
        }

        #pragma omp parallel for schedule(guided) reduction(+:matchings) 
        for(int idx = 0; idx < worklist_size; idx+=1){ // UNMATCHED u VERTICES
            int i = worklist[idx];
            // SAFETY: this read is okay. Because no unmatched u's will be visited from searches starting on other unmatches u's. only matched u's will be found.
            int t_id = omp_get_thread_num();
            u_i_t *stack = &stacks[graphSize * t_id];
            // these are inlined, so get two monomorphised versions
            int res;
            if (order) res = dfs_la_ts(graph,graphOffset,graphSize,lookahead,i,visited,v_pair,  stack, true);
            else       res = dfs_la_ts(graph,graphOffset,graphSize,lookahead,i,visited,v_pair,  stack, false);
            if(res != 0) {
                u_matched[i >> 1] = 1;
                path_found = 1;
                matchings++;
            }
        }
    }
    // free(lookahead);
    // free(visited);
    // free(stacks);
    return matchings;
}

static void matchAndUpdate(int* graph, int* graphOffset, int* deg, int* pair, bool* visited, int u, int* matchingSize) {
    if (visited[u]) return;
    visited[u] = true;
    for (int j = graphOffset[u]; j < (graphOffset[u+1]); j++) {
        int v = graph[j];
        if (!visited[v]) {
            visited[v] = true;
            *matchingSize += 1;
            pair[u] = 1;
            pair[v] = u;
            for (int i = graphOffset[v]; i < graphOffset[v+1]; i++) {
                int w = graph[i];
                deg[w]--;
                if (deg[w] == 1) {
                    matchAndUpdate(graph, graphOffset, deg, pair, visited, w, matchingSize);
                }
            }
            break;
        }
    }
}

static int karpSipser(int* graph, int graphSize, int* graphOffset, int* pair){
    int* deg = (int*)malloc(graphSize * sizeof(int));
    for (int i = 0; i < graphSize; i++) {
        deg[i] = graphOffset[i+1] - graphOffset[i];
    }
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
        matchAndUpdate(graph, graphOffset, deg, pair, visited, deg1[i], &matchingSize);
    }
    for(int i = 0; i < deg2Size; i++) {
        matchAndUpdate(graph, graphOffset, deg, pair, visited, deg2[i], &matchingSize);
    }
    return matchingSize;
}


static int maximumBipartiteMatching(int* graph, int graphSize, int* graphOffset) {
    graphSize++;
    static bool inited = false;
    static int *v_pair;
    static bool *u_matched;
    if (!inited) {
        v_pair = malloc((graphSize) * sizeof(*v_pair));
        u_matched = malloc((graphSize) * sizeof(*u_matched));
    }
    memset(v_pair, 0, (graphSize) * sizeof(*v_pair));
    memset(u_matched, 0, (graphSize) * sizeof(*u_matched));

    clock_t start_time = clock();
    int init_matchings = 0;

    // Stop the timer after karpSipser and before hopcroftKarp
    clock_t end_time_1 = clock();
    double time_karpSipser = (double)(end_time_1 - start_time) / CLOCKS_PER_SEC;
    //printf("Time taken by karpSipser: %f seconds\n", time_karpSipser);

    // Now, call hopcroftKarp
    start_time = clock();
    int main_matchings = parallel_pothen_fan(graph, graphSize, graphOffset,v_pair, u_matched);
    //int main_matchings = hopcroftKarp(graph, graphSize, graphColSize, pairA, pairB);

    // Stop the timer after hopcroftKarp
    clock_t end_time_2 = clock();
    double time_hopcroftKarp = (double)(end_time_2 - start_time) / CLOCKS_PER_SEC;
    //printf("Time taken by hopcroftKarp: %f seconds\n", time_hopcroftKarp);
    //printf("%d %d\n", init_matchings, main_matchings);
    return init_matchings + main_matchings;
}

size_t matching(size_t n, size_t m, xedge_t e[]) {
    static bool inited = false;
    static int* edgeCounts;
    static int* currentIndices;
    static int* graphFlat;
    static int* graphOffset;

    if (!inited) {
        inited = true;
        edgeCounts = (int*)calloc(n + 1, sizeof(int));
        currentIndices = (int*)calloc(n + 1, sizeof(int));
        graphOffset = (int*)calloc(n + 2, sizeof(int)); // n+2 for safer bounds
    }

    memset(edgeCounts, 0, (n + 1) * sizeof(int));
    memset(currentIndices, 0, (n + 1) * sizeof(int));

    for (size_t i = 0; i < m; i++) {
        edgeCounts[e[i].u]++;
        edgeCounts[e[i].v]++;
    }

    graphOffset[0] = 0;
    for (size_t i = 1; i <= n; i++) {
        graphOffset[i] = graphOffset[i - 1] + edgeCounts[i - 1];
    }
    graphOffset[n + 1] = graphOffset[n] + edgeCounts[n];

    size_t totalEdges = graphOffset[n+1];
    graphFlat = (int*)realloc(graphFlat, totalEdges * sizeof(int));
    // Populate flat graph
    for (size_t i = 0; i < m; i++) {
        int u = e[i].u;
        int v = e[i].v;
        graphFlat[graphOffset[u] + currentIndices[u]++] = v;
        graphFlat[graphOffset[v] + currentIndices[v]++] = u;
    }


    size_t maxMatching = maximumBipartiteMatching(graphFlat, n, graphOffset);

    return maxMatching;
}
