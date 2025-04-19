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
static inline bool dfs_la_ts(const int* graph, const int* graphOffset,int start, visit_t* visited,int* v_pair, u_i_t *stack, bool order){
    int stack_at = 0;
    int u = start;
    WORK_ON_NEW_U:;
    int i = order ? graphOffset[u] : (graphOffset[u+1] - 1);
    // take from stack
    assert(u != 0); // should not be 0
    for (int j = graphOffset[u]; j < graphOffset[u+1]; j++) {
        int v = graph[j];
        if(v_pair[v] == 0){ // unmatched v. Can match them together
            if (!atomic_load_explicit(&visited[v].b, memory_order_relaxed)) {
                if (atomic_exchange_explicit(&visited[v].b,true, memory_order_relaxed) == 0){
                    // found unmatched v, that wasn't visited. Let's match with this one.
                    v_pair[v] = u;
                    // go backwards through stack and match things
                    while (stack_at > 0) {
                        stack_at--;
                        u = stack[stack_at].u;
                        i = stack[stack_at].i;
                        int v = graph[i];
                        assert(u != 0);
                        v_pair[v] = u;
                    }
                    return 1;
                }
            }
        }
    }
    CONTINUE_ON_SECOND_LOOP:;
    for (; order ? (i < graphOffset[u+1]) : (i >= graphOffset[u]); order ? i++ : i--) {
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
        return 0; // nothing left
    }
    stack_at--;
    u = stack[stack_at].u;
    i = stack[stack_at].i + (order ? (1) : -1); // after ourselves
    goto CONTINUE_ON_SECOND_LOOP;
}

// pothen-fan threads
#define DFS_NR_THREADS 32
// graph-creation threads. These need to be fewer, want one per real core on a NUMA node
#define CREATION_NR_THREADS 4

static size_t parallel_pothen_fan(const int* graph, int graphSize, const int* graphOffset, int* v_pair, bool *u_matched){
    static u_i_t *stacks;
    static visit_t* visited;
    static int* worklist;
    if (stacks == NULL) {
        stacks = calloc(graphSize * DFS_NR_THREADS, sizeof(*stacks));
        visited = calloc(graphSize, sizeof(*visited));
        worklist = malloc(graphSize * sizeof(int));
    }

    // allocate stacks
    bool path_found = 1;
    size_t matchings = 0;
    bool order = false;
    // FOR V NODES
    while(path_found){
        memset(visited, 0, sizeof(*visited)*graphSize);
        path_found = 0;
        order = !order;
        size_t worklist_size = 0;
        for (int i = 0; i < (graphSize >> 1); ++i) {
            if (!u_matched[i]) {
                worklist[worklist_size++] = 1 + (i << 1);
            }
        }

        #pragma omp parallel for schedule(guided) reduction(+:matchings) reduction(||:path_found)   num_threads(DFS_NR_THREADS)
        for(size_t idx = 0; idx < worklist_size; idx+=1){ // UNMATCHED u VERTICES
            int i = worklist[idx];
            // SAFETY: this read is okay. Because no unmatched u's will be visited from searches starting on other unmatches u's. only matched u's will be found.
            int t_id = omp_get_thread_num();
            u_i_t *stack = &stacks[graphSize * t_id];
            // these are inlined, so get two monomorphised versions
            bool res = false;
            if (order) res = dfs_la_ts(graph,graphOffset,i,visited,v_pair,  stack, true);
            else       res = dfs_la_ts(graph,graphOffset,i,visited,v_pair,  stack, false);
            if(res) {
                u_matched[i >> 1] = true;
                path_found = true;
                matchings++;
            }
        }
    }
    return matchings;
}

static size_t maximumBipartiteMatching(const int* graph, int graphSize, int* graphOffset) {
    static int *v_pair;
    static bool *u_matched;
    if (v_pair == NULL) {
        v_pair = malloc((graphSize) * sizeof(*v_pair));
        u_matched = malloc((graphSize) * sizeof(*u_matched));
    }
    memset(v_pair, 0, (graphSize) * sizeof(*v_pair));
    memset(u_matched, 0, (graphSize) * sizeof(*u_matched));

    // Now, call hopcroftKarp
    return parallel_pothen_fan(graph, graphSize, graphOffset,v_pair, u_matched);
}

size_t matching(size_t n, size_t m, xedge_t e[]) {
    static int* edgeCounts;
    static int* graphFlat;
    static int* graphOffset;
    static int *t_offsets[CREATION_NR_THREADS];

    n += 1; // n goes from 1M to  (1M + 1), as node indexing is [1,1M]

    if (graphFlat == NULL) {
        edgeCounts = (int*)calloc(n , sizeof(int));
        graphOffset = (int*)calloc(n, sizeof(int));
        graphFlat = (int*)malloc((4*m) * sizeof(int));
        for (int i = 0; i < CREATION_NR_THREADS; ++i) {
            t_offsets[i] = calloc(n, sizeof(*t_offsets));
        }
    }
    memset(edgeCounts, 0, n  * sizeof(int));
    memset(graphOffset, 0, n*sizeof(int));

    #pragma omp parallel for schedule(static) num_threads(CREATION_NR_THREADS)
    for (int i = 0; i < CREATION_NR_THREADS; ++i)
        memset(t_offsets[i], 0, (n * sizeof(*t_offsets)));

    // each thread has offsets for each node

    // first count, for each thread, how many nodes each of them have
    // (disjoint edges, can have same nodes)
    // basically local degrees for nodes
    #pragma omp parallel for schedule(static) num_threads(CREATION_NR_THREADS)
    for (size_t i = 0; i < m; ++i) {
        const int tid = omp_get_thread_num();
        t_offsets[tid][e[i].u]++;
        t_offsets[tid][e[i].v]++;
    }
    // convert to thread local offsets and global degrees
    #pragma omp parallel for schedule(static) num_threads(CREATION_NR_THREADS)
    for (size_t i = 1; i < n; ++i) {
        for (int j = 0; j < CREATION_NR_THREADS; ++j)
            graphOffset[i+1] += t_offsets[j][i];
        for (size_t j = 1; j < CREATION_NR_THREADS; ++j)
            t_offsets[j][i] += t_offsets[j - 1][i];
    }
    // convert to global offsets
    #pragma omp single
    {
        graphOffset[0] = 0;
        for (size_t i = 1; i < n; ++i)
            graphOffset[i] += graphOffset[i - 1];
    }
    // fill the edges
    #pragma omp parallel for schedule(static) num_threads(CREATION_NR_THREADS)
    for (size_t i = 0 ; i < m; ++i) {
        const int tid = omp_get_thread_num();
        xedge_t ed = e[i];
        graphFlat[graphOffset[ed.u] + --t_offsets[tid][ed.u]] = ed.v;
        graphFlat[graphOffset[ed.v] + --t_offsets[tid][ed.v]] = ed.u;
    }

    return maximumBipartiteMatching(graphFlat, n, graphOffset);
}
