#define _GNU_SOURCE
#include <pthread.h>
#include <sched.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <sys/types.h>

#define PRINT 0 /* enable/disable prints. */
#define KILL_ALL_PRINTS 1

#if KILL_ALL_PRINTS == 1
#define fprintf(a, ...)
#endif

#if PRINT
#define pr(...)                                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        fprintf(stderr, __VA_ARGS__);                                                                                  \
    } while (0)
#else
#define pr(...) /* no effect at all */
#endif

static char *progname;
#if PRINT
#endif
void error(const char *fmt, ...)
{
    va_list ap;
    char buf[BUFSIZ];
    va_start(ap, fmt);
    vsprintf(buf, fmt, ap);
    if (progname != NULL)
    {
        fprintf(stderr, "%s: ", progname);
    }
    fprintf(stderr, "error: %s\n", buf);
    exit(1);
}
#define MIN(a, b) (((a) <= (b)) ? (a) : (b))
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef uint32_t u32;
typedef uint64_t u64;
typedef _Atomic i16 atomic_i16;
typedef _Atomic i32 atomic_i32;
typedef _Atomic i64 atomic_i64;
typedef _Atomic u64 atomic_u64;
#define Relaxed memory_order_relaxed
#define Acquire memory_order_acquire
#define Release memory_order_release
#define AcqRel memory_order_acq_rel
#define SeqCst memory_order_seq_cst

typedef i32 NodeIdx;
typedef i32 EdgeIdx;

typedef struct
{
    int u;
    int v;
    int c;
} xedge_t;

typedef struct
{
    NodeIdx other;
    EdgeIdx ei;
} Edge;

// how many adj lists a node has
#define ADJ_LISTS 9

typedef struct
{
    i32 max_cores;
    pthread_barrier_t round_barrier;
    pthread_barrier_t concat_barrier;
    pthread_t *the_threads;
    NodeIdx s;
    NodeIdx t;
    i32 nr_nodes; /* nodes.			*/
    i32 nr_edges; /* edges.			*/
    i32 node_cap; // previous max, use to know if we should allocate more
    i32 edge_cap;

    i32 *worklist;

    // following index with node idx or edge idx
    i32 *node_e;
    i32 *node_h;

    i32 *edges_flow;

    // ONLY for THREAD 0
    i32 *thread0_node_array;
    bool *thread0_nodes_visited;

    Edge **edges_edge[ADJ_LISTS]; // each `nr_nodes` we get new node with ei, len=nodes_edges_len

    i32 *nodes_edges_len[ADJ_LISTS]; // for node, how many edges it has

    i32 *nodes_edges_cap[ADJ_LISTS]; // for node, capacity of vec

    // has 'nr_nodes' length, index with height to get nbr of nodes with that height
    i32 *height_counts;

    xedge_t *e;

    i32 xedge_per_thread;
    i32 nodes_per_thread;

    i32 workers_working; // how many are working for this graph
} Graph;
// THE GLOBAL GRAPH!
static Graph G;

static void global_relabel_initial()
{
    NodeIdx *node_array = G.thread0_node_array;
    // is 0 at start bc xcalloc, have to set to 0 again after this
    bool *nodes_visited = G.thread0_nodes_visited;
    // populate with first element: sink
    nodes_visited[G.t] = true;
    //nodes_visited[G.s] = true;
    i32 shell_timer = 1;
    i32 distance = 1;

    NodeIdx *start = &node_array[0];
    NodeIdx *head = start;
    NodeIdx *tail = start;
    NodeIdx *end = &node_array[G.nr_nodes];

    *tail = G.t;
    if (++tail >= end)
        tail = start;

    while (head != tail)
    {
        const NodeIdx node = *head;
        if (++head >= end)
            head = start;
        const i32 len = G.nodes_edges_len[0][node];
        for (i32 i = 0; i < len; ++i)
        {
            const NodeIdx other = G.edges_edge[0][node][i].other;
            if (nodes_visited[other])
                continue;
            nodes_visited[other] = true;
            if (other != G.s)
                G.node_h[other] = distance, ++G.height_counts[distance];
            // add to node_array
            *tail = other;
            if (++tail >= end)
                tail = start;
        }
        if (--shell_timer == 0)
        {
            shell_timer = (i64)(head < tail ? ((i64)tail - (i64)head) / sizeof(NodeIdx)
                                            : ((i64)tail - (i64)start + (i64)end - (i64)head) / sizeof(NodeIdx));
            ++distance;
        }
    }
}

#define MA(name, len) name = malloc(len * sizeof(name[0]))
#define RE(name, len) name = realloc(name, len * sizeof(name[0]))
#define ZEROFY(name, len) name = memset(name, 0, len * sizeof(name[0]))
static void concat_adj_lists(const i32 id)
{
    const i32 lower = id * G.nodes_per_thread;
    const i32 higher = MIN((id + 1) * G.nodes_per_thread, G.nr_nodes);
    for (i32 adjlist = 1; adjlist < ADJ_LISTS; ++adjlist)
    {
        for (i32 node = lower; node < higher; ++node) {
            const i32 add_len = G.nodes_edges_len[adjlist][node];
            const i32 old_len = G.nodes_edges_len[0][node];
            const i32 new_len = old_len + add_len;
            G.nodes_edges_len[0][node] = new_len;
            if (new_len > G.nodes_edges_cap[0][node])
            {
                G.nodes_edges_cap[0][node] = new_len;
                RE(G.edges_edge[0][node], G.nodes_edges_cap[0][node]);
            }
            memcpy(&G.edges_edge[0][node][old_len], G.edges_edge[adjlist][node], add_len * sizeof(G.edges_edge[adjlist][node][0]));
        }
    }
}
static void xedge_parsing(const i32 id)
{
    i32 me = id;
    i32 *adj_list_len = G.nodes_edges_len[me];
    i32 *adj_list_cap = G.nodes_edges_cap[me];
    Edge **adj_list = G.edges_edge[me];
    const i32 lower = id * G.xedge_per_thread;
    const i32 higher = MIN((id + 1) * G.xedge_per_thread, G.nr_edges);
    const i32 nr_nodes = G.nr_nodes;
    for (i32 i = lower; i < higher; ++i)
    {
        const xedge_t edge = G.e[i];
        const NodeIdx u = edge.u;
        const NodeIdx v = edge.v;
        const i32 u_len = adj_list_len[u]++;
        if (u_len + 1 > adj_list_cap[u])
        {
            adj_list_cap[u] *= 2;
            RE(adj_list[u], adj_list_cap[u]);
        }
        adj_list[u][u_len] = (Edge){.ei = i, .other = v};

        const i32 v_len = adj_list_len[v]++;
        if (v_len + 1 > adj_list_cap[v])
        {
            adj_list_cap[v] *= 2;
            RE(adj_list[v], adj_list_cap[v]);
        }
        adj_list[v][v_len] = (Edge){.ei = i, .other = u};
    }
}

// returns whether node is still active
static bool relabel(const i32 node)
{
    const i32 old_height = G.node_h[node]++;
    if (old_height + 1 >= G.nr_nodes)
        return G.thread0_nodes_visited[node] = true, false;

    if (--G.height_counts[old_height] == 0)
    {
        // delete this node and all nodes with greater height, in this loop
        for (i64 i = old_height + 1; i < G.nr_nodes; ++i)
            G.height_counts[i] = 0;
        for (i64 i = 0; i < G.nr_nodes; ++i)
        {
            if (G.node_h[i] > old_height)
            {
                // delete this node
                G.node_h[i] = G.nr_nodes + 1;
                G.thread0_nodes_visited[i] = true; // so we don't add it to worklist
            }
        }
        return false;
    }
    else
    {
        ++G.height_counts[old_height + 1];
        return true;
    }
}

static void work()
{
    global_relabel_initial();
    for (i32 i = 0; i < G.nr_nodes; ++i)
        G.thread0_nodes_visited[i] = false;
    const i32 nr_nodes = G.nr_nodes;
    i32 *const work = G.worklist;
    i32 *const start = &work[0];
    i32 *const end = &work[nr_nodes];
    i32 *head = start;
    i32 *tail = start;

    G.thread0_nodes_visited[G.s] = true;
    G.thread0_nodes_visited[G.t] = true;

    const i32 source_len = G.nodes_edges_len[0][G.s];
    for (int i = 0; i < source_len; ++i)
    {
        const Edge edge = G.edges_edge[0][G.s][i];
        const NodeIdx other = edge.other;
        const i32 cap = G.e[edge.ei].c;
        G.node_e[G.s] -= cap;
        G.edges_flow[edge.ei] += cap * (G.s < other ? 1 : -1);
        G.node_e[other] += cap;
        pr("initial push to: %d, with flow: %d\n", other, cap);
        if (!G.thread0_nodes_visited[other])
        {
            G.thread0_nodes_visited[other] = true;
            *tail = other;
            ++tail; // this can't go over end
        }
    }

    while (head != tail)
    {
        pr("head: %d, tail: %d\n", head - start, tail - start);
        const i32 node = *head;
        if (++head == end)
            head = start;
        pr("TOOK: %d\n", node);
    GOTO_SAME_NODE:;
        const i32 height = G.node_h[node];
        if (height >= nr_nodes) // due to gap heuristic
            continue;
        G.thread0_nodes_visited[node] = false; // should be after the height-continue thing
        i32 excess = G.node_e[node];

        const i32 edge_len = G.nodes_edges_len[0][node];
        for (i32 j = 0; j < edge_len; ++j)
        {
            const Edge edge = G.edges_edge[0][node][j];
            const i32 other = edge.other;
            const i32 ei = edge.ei;
            const i32 other_height = G.node_h[other];
            if (height == other_height + 1)
            {
                // we are one above them
                const bool c = node < other;
                const i32 e_flow = G.edges_flow[ei];
                const i32 can_push = MIN(excess, G.e[ei].c + (c ? (-e_flow) : e_flow));
                if (can_push > 0)
                {
                    pr("Push from %d to %d, flow: %d\n", node, other, can_push);
                    G.edges_flow[ei] += c ? can_push : (-can_push);
                    excess -= can_push;
                    G.node_e[other] += can_push;
                    // add to queue for p2
                    if (!G.thread0_nodes_visited[other])
                    {
                        pr("registered!\n");
                        G.thread0_nodes_visited[other] = true;
                        *tail = other;
                        ++tail;
                        if (tail == end)
                            tail = start;
                    }
                    if (excess == 0)
                        break;
                }
            }
        }
        G.node_e[node] = excess; // update global
        // if we still have excess, we maxed all the edges and thus have to relabel
        // add to relabel queue, 0 is sentinel for relabel
        if (excess > 0)
        {
            pr("relabel %d\n", node);
            if (relabel(node))
                goto GOTO_SAME_NODE; // or add to end of worklist possibly
        }
    }
}

static i64 GRAPH_NUM = 0;

static i32 xpreflow()
{
    // start by pushing as much as possible (limited by the edge capacity) from the source to its neighbors.
    // work();
    return G.node_e[G.t];
}
static void new_graph(const i32 id)
{
    const i32 nr_nodes = G.nr_nodes;
    const i32 nr_edges = G.nr_edges;
    // instead of setting to new maxima, we could double G.max_nr_nodes instead
    if (nr_nodes > G.node_cap)
    {
        i32 prev_node_cap = G.node_cap;
        G.node_cap = nr_nodes;

        RE(G.node_e, G.node_cap);
        RE(G.node_h, G.node_cap);
        // maybe opt: if not 0, lower cap
        for (i32 adjlist = 0; adjlist < ADJ_LISTS; ++adjlist)
        {
            RE(G.nodes_edges_len[adjlist], G.node_cap);
            RE(G.edges_edge[adjlist], G.node_cap);
            RE(G.nodes_edges_cap[adjlist], G.node_cap);
            for (i32 i = prev_node_cap; i < G.node_cap; ++i)
            {
                G.nodes_edges_cap[adjlist][i] = adjlist == 0 ? 256 : 32; // the best
                MA(G.edges_edge[adjlist][i], G.nodes_edges_cap[adjlist][i]);
            }
        }
        RE(G.height_counts, G.node_cap);
        // set new ones to 128

        RE(G.thread0_node_array, G.node_cap);
        RE(G.thread0_nodes_visited, G.node_cap);
        RE(G.worklist, G.node_cap);
    }
    // set things to zero
    // isn't faster to just zero the v, memset is fastest.
    for (i32 adjlist = 0; adjlist < ADJ_LISTS; ++adjlist)
    {
        ZEROFY(G.nodes_edges_len[adjlist], nr_nodes);
    }

    { // thread/worker stuff
        fprintf(stderr, "nr_nodes: %d, nr_edges: %d\n", G.nr_nodes, G.nr_edges);
    }

    { // set how much each thread shall do xedges adj creation

        // G.xedge_nodes_per_thread = ((double)4000 / (double)G.nr_edges);
        //  much faster for some reason. Goes from 7700 with above to 8300 with below
        // G.xedge_per_thread = ((40000.0) / (double)sqrt(G.nr_edges));
        //  minimum work for a thread
        G.xedge_per_thread = 500; // best

        // make sure EVERY edge get's a thread.
        G.xedge_per_thread = MAX(G.xedge_per_thread, (G.nr_edges - 1) / ADJ_LISTS + 1); // round up
        fprintf(stderr, "xedge_per_thread: %d\n", G.xedge_per_thread);
        // how many will be working now.
        G.workers_working = MIN((nr_edges - 1) / G.xedge_per_thread + 1, G.max_cores);

        G.nodes_per_thread = 100; // TODO set right
        G.nodes_per_thread = MAX(G.nodes_per_thread, (G.nr_nodes - 1) / (G.max_cores) + 1);
    }
    // xedge_t parsing begins here !!!!!!!!!!!!!!! XEDGE
    pthread_barrier_wait(&G.round_barrier);
    if (nr_edges > G.edge_cap)
    {
        if (G.edge_cap == 0)
            G.edge_cap = 2048;
        while (G.edge_cap < nr_edges)
            G.edge_cap *= 2;
        RE(G.edges_flow, G.edge_cap);
    }
    ZEROFY(G.edges_flow, nr_edges);
    ZEROFY(G.thread0_nodes_visited, nr_nodes);
    ZEROFY(G.node_e, nr_nodes);
    ZEROFY(G.height_counts, nr_nodes);

    // second phase of xedge_parsing, TODO optimize away later
    pthread_barrier_wait(&G.round_barrier); // XEDGE END !!!!!!!!!!!!!

    { // optimization based on switching ends if initial capacitites are smaller at sink
        i64 source_caps = 0;
        const i32 source_len = G.nodes_edges_len[0][G.s];
        for (i32 i = 0; i < source_len; ++i)
        {
            EdgeIdx ei = G.edges_edge[0][G.s][i].ei;
            source_caps += (i64)G.e[ei].c;
        }
        i64 sink_caps = 0;
        const i32 sink_len = G.nodes_edges_len[0][G.t];
        for (i32 i = 0; i < sink_len; ++i)
        {
            EdgeIdx ei = G.edges_edge[0][G.t][i].ei;
            sink_caps += (i64)G.e[ei].c;
        }
        if (sink_caps < source_caps)
        {
            i32 tmp = G.s;
            G.s = G.t;
            G.t = tmp;
        }
    }
    G.node_h[G.s] = G.nr_nodes;
    G.node_h[G.t] = 0;
}

static void *green_room(void *void_data)
{
    const i32 id = (i32)(i64)void_data;
    for (;;)
    {
        // // wait for structuring of xedge_t starts
        pthread_barrier_wait(&G.round_barrier);
        const bool special = (GRAPH_NUM % G.max_cores == (id));

        if (special)
            new_graph(id);
        else
        {
            pthread_barrier_wait(&G.round_barrier);
            xedge_parsing(id);
            pthread_barrier_wait(&G.concat_barrier);
            concat_adj_lists(id);
            pthread_barrier_wait(&G.round_barrier);
        }

        pthread_barrier_wait(&G.round_barrier);
        if (special)
            work();
        pthread_barrier_wait(&G.round_barrier);
    }
    return NULL;
}
void first_init()
{
    { // main thread (THIS ONE) at zero (0)
        cpu_set_t main_cpuset;
        CPU_ZERO(&main_cpuset);
        CPU_SET(0, &main_cpuset); // Set the main thread to t0
        pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &main_cpuset);
    }

    G.max_cores = ADJ_LISTS; // MAX CORES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // max_cores + 1, so main thread can start the work
    pthread_barrier_init(&G.round_barrier, NULL, G.max_cores + 1);
    pthread_barrier_init(&G.concat_barrier, NULL, G.max_cores);
    G.the_threads = malloc(G.max_cores * sizeof(pthread_t));

    for (i64 i = 0; i < G.max_cores; ++i)
    {
        pthread_create(&G.the_threads[i], NULL, green_room, (void *)i);

        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        i32 idx = 4 + 4 * i;
        CPU_SET(idx, &cpuset);

        pthread_setaffinity_np(G.the_threads[i], sizeof(cpu_set_t), &cpuset);
    }
}

#ifndef MAIN
static bool FIRST_TIME = true;
static i64 SOME_COUNT = 0;
// CALLED BY FORSETE, NOT STATIC
int preflow(int n, int m, int s, int t, xedge_t *e)
{
    G.nr_nodes = n;
    G.nr_edges = m;
    G.s = s;
    G.t = t;
    G.e = e;
    if (FIRST_TIME)
        first_init(), FIRST_TIME = false;
    pthread_barrier_wait(&G.round_barrier);
    i32 spec_id = GRAPH_NUM % G.max_cores;

    pthread_barrier_wait(&G.round_barrier);
    xedge_parsing(spec_id);
    pthread_barrier_wait(&G.concat_barrier);
    concat_adj_lists(spec_id);
    pthread_barrier_wait(&G.round_barrier);

    pthread_barrier_wait(&G.round_barrier);
    pthread_barrier_wait(&G.round_barrier);

    i32 f = xpreflow();
    SOME_COUNT += m;
    if (SOME_COUNT > 100000000) {
        GRAPH_NUM++;
        SOME_COUNT = 0;
    }
    return f;
}
#else
int main(int argc, char *argv[])
{
    (void)argc;
    i32 f; /* output from preflow.		*/

    progname = argv[0]; /* name is a string in argv[0]. */

    // get from STDIN
    setup_mmap_input();
    i32 nr_nodes = next_int();
    i32 nr_edges = next_int();

    /* skip C and P from the 6railwayplanning lab in EDAF05 */
    next_int();
    next_int();
    i32 s = 0;
    i32 t = nr_nodes - 1;

    first_init(nr_nodes, nr_edges);

    new_graph(nr_nodes, nr_edges, s, t, NULL);

    FILE *in;   /* input file set to stdin	*/
    in = stdin; /* same as System.in in Java.	*/
    fclose(in);

    f = xpreflow();

    printf("f = %d\n", f);
    return 0;
}
#endif
