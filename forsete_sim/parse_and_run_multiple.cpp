#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>

#include <unordered_set>

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif
typedef struct {
    size_t        u;    /* odd.        */
    size_t        v;    /* even.    */
} xedge_t;
size_t matching(size_t n, size_t m, xedge_t e[]);
#ifdef __cplusplus /* If this is a C++ compiler, end C linkage */
}
#endif



static unsigned int generate_random(){
    return rand();
}

bool operator==(const xedge_t& lhs, const xedge_t& rhs) {
    return lhs.u == rhs.u && lhs.v == rhs.v;
}
struct xedge_hash {
    size_t operator()(const xedge_t& edge) const {
        // A simple and effective hash combination
        return std::hash<size_t>()(edge.u) ^ (std::hash<size_t>()(edge.v) << 1);
    }
};

#include <vector>

// create graph from 1..1M  where odd are u and even are v.
static void create_graph(xedge_t *edges, int n, int m) {

    size_t total_pairs = 500000UL * 500000UL;
    std::vector<bool> vec(total_pairs, false);

    int pre_m = 0;
    while (pre_m < m) {
        size_t u = (generate_random() % (n/2)) * 2 + 1;
        size_t v = ((generate_random() % (n/2)) + 1) * 2;
        xedge_t e = {u, v};
        size_t index = ((u -1) / 2) * 500000UL + (v-2)/2;
        if (!vec.at(index)) {
            vec[index] = true;
            edges[pre_m++] = e;
        }
    }
}

int main() {

    size_t m = 4000000;
    size_t n = 1000000;
    xedge_t *edges = (xedge_t *)malloc(m * sizeof(*edges));

    double total_time = 0;

    int i = 1;
    for (; ; ++i) {
        if (total_time >= 120) {
            printf(" >> TIME OVER! 120 s HAS PASSED! YOU CLEARED %d GRAPHS\n", i-1);
            break;
        }

        double start = omp_get_wtime();
        create_graph(edges, n, m);

        size_t res = matching(n, m, edges);
        double end = omp_get_wtime();
        total_time += end - start;
        // total_time += 1.078; // WASTE TIME
        printf(" << %d: GOT: %zd    took: %f sec. Total time: %f sec. >>\n", i, res, end - start, total_time);
    }
    printf("TOTAL TIME TAKEN: %lf\n", total_time);
}
