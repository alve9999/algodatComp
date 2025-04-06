#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <omp.h>

typedef struct {
    size_t        u;    /* odd.        */
    size_t        v;    /* even.    */
} xedge_t;

size_t matching(size_t n, size_t m, xedge_t e[]);

int main() {
    int u, v, capacityEdges = 100, m = 0;
    xedge_t *edges = malloc(capacityEdges * sizeof(*edges));

    int n = 0;

    while (scanf("%d %d", &u, &v) == 2) {
        if (m == capacityEdges) {
            capacityEdges *= 2;
            edges = realloc(edges, capacityEdges * sizeof(*edges));
        }
        u += 1; // just increase both by 1, to allow 0 as sentinel, cuz Skeppstedts input is like that. Does NOT hurt anything
        v += 1;

        edges[m].u = u;
        edges[m].v = v;
        m++;
        if (u > n) n = u;
        if (v > n) n = v;
    }

    double start = omp_get_wtime();

    size_t res = matching(n, m, edges);
    double end = omp_get_wtime();
    printf(" << GOT: %zd\n    took: %f seconds >>\n", res, end - start);
}
