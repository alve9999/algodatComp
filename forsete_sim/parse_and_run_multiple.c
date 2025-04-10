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

    size_t u, v, capacityEdges = 100, m = 0;
    size_t n = 0;
    xedge_t *edges = malloc(capacityEdges * sizeof(*edges));

    double total_time = 0;

    int i = 1;
    for (; ; ++i) {
        if (total_time >= 120) {
            printf(" >> TIME OVER! 120 s HAS PASSED! YOU CLEARED %d GRAPHS\n", i-1);
            break;
        }
        char filename[256];

        snprintf(filename, sizeof(filename), "data/%d_1M.txt", i);
        printf(" >> Opening file: %s   :: ", filename);

        FILE *fp = fopen(filename, "r");
        if (fp == NULL) {
            printf(" >> No more files! Stopped at %i\n", i);
            break;
        }

        m = 0;

        while (fscanf(fp, "%zd %zd", &u, &v) == 2) {
            if (m >= capacityEdges) {
                capacityEdges *= 2;
                edges = realloc(edges, capacityEdges * sizeof(*edges));
            }
            //u += 1; // just increase both by 1, to allow 0 as sentinel, cuz Skeppstedts input is like that. Does NOT hurt anything
            //v += 1;

            edges[m].u = u;
            edges[m].v = v;
            m++;
            if (u > n) n = u;
            if (v > n) n = v;
        }
        fclose(fp);

        double start = omp_get_wtime();

        size_t res = matching(n, m, edges);
        double end = omp_get_wtime();
        printf(" << GOT: %zd    took: %f seconds >>\n", res, end - start);
        total_time += end - start;
    }
    printf("TOTAL TIME TAKEN: %lf\n", total_time);
}
