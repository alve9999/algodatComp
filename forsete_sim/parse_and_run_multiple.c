#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <ctype.h>

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


        int fd = open(filename, O_RDONLY);
        if (fd < 0) {
            printf(" >> No more files! Stopped.\n");
            break;
        }
        // Get the size of the file.
        struct stat st;
        fstat(fd, &st);
        size_t filesize = st.st_size;

        // Memory-map the file.
        char *data = mmap(NULL, filesize, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);  // Close FD after mapping; it is no longer needed.
        if (data == MAP_FAILED) {
            perror("mmap failed");
            return 1;
        }

        m = 0; n = 0;
        // Pointer for iterating through file contents.
        char *p = data;
        char *end_file = data + filesize;
        while (p < end_file) {
            // Skip any whitespace.
            while (p < end_file && isspace((unsigned char)*p)) {
                p++;
            }
            if (p >= end_file) break;

            // Parse first number (u).
            char *next;
            size_t u = strtoul(p, &next, 10);
            p = next;

            // Skip whitespace between numbers.
            while (p < end_file && isspace((unsigned char)*p)) {
                p++;
            }

            // Parse second number (v).
            size_t v = strtoul(p, &next, 10);
            p = next;

            // Resize array if necessary.
            if (m >= capacityEdges) {
                capacityEdges *= 2;
                edges = realloc(edges, capacityEdges * sizeof(*edges));
            }

            edges[m].u = u;
            edges[m].v = v;
            m++;
            if (u > n) n = u;
            if (v > n) n = v;
        }

        // Clean up.
        munmap(data, filesize);

        double start = omp_get_wtime();

        size_t res = matching(n, m, edges);
        double end = omp_get_wtime();
        total_time += end - start;
        total_time += 1.078; // WASTE TIME
        printf(" << GOT: %zd    took: %f sec. Total time: %f sec. >>\n", res, end - start, total_time);
    }
    printf("TOTAL TIME TAKEN: %lf\n", total_time);
}
