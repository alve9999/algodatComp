
#include <stdio.h>
#include <pthread.h>
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


#define MIN(a, b) (((a) <= (b)) ? (a) : (b))
#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef uint32_t u32;
typedef uint64_t u64;

typedef struct xedge_t	xedge_t;

struct xedge_t {
	int		u;	/* one of the two nodes.	*/
	int		v;	/* the other. 			*/
	int		c;	/* capacity.			*/
};


static char *input_ptr = NULL;
static char *input_start = NULL;
static size_t input_size = 0;

static void setup_mmap_input(const char *file_name)
{
    int fd = open(file_name, 0);
    //int fd = fileno(stdin);
    struct stat st;
    int err = fstat(fd, &st);
    size_t size = st.st_size;
    input_size = size;

    if (size == 0)
        fprintf(stderr, "ERROR: STDIN has size 0, meaning you didn't pipe in a file correctly!\n");
    if (err)
        fprintf(stderr, "ERRONEOUS FSTAT!!!!!\n");

    fprintf(stderr, "want to mmap file with size: %zu\n", size);
    input_ptr = mmap(NULL, size, PROT_READ, MAP_PRIVATE, fd, 0);
    input_start = input_ptr;
}

static inline i32 next_int()
{
    i32 x = 0;
    while (*input_ptr >= '0')
        x = 10 * x + (*(input_ptr++) - '0');
    ++input_ptr;
    return x;
}

static void *xmalloc(size_t s)
{
    void *p = malloc(s);
    return p;
}
static void *xcalloc(size_t n, size_t s)
{
    void *p = xmalloc(n * s);
    memset(p, 0, n * s);
    return p;
}

int preflow(int n, int m, int s, int t, xedge_t *e);

int main(void) {
    const char *paths[2] = { "../data/10k.txt", "../data/1M.txt"};
    for (int nbr=0; nbr < sizeof(paths)/sizeof(paths[0]); ++nbr) {
        setup_mmap_input(paths[nbr]);

        i32 nr_edges = 0;
        i32 edges_cap = 2;
        i32 nr_nodes = 0; // first record largest index, then increment by one later
        xedge_t *e = malloc(edges_cap * sizeof(e[0]));
        for (i32 i = 0;; ++i)
        {
            nr_edges++;
            if (nr_edges > edges_cap) {
                edges_cap *= 2;
                e = realloc(e, edges_cap * sizeof(e[0]));
            }
            i32 u = next_int();
            i32 v = next_int();
            if (u > nr_nodes)
                nr_nodes = u;
            if (v > nr_nodes)
                nr_nodes = v;
            e[i].c = 1; // all has cap 1
            e[i].u = u;
            e[i].v = v;

            if (input_ptr + 0 >= input_start + input_size) {
                break;
            }
        }
        nr_nodes += 1;
        printf("found %d edges in total!\n", nr_edges);
        printf("found %d nodes in total!\n", nr_nodes);
        int left = nr_nodes/2;
        int right = nr_nodes/2;
        int s = nr_nodes++;
        int t = nr_nodes++;
        int offset = nr_edges;
        for (int i = 0; i < left; ++i) {
            nr_edges++;
            if (nr_edges > edges_cap) {
                edges_cap *= 2;
                e = realloc(e, edges_cap * sizeof(e[0]));
            }
            e[offset].c = 1;
            e[offset].u = s;
            e[offset].v = i;
            offset += 1;
        }
        for (int i = 0; i < right; ++i) {
            nr_edges++;
            if (nr_edges > edges_cap) {
                edges_cap *= 2;
                e = realloc(e, edges_cap * sizeof(e[0]));
            }
            e[offset].c = 1;
            e[offset].u = left + i;
            e[offset].v = t;
            offset += 1;
        }

        for (int i = 0; i < 1;++i) {
            int f = preflow(nr_nodes, nr_edges, s, t, e);
            fprintf(stderr, "%s: le flow: %d\n",paths[nbr], f);
        }

    }
}
