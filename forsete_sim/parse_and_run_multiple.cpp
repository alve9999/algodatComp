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

#include <vector>
#include <algorithm>

struct xedge_hash {
    size_t operator()(const xedge_t& edge) const {
        // A simple and effective hash combination
        return std::hash<size_t>()(edge.u) ^ (std::hash<size_t>()(edge.v) << 1);
    }
};

// Define the equality operator for xedge_t to support std::unique.
bool operator==(const xedge_t& lhs, const xedge_t& rhs) {
    return lhs.u == rhs.u && lhs.v == rhs.v;
}

// Helper: comparison lambda for sorting xedge_t objects.
static bool less_xedge(const xedge_t &a, const xedge_t &b) {
    return (a.u < b.u) || ((a.u == b.u) && (a.v < b.v));
}

// Simple random generator function.
static unsigned int generate_random(){
    return rand();
}

// Create graph with edges having odd u and even v in the range 1..n.
// We aim to generate m unique edges and store them in the provided array.
static void create_graph(xedge_t *edges, int n, int m) {
    // We'll use a vector to accumulate candidate edges.
    std::vector<xedge_t> edgeVec;
    // Reserve the initial target count.
    edgeVec.reserve(m);

    // First pass: generate m candidate edges.
    for (int i = 0; i < m; i++) {
        size_t u = (generate_random() % (n/2)) * 2 + 1;         // odd in 1..n
        size_t v = ((generate_random() % (n/2)) + 1) * 2;         // even in 1..n
        edgeVec.push_back({u, v});
    }

    // Sort the candidate edges.
    std::sort(edgeVec.begin(), edgeVec.end(), less_xedge);
    // Remove duplicates.
    auto unique_end = std::unique(edgeVec.begin(), edgeVec.end());
    edgeVec.erase(unique_end, edgeVec.end());

    // If we have fewer than m unique edges, generate additional ones.
    while (static_cast<int>(edgeVec.size()) < m) {
        int needed = m - static_cast<int>(edgeVec.size());
        std::vector<xedge_t> extra;
        extra.reserve(needed);
        for (int i = 0; i < needed; i++){
            size_t u = (generate_random() % (n/2)) * 2 + 1;
            size_t v = ((generate_random() % (n/2)) + 1) * 2;
            extra.push_back({u, v});
        }
        // Sort and remove duplicates among the new edges.
        std::sort(extra.begin(), extra.end(), less_xedge);
        auto extra_unique_end = std::unique(extra.begin(), extra.end());
        extra.erase(extra_unique_end, extra.end());

        // Merge the extra edges into our main vector.
        edgeVec.insert(edgeVec.end(), extra.begin(), extra.end());
        // Sort and deduplicate the merged vector.
        std::sort(edgeVec.begin(), edgeVec.end(), less_xedge);
        auto merged_unique_end = std::unique(edgeVec.begin(), edgeVec.end());
        edgeVec.erase(merged_unique_end, edgeVec.end());
    }

    // Finally, copy the first m unique edges into the output array.
    for (int i = 0; i < m; i++){
         edges[i] = edgeVec[i];
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
            printf(" >> TIME OVER! 120 s HAS PASSED! YOU CLEARED %d GRAPHS\n", i-2);
            break;
        }

        double start = omp_get_wtime();
	srand(i);
        create_graph(edges, n, m);

        size_t res = matching(n, m, edges);
        double end = omp_get_wtime();
        total_time += end - start;
        total_time += 0.046; // WASTE TIME
        printf(" << %d: GOT: %zd    took: %f sec. Total time: %f sec. >>\n", i, res, end - start, total_time);
    }
    printf("TOTAL TIME TAKEN: %lf\n", total_time);
}
