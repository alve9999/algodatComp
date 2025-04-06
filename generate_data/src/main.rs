use std::hash::RandomState;
use std::hash::{BuildHasher, Hasher};

// WARNING: WILL GENERATE DUPLICATE EDGES, LOOK AT COMMENT BELOW
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 1000;
    if args.len() == 1 {
        eprintln!("No number of nodes, default is 1000");
    } else {
        n = args[1].parse().unwrap();
    }
    let (_nr_nodes, edges) = generate(n);
    print_edges(&edges);
}

// takes number of nodes (even)
// what this will do: left and right side will have n nodes
// from [1 ... n]. left (u) will be odd [1, 3, 5 ... n-1], right (v) will be even [2, 4, 6 ... n]
// it generates n * 4 edges between these two sets. Could contain duplicates, see comment below
fn generate(n: u64) -> (u64, Vec<(u64, u64)>) {
    assert!(n % 2 == 0);
    let total_edges = n * 4;
    let mut edges = Vec::with_capacity(total_edges as usize);
    let random_state = RandomState::new();
    let mut hasher = random_state.build_hasher();
    let mut counter = 0u64;
    for _ in 0..total_edges {
        // choose u and v for this edge!
        // u must be odd, v must be even
        hasher.write_u64(counter);
        counter += 1;
        let ra = hasher.finish();
        // if n = 1000:
        // (num: [0, 499]) -> (even num: [0, 998]) -> (odd num: [1, 999])
        let u = (ra % (n / 2)) * 2 + 1;
        hasher.write_u64(counter);
        counter += 1;
        let ra = hasher.finish();
        // (num: [0, 499]) -> (num: [1, 500]) -> (even num: [2, 1000])
        let v = ((ra % (n / 2)) + 1) * 2;
        edges.push((u, v));
    }

    // THIS REMOVES DUPLICATES, BUT IT IS O(N^2), AND EVERYTHING SEEMS TO WORK EVEN WITH DUPLICATES, SO....
    if false {
        let mut i = 0;
        while i + 1 < edges.len() {
            let mut j = i + 1;
            while j < edges.len() {
                if edges[i].0 == edges[j].0 && edges[i].1 == edges[j].1 {
                    // remove dup
                    edges.remove(j);
                } else {
                    j += 1;
                }
            }
            i += 1;
        }
    }
    (n, edges)
}

fn print_edges(edges: &[(u64, u64)]) {
    for &edge in edges.iter() {
        println!("{} {}", edge.0, edge.1);
    }
}
