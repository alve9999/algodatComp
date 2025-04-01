use std::hash::RandomState;
use std::hash::{BuildHasher, Hasher};

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

// takes number of nodes on left and right side, plus number of edges
fn generate(n: u64) -> (u64, Vec<(u64, u64)>) {
    let mut edges = Vec::new();
    let random_state = RandomState::new();
    let mut hasher = random_state.build_hasher();
    let mut counter = 0u64;
    for i in 0..n {
        // for each node, generate some edges
        hasher.write_u64(counter);
        counter += 1;
        let ra = hasher.finish();
        let edges_to_add = (ra % 10).min(n - 1) + 1;
        let mut went_to_already = Vec::new();
        for _ in 0..edges_to_add {
            // choose node on other side to draw edge to
            loop {
                hasher.write_u64(counter);
                counter += 1;
                let ra = hasher.finish();
                let edge_to = (ra % n) + n; // other side, has +n
                if went_to_already.contains(&edge_to) {
                    continue;
                }
                went_to_already.push(edge_to);
                // also add the edge
                edges.push((i, edge_to));
                break;
            }
        }
        // added nr of edges we wanted from n to somebody else
    }
    (n * 2, edges)
}

fn print_edges(edges: &[(u64, u64)]) {
    for &edge in edges.iter() {
        println!("{} {}", edge.0, edge.1);
    }
}
