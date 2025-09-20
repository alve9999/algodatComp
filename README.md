# Code for 2025 Algorithms and data structures LTH course competition. 
We won. 

### What we do? 

We do maximum bipartite graph matching. We use the PPF+ algorithm (parallel pothen fan+). It's really fast for our graphs (1 million nodes, 4 million edges). It's super multi-threaded, both graph creation and the actual algorithm. It stores the adjacency lists contiguously in memory. We use OpenMP for threading. We do manual recursion to avoid stack overflow (it's also faster). Fairness works really great. 'grok.c' is also okay solution. 

### Structure: 

- aitesting: here are all the good solutions. It started out as just the ai generated bad solutions, but now it keeps our premier winner solution: 'alve_experimental.c'. Do 'make' to run some tests there if you want. But I prefer ./forsete_sim for testing. 

- forsete_sim: the greatest folder. Do 'run.sh myfile.c' to test your solution. This simulates forsete very well. Because it accurately creates the graphs live just like forsete, which makes the runtime much more deterministic than if we just sleep some time. Some C++ code to create the graph efficiently without duplicates.

- preflow_sol: nothing of value. Does not work as it's for undirected graphs. It's also slower so...

- generate_data: Old rust code for generating data files. run 'make' to create data, which will be copied to a ./data folder. 
