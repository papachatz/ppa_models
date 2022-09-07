import numpy as np

#paths given by DFS algorithm
#indexes of prefix nodes are stored
paths = [];
paths.append(["11:10","11:8","11:0"]);
paths.append(["9:8","11:8","11:0"]);
paths.append(["7:6","7:4","7:0","11:0"]);
paths.append(["5:4","7:4","7:0","11:0"]);
paths.append(["3:2","3:0","7:0","11:0"]);
paths.append(["1:0","3:0","7:0","11:0"]);

# G matrix construction
G = [];     # prefix node matrix
Gd = [];    # depth prefix node matrix
for path in paths:
    for node in reversed(path):
        if(node not in G):
            G.append(node);
            Gd.append(path.index(node));
        if((node in G) and Gd[G.index(node)]<path.index(node)):
            Gd[G.index(node)] = path.index(node);

# Sorting based on Gd
G=zip(G, Gd);
print G;
G = sorted(G, key=lambda tup: tup[1]);
G, Gd = zip(*G);
G = list(G);
print G;

# Initialize A matrix
A = [[0 for j in range(len(G))] for i in range(len(paths))];

# A matrix construction
for path in paths:
    for node in G:
        if(node in path):
            A[paths.index(path)][G.index(node)] = 1;

#print A matrix
for row in A:
    print row;
