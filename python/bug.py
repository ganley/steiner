import collections
import copy
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import sys

import steinlib     # local


if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        G = steinlib.parse(f)

    apsp = dict(nx.algorithms.shortest_paths.all_pairs_dijkstra_path_length(G))

    print(f"{apsp[22][24]=}")

    path = nx.algorithms.shortest_paths.shortest_path(G, 22, 24)
    print(f"{path=}")
    for a,b in zip(path, path[1:]):
        print(f"{a=} {b=} {G[a][b]['weight']=}")

    print("----")
    print(f"{apsp[24][22]=}")

    path = nx.algorithms.shortest_paths.shortest_path(G, 24, 22)
    print(f"{path=}")
    for a,b in zip(path, path[1:]):
        print(f"{a=} {b=} {G[a][b]['weight']=}")

    print("----")
    print(f"{G[24][28]['weight']=}")
    print(f"{G[28][18]['weight']=}")
    print(f"{G[18][43]['weight']=}")
    print(f"{G[43][22]['weight']=}")
