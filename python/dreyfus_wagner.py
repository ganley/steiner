import collections
import copy
import itertools
import matplotlib.pyplot as plt
import networkx as nx
import sys

import steinlib     # local


def _powerset(s):
    return (c for r in range(1, len(s)) for c in itertools.combinations(s, r))


# memo has key (u, list of terminals) and value (v, subset, weight)
# TODO change those to named tuples?
def _dw_solve(G: nx.Graph, terminals, u, apsp, memo={}):
    #pmrint(f"_dw_solve {u=} {terminals=}")
    if (u, terminals) in memo:
        # print(f"{memo[(u,terminals)]=}")
        return memo[(u, terminals)][2]

    if len(terminals) >= len(G.graph["terminals"]) - 2:
        print(f"{u=} {terminals=}")

    maxweight = G.size(weight="weight")

    if len(terminals) == 1:
        bestv = None
        bestsubset = None
        bestweight = apsp[u][terminals[0]]
        #print(f"degenerate: {weight=}")
    else:
        assert len(terminals) > 1
        bestweight = maxweight
        for v in G.nodes:
            # print(f"{v=}")
            xvert_set = set(terminals)
            # print(f"{_powerset(terminals)=}")
            for x in _powerset(terminals):
                xset = set(x)
                xpset = xvert_set.difference(xset)
                w1 = apsp[u][v]
                if w1 >= bestweight:
                    continue
                w2 = _dw_solve(G, tuple(sorted(list(xset))), v, apsp, memo)
                if w1 + w2 >= bestweight:
                    continue
                w3 = _dw_solve(G, tuple(sorted(list(xpset))), v, apsp, memo)
                w = w1 + w2 + w3
                #print(f"{u=} {v=} {xset=} {xpset=} {w1=} {w2=} {w3=} {w=}")
                if w < bestweight:
                    bestv = v
                    bestsubset = copy.copy(x)
                    bestweight = w

    assert bestweight < maxweight
    memo[(u, terminals)] = (bestv, bestsubset, bestweight)
    return bestweight


def _dw_build_graph(T, terminals, u, G, memo):
    if len(terminals) == 1:
        print(f"{terminals[0]=} {u=}")  # {G[terminals[0]][u]['weight']=}")
        if terminals[0] != u:
            T.add_edge(terminals[0], u, weight=G[terminals[0]][u]["weight"])
        return
    bestv, bestsubset, _ = memo[(u, terminals)]
    _dw_build_graph(T, bestsubset, bestv, G, memo)
    term_set = set(terminals)
    xpset = term_set.difference(bestsubset)
    _dw_build_graph(T, tuple(sorted(list(xpset))), bestv, G, memo)


# This is the only external function. It takes a graph and returns:
# - A NetworkX graph containing the optimal Steiner tree if weight_only = False
# - The total weight of the optimal Steiner tree if weight_only = True
def dreyfus_wagner(G: nx.Graph, weight_only=False):
    apsp = dict(nx.algorithms.shortest_paths.all_pairs_dijkstra_path_length(G))
    cache = {}
    terminals = sorted(G.graph["terminals"])
    u0 = terminals[0]
    t0 = tuple(terminals[1:])
    weight = _dw_solve(G, t0, u0, apsp, cache)

    if weight_only:
        return weight

    T = nx.Graph(Remark="Optimal Steiner tree via Dreyfus-Wagner")
    _dw_build_graph(T, tuple(terminals[1:]), terminals[0], G, cache)

    bestv, _, _ = cache[(u0, t0)]
    print(f"{u0=} {bestv=}")
    T.add_edge(u0, bestv, weight=G[u0][bestv]["weight"])

    return T


if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        G = steinlib.parse(f)

    print(G)
    print(G.graph)
    print(G.nodes)
    print(G.nodes.data())
    print(G.edges)

    T = dreyfus_wagner(G)
    print(T.size(weight="weight"))

    nx.draw(T)
    plt.savefig("steinlib.png")
