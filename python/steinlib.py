# steinlib.py
#
# Read a Steiner tree test instance in STP (SteinLib) format into a NetworkX
# graph. It does not handle directed graphs, obstacles, rooted or degree
# constrained instances.

import matplotlib.pyplot as plt
import networkx as nx
import re
import sys


# This is largely untested, as all of my testcases are valid STP
class ParseError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


# This is the main function - it takes a file object and returns a NetworkX
# Graph
def parse(f):
    section = None
    for rawline in f:
        if line := rawline.strip():
            first, *rest = line.upper().split()

            if first == "EOF":
                return G

            if first in ["OBSTACLES", "ARCS", "A", "ROOTP", "ROOT", "TP", "MD"]:
                raise ParseError(line, f"Unsupported feature {first}")

            if first == "33D32945":    # STP magic number
                G = nx.Graph(STP_version=rest[-1])
                continue

            if first == "SECTION":
                section = rest[0]
                if section == "MAXIMUMDEGREES":
                    raise ParseError(line, "Section type {section} unsupported")
            elif first == "END":
                section = None
            elif first in ["NODES", "EDGES", "TERMINALS"]:
                # we don't need these and ignore them
                continue
            elif first == "E":
                if section != "GRAPH":
                    raise ParseError(line, "Only valid in Graph section")
                G.add_edge(int(rest[0]), int(rest[1]), weight=int(rest[2]))
            elif first == "T":
                if section != "TERMINALS":
                    raise ParseError(line, "Only valid in Terminals section")
                G.nodes[int(rest[0])]["terminal"] = True
                G.graph["terminals"] = G.graph.get(
                    "terminals", []) + [int(rest[0])]
            elif re.match("D+$", first):
                if section != "COORDINATES":
                    raise ParseError(line, "Only valid in Coordinates section")
                if len(rest) != len(first) + 1:
                    raise ParseError(line, "Incorrect number of dimensions")
                G.nodes[int(rest[0])]["pos"] = tuple([int(x)
                                                      for x in rest[1:]])
            elif section == "COMMENT":
                G.graph[line.split()[0]] = line.split("\"")[1].strip("\"")
            else:
                raise ParseError(line, "Unrecognized syntax")

    return G


if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        G = parse(f)

    print(G)
    print(G.graph)
    print(G.nodes)
    print(G.nodes.data())
    print(G.edges)

    #nx.draw(G, with_labels=True)
    #nx.draw_networkx(G, with_labels=True)
    nx.draw_kamada_kawai(G, with_labels=True)
    #nx.draw_spring(G, with_labels=True)
    #nx.draw_spectral(G, with_labels=True)
    plt.savefig("steinlib.png")
