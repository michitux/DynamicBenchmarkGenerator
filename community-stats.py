#!/usr/bin/python3

from networkit import *
import seaborn
import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt 

parser = argparse.ArgumentParser(description="Analyze ground truth and found community structure, optionally calling moses")
parser.add_argument("graph", help="The input graph file as tab-separated edgelist with zero-based node ids")
parser.add_argument("communities", help="The ground truth community structure")
parser.add_argument("--moses", help="Path to moses binary, if provided the graph will be clustered using the provided binary")
parser.add_argument("--found-communities", dest="found", help="Path to the found community structure, this file will be written if moses is given.", required=True)
parser.add_argument("--plot-path", dest="plot", help="The path of the plot to generate", required=True)
parser.add_argument("--separator", help="The separate of the nodes in the edge list", default="tab", choices=["tab", "space"])
parser.add_argument("--truth-format", dest="format", help="The format of the ground truth cover", default="node-community", choices=["node-community", "community-nodes"])

args = parser.parse_args()

if args.separator == "space":
    G = readGraph(args.graph, Format.EdgeListSpaceZero)
else:
    G = readGraph(args.graph, Format.EdgeListTabZero)

if args.format == "node-community":
    C = graphio.EdgeListCoverReader(0).read(args.communities, G)
else:
    C = graphio.CoverReader().read(args.communities, G)

if not args.moses is None:
    import subprocess
    subprocess.call([args.moses, args.graph, args.found])

mosesC = graphio.CoverReader(0).read(args.found, G)

print("Found {found} communities of {truth} ground truth communities.".format(truth=C.numberOfSubsets(), found=mosesC.numberOfSubsets()))

foundCommunities = { i : mosesC.getMembers(i) for i in  mosesC.getSubsetIds() }
trueCommunities = { i : C.getMembers(i) for i in C.getSubsetIds() }
comp = scd.SCDGroundTruthComparison(G, mosesC, trueCommunities, True).run()

trueComms = list(trueCommunities.keys())
trueSizes = [len(trueCommunities[i]) for i in trueComms]
foundF1ofTrue = comp.getIndividualF1()
trueF1s = [foundF1ofTrue[i] for i in trueComms]
avgNumComms = [ sum([len(C.subsetsOf(u)) for u in trueCommunities[ci]]) / len(trueCommunities[ci]) for ci in trueComms ]
avgDeg = [ sum([G.degree(u) for u in trueCommunities[ci]]) / len(trueCommunities[ci]) for ci in trueComms ]

recalls = [0 for i in trueComms]
matches = [-1 for i in trueComms]
precisions = [0 for i in trueComms]
jaccards = [0 for i in trueComms]

for i, com in enumerate(trueComms):
    bestMatch = -1
    bestRecall = 0
    bestPrecision = 0
    bestJaccard = 0
    comSet = set(trueCommunities[com])
    for j, biggerCom in enumerate(trueComms):
        if i == j:
            continue

        biggerComSet = set(trueCommunities[biggerCom])

        intersectionSize = len(comSet.intersection(biggerComSet))
        recall = intersectionSize / len(comSet)
        precision = intersectionSize / len(biggerComSet)
        if recall > bestRecall or (recall == bestRecall and precision > bestPrecision):
            bestMatch = j
            bestRecall = recall
            bestPrecision = precision

        jaccard = intersectionSize / (len(comSet) + len(biggerComSet) - intersectionSize)

        if jaccard > bestJaccard:
            bestJaccard = jaccard

    recalls[i] = bestRecall
    precisions[i] = bestPrecision
    matches[i] = bestMatch
    jaccards[i] = bestJaccard

df = pd.DataFrame(np.array([trueComms, trueSizes, trueF1s, recalls, precisions, avgNumComms, avgDeg]).T, columns=["Id", "Size", "Found F1", "Contained", "Fraction of parent", "Avg. #memberships", "Avg. node degree"])

g = seaborn.PairGrid(df, diag_sharey=False)
g.map_lower(seaborn.kdeplot, cmap="Blues_d")
g.map_upper(plt.scatter, s=1)
g.map_diag(seaborn.kdeplot, lw=3)
plt.savefig(args.plot)

#foundRecallofTrue = comp.getIndividualRecall()
#trueRecalls = [foundRecallofTrue[i] for i in trueComms]
#
#foundPrecisionofTrue = comp.getIndividualPrecision()
#truePrecisions = [foundPrecisionofTrue[i] for i in trueComms]

duplicate_communities = [i for i, j in enumerate(jaccards) if j == 1]
perfect_recalls = [(i, r, p, m) for i, (r, p, m) in enumerate(zip(recalls, precisions, matches)) if r == 1.0]
print("{} communities are perfect subsets of other communities, these are {}% of all communities".format(len(perfect_recalls), 100*len(perfect_recalls)/len(trueCommunities)))
print("{} communities are duplicates".format(len(duplicate_communities)))
print("List of these communities:")
print(perfect_recalls)

badCommunities = [i for i, f in enumerate(trueF1s) if f < 0.9]
badSubsetCommunities = [i for i in badCommunities if recalls[i] == 1]

print("{} communities have a F1 score of less than 0.9. {} of these communities are subsets of each other. These are {}%.".format(len(badCommunities), len(badSubsetCommunities), 100*len(badSubsetCommunities)/len(badCommunities)))
