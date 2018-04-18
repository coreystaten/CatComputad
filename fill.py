from ontology import *
from decomp import *
from primfamily import *
import time
# Returns True if the generated prim (found from a source) is acceptable, False if not.
# Used to eliminate degenerate or expanding primitives.

# TODO: Include middle-identity sites?  I.e. (idSource, idSource, x) and (x, idTarget, idTarget).  Is this necessary for any patterns?
# Returns a list of 3-tuples, which are 3-fold tensor decompositions of mol0.  Includes trivial decompositions.
def sites(x, idSource, idTarget, decomp):
    decomps = [(idSource, x, idTarget)]
    d1 = list(decomp(x))
    decomps.extend([(idSource, l, r) for (l, r) in d1])
    decomps.extend([(l, r, idTarget) for (l, r) in d1])
    for (l,r) in d1:
        d2 = list(decomp(l))
        decomps.extend([(l1,r1,r) for (l1,r1) in d2])
    return decomps


def tensorSitesMol0(mol0):
    return sites(mol0, fMol0(()), fMol0(()), tensorDecompMol0)

def horizontalSitesEqMol1(eqMol1):
    return sites(eqMol1, fEqMol1(fIdMol1(eqMol1.source)), fEqMol1(fIdMol1(eqMol1.target)), horizontalDecompEqMol1)

def tensorSitesEqMol1(eqMol1):
    return sites(eqMol1, fEqMol1(fIdMol1(fMol0(()))), fEqMol1(fIdMol1(fMol0(()))), tensorDecompEqMol1)

def verticalSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(eqAEMol2.source)), fEqAEMol2(fAEIdMol2(eqAEMol2.target)), verticalDecompEqAEMol2)

def horizontalSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(eqAEMol2.source.source)))), fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(eqAEMol2.target.target)))), horizontalDecompEqAEMol2)

def tensorSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(fMol0(()))))), fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(fMol0(()))))), tensorDecompEqAEMol2)


def findPaths(start, keepCond, endCond, primFamilies, cellsAwayFunc, collateFunc, memory={}):
    # Maps 0-cells to the 1-cells from primFamilies having them as sources.
    cellGraph = {}

    path = {"cells": (), "seen": frozenset([start]), "current": start}
    # If start == end, we don't include it in seen so that we don't preclude all paths.
    # TODO: Got rid of end parameter, does this formulation cause any problems?
    if endCond(path) and keepCond(path):
        path["seen"] = frozenset()

    (finishedPaths, pathCount, reportSum, lastReportSum) = _stepPathNarrow(path, cellGraph, keepCond, endCond, primFamilies, cellsAwayFunc, memory, 0.0, 0.0, 1.0, 0)
    return list(set([collateFunc(tuple(p["cells"])) for p in finishedPaths]))

def _stepPathNarrow(path, cellGraph, keepCond, endCond, primFamilies, cellsAwayFunc, memory, reportSum, lastReportSum, reportWeight, pathCount):
#print(list(map(str, paths[-1]["cells"])))
#        if len(paths) == 8:
#            for p in paths:
#                print(list(map(str, p["cells"])))
    #if len(paths[0]["cells"]) > 10:
        #print(str(paths[0]["cells"][2]))
    if path["current"] in cellGraph:
        cellsAway = cellGraph[path["current"]]
    else:
        cellsAway = cellsAwayFunc(path["current"], primFamilies, memory)
        cellGraph[path["current"]] = cellsAway
    (finished, new) = _stepPath(path, cellsAway, keepCond, endCond, primFamilies, cellsAwayFunc, memory)
    if len(new) == 0:
        returnPc = 1
        reportSum += reportWeight
    else:
        returnPc = 0

    for n in new:
        subWeight = reportWeight / len(new)
        (newFinished, newPc, reportSum, lastReportSum) = _stepPathNarrow(n, cellGraph, keepCond, endCond, primFamilies, cellsAwayFunc, memory, reportSum, lastReportSum, subWeight, pathCount)
        finished.extend(newFinished)
        pathCount += newPc
        returnPc += newPc
        if int(reportSum  * 100) > int(lastReportSum * 100):
            #print("%d%% progress after %d paths:" % (int(reportSum * 100), pathCount))
            lastReportSum = reportSum

    return (finished, returnPc, reportSum, lastReportSum)

# TODO: Add adjoint prim matching.
# keepCond and endCond should be True/False functions on paths.
# If keepCond is True, a record of the path will be returned.  If endCond is True, that path will no longer be explored.
def findPathsBroad(start, keepCond, endCond, primFamilies, cellsAwayFunc, collateFunc, memory={}):
    # Maps 0-cells to the 1-cells from primFamilies having them as sources.
    cellGraph = {}

    paths = [{"cells": (), "seen": frozenset([start]), "current": start}]
    # If start == end, we don't include it in seen so that we don't preclude all paths.
    # TODO: Got rid of end parameter, does this formulation cause any problems?
    if endCond(paths[0]) and keepCond(paths[0]):
        paths[0]["seen"] = frozenset()

    finishedPaths = []
    while paths:
        (newFinished, newPaths) = _stepPaths(paths, cellGraph, keepCond, endCond, primFamilies, cellsAwayFunc, memory)
        finishedPaths.extend(newFinished)
        paths = newPaths
    return list(set([collateFunc(tuple(p["cells"])) for p in finishedPaths]))

def _stepPathsBroad(paths, cellGraph, keepCond, endCond, primFamilies, cellsAwayFunc, memory):
    if len(paths) > 10000:
         print("TOO MANY PATHS.  Quitting search.")
         return ([], [])
        #print(list(map(str, paths[-1]["cells"])))
    print("(_stepPaths: %d)" % (len(paths),))
#        if len(paths) == 8:
#            for p in paths:
#                print(list(map(str, p["cells"])))
    #if len(paths[0]["cells"]) > 10:
        #print(str(paths[0]["cells"][2]))
    finishedPaths = []
    newPaths = []
    for path in paths:
        if path["current"] in cellGraph:
            cellsAway = cellGraph[path["current"]]
        else:
            cellsAway = cellsAwayFunc(path["current"], primFamilies, memory)
            cellGraph[path["current"]] = cellsAway
        (finished, new) = _stepPath(path, cellsAway, keepCond, endCond, primFamilies, cellsAwayFunc, memory)
        finishedPaths.extend(finished)
        newPaths.extend(new)
    return (finishedPaths, newPaths)

def _stepPath(path, cellsAway, keepCond, endCond, primFamilies, cellsAwayFunc, memory):
    finishedPaths = []
    newPaths = []

    for cell in cellsAway:
        if cell.target in path["seen"]:
            continue
        q = _checkSupp(path,cell)
        p = _buildNewPath(path, cell)

        if keepCond(p):
            finishedPaths.append(p)
        if not(endCond(p)):
            newPaths.append(p)
    return (finishedPaths, newPaths)

def _checkSupp(path, cell):
    p = dict(path)
    p["cells"] = p["cells"] + (cell,)
    p["seen"] = p["seen"].union(frozenset([cell.target]))
    return p

def _buildNewPath(path, cell):
    p = dict(path)
    p["cells"] = p["cells"] + (cell,)
    p["seen"] = p["seen"].union(frozenset([cell.target]))
    p["current"] = cell.target
    return p

def cellsAway1(mol0, primFamilies, memory={}):
    # Path-cache is used to keep track of path graphs generated while looking under FunctorPrims
    if "pc" not in memory:
        memory["pc"] = {}

    cellsAway = []
    for site in tensorSitesMol0(mol0):
        sitePrims = []
        for fam in primFamilies:
            matches = fam.sourceAST.match(site[1], [ASTMatch()])
            matchPrims = [fam.buildPrim(indexDictToList(m.varMatch), m.functorMatches) for m in matches if not(fam.isDegen(m.varMatch))]
            sitePrims.extend(matchPrims)
        q = asPrim0(site[1])
        if q is not None and isinstance(q, FunctorPrim0):
            if q.mol0 in memory["pc"]:
                cells = memory["pc"][q.mol0]
            else:
                cells = findPaths(q.mol0, lambda x: True, lambda x: False, primFamilies, cellsAway1, collate1)
                memory["pc"][q.mol0] = cells
            for cell in cells:
                sitePrims.append(fFunctorPrim1(q.functor, cell))
        siteCells = [comp0s(site[0], p, site[2]) for p in sitePrims]
        cellsAway.extend(siteCells)
    return cellsAway

def collate1(cells):
    return fEqMol1(fNonIdMol1(cells))

def targetEndCond(target):
    return lambda x: x["current"] == target

# TODO: Add limited use prim families.
def findPaths1(start, end, primFamilies):
    return findPaths(start, targetEndCond(end), targetEndCond(end), primFamilies, cellsAway1, collate1)

def cellsAway2(eqMol1, primFamilies, memory={}):
    if "pc" not in memory:
        memory["pc"] = {}

    cellsAway = []
    for hSite in horizontalSitesEqMol1(eqMol1):
        for tSite in tensorSitesEqMol1(hSite[1]):
            sitePrims = []
            for fam in primFamilies:
                matches = fam.sourceAST.match(tSite[1], [ASTMatch()])
                matchPrims = [fam.buildPrim(indexDictToList(m.varMatch), m.functorMatches) for m in matches if not(fam.isDegen(m.varMatch))]
                sitePrims.extend(matchPrims)
            q = asPrim1(tSite[1])
            if q is not None and isinstance(q, FunctorPrim1):
                if q.eqMol1 in memory["pc"]:
                    cells = memory["pc"][q.eqMol1]
                else:
                    cells = findPaths(q.eqMol1, lambda x: True, lambda x: False, primFamilies, cellsAway2, collate2)
                    memory["pc"][q.eqMol1] = cells
                for cell in cells:
                    sitePrims.append(fFunctorPrim2(q.functor, cell))
            siteCells = [comp1s(hSite[0], comp0s(tSite[0], p, tSite[2]), hSite[2]) for p in sitePrims]
            cellsAway.extend(siteCells)
    return cellsAway

def collate2(cells):
    return fEqAEMol2(fAENonIdMol2(cells))


def findPaths2(start, end, primFamilies):
    return findPaths(start, targetEndCond(end), targetEndCond(end), primFamilies, cellsAway2, collate2)

# TODO: What's the justification for iterative decomposition in the order we do it?
# Should be something like: every AST describe a 2-cell from primitives can be reduced to a disjunctive normal form; check that we can move higher dimensional comps up the AST.
def cellsAway3(eqAEMol2, primFamilies, memory={}):
    if "pc" not in memory:
        memory["pc"] = {}

    cellsAway = []
    for vSite in verticalSitesEqAEMol2(eqAEMol2):
        for hSite in horizontalSitesEqAEMol2(vSite[1]):
            for tSite in tensorSitesEqAEMol2(hSite[1]):
                sitePrims = []
                for fam in primFamilies:
                    matches = fam.sourceAST.match(tSite[1], [ASTMatch()])
                    matchPrims = [fam.buildPrim(indexDictToList(m.varMatch), m.functorMatches) for m in matches if not(fam.isDegen(m.varMatch))]
                    sitePrims.extend(matchPrims)
                q = asPrim2(tSite[1])
                if q is not None and isinstance(q, FunctorPrim2):
                    if q.eqAEMol2 in memory["pc"]:
                        cells = memory["pc"][p.eqAEMol2]
                    else:
                        cells = findPaths(q.eqAEMol2, lambda x: True, lambda x: False, primFamilies, cellsAway3, collate3)
                        memory["pc"][q.eqAEMol2] = cells
                    for cell in cells:
                        sitePrims.append(fFunctorPrim3(q.functor, cell))
                siteCells = [StepCell3(vSite[0], hSite[0], tSite[0], p, tSite[2], hSite[2], vSite[2]) for p in sitePrims]
                cellsAway.extend(siteCells)
    return cellsAway

def collate3(cells):
    return CompCell3(cells)

def findPaths3(start, end, primFamilies):
    return findPaths(start, targetEndCond(end), targetEndCond(end), primFamilies, cellsAway3, collate3)




# Takes a starting and ending mol0, returns a dictionary mapping mol0s to a list of eqMol1s coming out of it from the prim families.
def buildPathGraph1DEFUNCT(start, primFamilies):
    graph = {}
    known = set([start])
    unexplored = [start]
    while unexplored:
        nextStart = unexplored[0]
        cells = []
        for site in tensorSitesMol0(nextStart):
            sitePrims = []
            for fam in primFamilies:
                matches = fam.sourceAST.match(site[1], [ASTMatch()])
                matchPrims = [fam.buildPrim(m) for m in matches]
                sitePrims.extend(matchPrims)
            siteCells = [comp0s(site[0], p, site[2]) for p in sitePrims]
            cells.extend(siteCells)
        graph[nextStart] = cells
        toExplore = set([c.target for c in cells]).difference(known)
        unexplored.extend(list(toExplore))
        known.update(toExplore)
    return graph


# Paths have the following data:
#  cells: a list of the cells followed so far
#  seen: a set of source/targets that have already occurred in the path (shouldn't revisit).  Exception: if start == end, don't count start as being seen.
#  current: the current source/target
# TODO: Do we have any primitives with the same source and target?  Those might not show up, even when limited.
# TODO: Add limitations on using of certain prim families.
def findPaths1DEFUNCT(start, end, pathGraph):
    if start == end:
        paths = [{"cells": [], "seen": set(), "current": start}]
    else:
        paths = [{"cells": [], "seen": set(start), "current": start}]

    finishedPaths = []
    while paths:
        newPaths = []
        for path in paths:
            cellsAway = pathGraph[path[current]]
            for cell in cellsAway:
                if cell.target in path["seen"]:
                    continue
                p = dict(path)
                p["cells"] += cell
                p["seen"].add(cell.target)
                p["current"] = cell.target
                if cell.target == end:
                    finishedPaths.append(p)
                else:
                    newPaths.append(p)

# Takes a list of lists, each inner list representing a range of options.
# Returns a list of lists, each inner list representing a possible selection of such options.
# The length of each returned inner list is the same as the length of the initial outer list.
def probConcat(ls):
    if len(ls) == 0:
        return [[]]
    elif len(ls) == 1:
        return [[y] for y in ls[0]]
    else:
        return [[y] + ys for y in ls[0] for ys in probConcat(ls[1:])]

# For each family in adjFamilies, chooses either the family or its adjoint.  Tries search for all possible such selections (2^len(adjFamilies)).
def searchForPathPairs2AdjFam(paths, adjFamilies, primFamilies):
    l = [[f, f.adj] for f in adjFamilies]
    choices = probConcat(l)
    for c in choices:
        print("Search using %s" % (str(list(map(str, c))),))
        searchForPathPairs2(paths, c + primFamilies)

# Takes a list of 1-cell paths.  Pairwise tries to find paths between them using given primitive families.  Prints results.
def searchForPathPairs2(paths, primFamilies):
    kk = len(paths)
    for ii in range(kk):
        for jj in range(kk):
            if ii == jj:
                continue
            #pathsFound = findPaths(paths[jj], targetEndCond(paths[ii]), lambda x: len(x["cells"]) > 4, primFamilies, cellsAway2, collate2)
            pathsFound = findPaths2(paths[jj], paths[ii], primFamilies)
            print("(%d, %d)" % (jj, ii))
            if len(pathsFound) >= 2:
                print("%s -- %s" % (str(paths[jj]), str(paths[ii])))
                print("%d paths found." % (len(pathsFound),))
                for path in pathsFound:
                    print("\t" + str(path))

def searchForPathPairs3AdjFam(group1, group2, adjFamilies, primFamilies):
    l = [[f, f.adj] for f in adjFamilies]
    choices = probConcat(l)
    for c in choices:
        print("Search using %s" % (str(list(map(str, c))),))
        searchForPathPairs3(group1, group2,  c + primFamilies)

def searchForPathPairs3(group1, group2, primFamilies):
    for ii in range(len(group1)):
        for jj in range(len(group2)):
            if group1[ii] == group2[jj]:
                continue
            #pathsFound = findPaths(paths[jj], targetEndCond(paths[ii]), lambda x: len(x["cells"]) > 4, primFamilies, cellsAway2, collate2)
            pathsFound = findPaths3(group1[ii], group2[jj], primFamilies)
            print("(%d, %d)" % (ii, jj))
            if len(pathsFound) >= 2:
                print("%s -- %s" % (str(group1[ii]), str(group2[jj])))
                print("%d paths found." % (len(pathsFound),))
                diffPrimSets = set(map(lambda x: frozenset(primSetCompCell3(x)), pathsFound))
                if len(diffPrimSets) > 1:
                    for path in pathsFound:
                        print("\t" + str(path))
                else:
                    print("Not distinct.  Primset: %s" % (str(lmap(str,(list(diffPrimSets)[0])))))

def searchForPathPairs3AdjFam(group1, group2, adjFamilies, primFamilies):
    l = [[f, f.adj] for f in adjFamilies]
    choices = probConcat(l)
    for c in choices:
        print("Search using %s" % (str(list(map(str, c))),))
        searchForPathPairs3(group1, group2,  c + primFamilies)

def exploreFromToAdjCells2(baseCells, seeking, adjFamilies, primFamilies, countFamilies):
    for b in baseCells:
        print("-----------------------------------")
        print("START CELL: " + str(b))
        print("-----------------------------------")
        seeking2 = set(seeking)
        if b in seeking2:
            seeking2.remove(b)
        exploreFromToAdj2(b, seeking2, adjFamilies, primFamilies, countFamilies)

def exploreFromToAdj2(base, seeking, adjFamilies, primFamilies, countFamilies):
    l = [[f, f.adj] for f in adjFamilies]
    choices = probConcat(l)
    for c in choices:
        print("Explore using %s" % (str(list(map(str, c))),))
        exploreFromTo2(base, seeking,  c + primFamilies, countFamilies)

def exploreFromToAdjCells3(baseCells, seeking, adjFamilies, primFamilies, countFamilies):
    for b in baseCells:
        print("-----------------------------------")
        print("START CELL: " + str(b))
        print("-----------------------------------")
        seeking2 = set(seeking)
        if b in seeking2:
            seeking2.remove(b)
        exploreFromToAdj3(b, seeking2, adjFamilies, primFamilies, countFamilies)

def exploreFromToAdj3(base, seeking, adjFamilies, primFamilies, countFamilies):
    l = [[f, f.adj] for f in adjFamilies]
    choices = probConcat(l)
    for c in choices:
        print("Explore using %s" % (str(list(map(str, c))),))
        exploreFromTo3(base, seeking,  c + primFamilies, countFamilies)

def exploreFromTo2(base, seeking, primFamilies, countFamilies):
    paths = findPaths(base, lambda x: ensureEqMol1(x["cells"][-1].target) in seeking, lambda x: False, primFamilies, cellsAway2, collate2, memory={})
    print("%d paths found." % len(paths))
    primCountsByTarget = {}
    countFamilyStrs = set(map(str, countFamilies))
    for p in paths:
        key = ensureEqMol1(p.target)
        if key not in primCountsByTarget:
            primCountsByTarget[key] = set()
        pcount = primCount2(p)
        for (name,t) in list(pcount.keys()):
            if name not in countFamilyStrs:
                del pcount[(name, t)]
        primCountsByTarget[key].add(frozenset(pcount.items()))
    for (t,p) in primCountsByTarget.items():
        if len(p) > 1:
            for p1 in paths:
                for p2 in paths:
                    if ensureEqMol1(p1.source) == ensureEqMol1(p2.source) and ensureEqMol1(p1.target) == ensureEqMol1(p2.target) and frozenset(primCount2(p1).items()) != frozenset(primCount2(p2).items()):
                        print("Dual paths detected:")
                        print("FROM:")
                        print(p1.source)
                        print("--")
                        print(p1.target)
                        print("PRIMSET 1 IS:")
                        print(lmap(str, primSetEqAEMol2(p1)))
                        print("PRIMSET 2 IS:")
                        print(lmap(str, primSetEqAEMol2(p2)))
                        print("PATH 1 IS:")
                        print(str(p1))
                        print("PATH 2 IS:")
                        print(str(p2))
            break

def exploreFromTo3(base, seeking, primFamilies, countFamilies):
    paths = findPaths(base, lambda x: ensureEqAEMol2(x["cells"][-1].target) in seeking, lambda x: False, primFamilies, cellsAway3, collate3, memory={})
    print("%d paths found." % len(paths))
    primCountsByTarget = {}
    countFamilyStrs = set(map(str, countFamilies))
    for p in paths:
        key = ensureEqAEMol2(p.target)
        if key not in primCountsByTarget:
            primCountsByTarget[key] = set()
        pcount = primCount3(p)
        for (name,t) in list(pcount.keys()):
            if name not in countFamilyStrs:
                del pcount[(name, t)]
        primCountsByTarget[key].add(frozenset(pcount.items()))
    for (t,p) in primCountsByTarget.items():
        if len(p) > 1:
            for p1 in paths:
                for p2 in paths:
                    if ensureEqAEMol2(p1.source) == ensureEqAEMol2(p2.source) and ensureEqAEMol2(p1.target) == ensureEqAEMol2(p2.target) and frozenset(primCount3(p1).items()) != frozenset(primCount3(p2).items()):
                        print("Dual paths detected:")
                        print("FROM:")
                        print(p1.source)
                        print("--")
                        print(p1.target)
                        print("PRIMSET 1 IS:")
                        print(lmap(str, primSetCompCell3(p1)))
                        print("PRIMSET 2 IS:")
                        print(lmap(str, primSetCompCell3(p2)))
                        print("PATH 1 IS:")
                        print(str(p1))
                        print("PATH 2 IS:")
                        print(str(p2))
            break



# Alternative algorithm for finding path pairs:
#   Do a simultaneous meet in the middle search for each path:
#     Take paths of length $n$ away from each path, separately using the families and their adjoints.
#     Check if there is any codomain overlap between any adjoint paths and non-adjoint paths; if so, we can turn this into a desired path.
def meetInMiddleSearchForPathPairs3(paths, primFamilies):
    print(list(map(len, paths)))
    adjFamilies = [p.adj for p in primFamilies]
    baseTime = time.time()
    print("Starting with %d 1-paths" % (len(paths),))
    soFar = {}
    work = 0
    for p in paths:
        awayNorm = list(map(lambda x: [x], cellsAway3(p, primFamilies)))
        awayAdj = list(map(lambda x: [x], cellsAway3(p, adjFamilies)))
        soFar[p] = (awayNorm, awayAdj)
        work = work + 1
        if work % 100 == 0:
            print(str(work) + ": " + str(int(time.time() - baseTime)) + " seconds")
    length = 2
    alternate = 0
    while True:
        print("Current length: %d" % (length,))
        print("Current # paths: " + str((sum(map(len, [p[0] for p in soFar.values()])), sum(map(len, [p[0] for p in soFar.values()])))))
        normCods = set()
        for x in soFar.values():
            for y in x[0]:
                normCods.add(ensureEqAEMol2(y[-1].target))
        adjCods = set()
        for x in soFar.values():
            for y in x[0]:
                adjCods.add(ensureEqAEMol2(y[-1].target))

        intersection = normCods.intersection(adjCods)
        if len(intersection) > 0:
            normByCod = {}
            adjByCod = {}
            normMatches = [p for p in x[0] if (ensureEqAEMol2(p[-1].target) in intersection) for x in soFar.values() ]
            adjMatches = [p for p in x[1] if (ensureEqAEMol2(p[-1].target) in intersection) for x in soFar.values() ]

            for ii in intersection:
                normByCod[ii] = []
                adjByCod[ii] = []
                for p in normMatches:
                    if ensureEqAEMol2(p[-1].target) == ii:
                        normByCod[ii].append(p)
                for p in adjMatches:
                    if ensureEqAEMol2(p[-1].target) == ii:
                        adjByCod[ii].append(p)
            print("Matches found at length %d:" % (length,))
            for ii in intersection:
                print("----------------------")
                for p in normByCod[ii]:
                    print("### " + str(lmap(str, p)))
                print("Can pair with")
                for p in adjByCod[ii]:
                    print("### " + str(lmap(str, p)))
        if alternate == 0:
            for p in paths:
                curPaths = soFar[p][0]
                newPaths = []
                for cur in curPaths:
                    away = cellsAway3(ensureEqAEMol2(cur[-1].target), primFamilies)
                    newPaths.extend([cur + [a] for a in away])
                    work = work + 1
                    if work % 100 == 0:
                        print(str(work) + ": " + str(int(time.time() - baseTime)) + " seconds")
                soFar[p] = (newPaths, soFar[p][1])
            alternate = 1
        elif alternate == 1:
            for p in paths:
                curPaths = soFar[p][1]
                newPaths = []
                for cur in curPaths:
                    away = cellsAway3(ensureEqAEMol2(cur[-1].target), adjFamilies)
                    newPaths.extend([cur + [a] for a in away])
                    work = work + 1
                    if work % 100 == 0:
                        print(str(work) + ": " + str(int(time.time() - baseTime)) + " seconds")
                soFar[p] = (soFar[p][0], newPaths)
            alternate = 0
        length += 1

            # Figure out which things can join.

def _incDict(tup1, tup2, d):
    if tup1 in d:
        d[tup1] += 1
    elif tup2 in d:
        d[tup2] -= 1
        if d[tup2] == 0:
            del d[tup2]
    else:
        d[tup1] = 1

def primCount2(eqAEMol2):
    return {(str(k1),k2):v for ((k1,k2),v) in _primCount2(eqAEMol2).items()}

def _primCount2(eqAEMol2):
    count = {}
    for eqAtom2 in next(iter(eqAEMol2.aeMol2s)).eqAtom2s:
        p = eqAtom2.righthand.p2
        if isinstance(p, FunctorPrim2):
            toAdd = _primCount2(p.eqAEMol2)
            for ((fam, par), v) in toAdd.items():
                for ii in range(v):
                    _incDict((fam, par), (fam.adj, par), count)
        else:
            tup1 = (p.family, tuple(p.params))
            tup2 = (p.family.adj, tuple(p.params))
            _incDict(tup1, tup2, count)
    return count

def primCount3(cc3):
    count = {}
    for stepCell in cc3.stepCell3s:
        tup1 = (str(stepCell.p.family), tuple(stepCell.p.params))
        tup2 = (str(stepCell.p.family.adj), tuple(stepCell.p.params))
        if tup1 in count:
            count[tup1] += 1
        elif tup2 in count:
            count[tup2] -= 1
            if count[tup2] == 0:
                del count[tup2]
        else:
            count[tup1] = 1
    return count


def circleAround3(startPath, primFamilies):
    paths = findPaths(startPath, lambda x: True, lambda x: len(x["cells"]) > 6, primFamilies, cellsAway3, collate3)
    pathsByCod = {}
    for p in paths:
        if p.target not in pathsByCod:
            pathsByCod[p.target] = []
        pathsByCod[p.target].append(p)
    for target in pathsByCod.keys():
        pcSets = set(map(primCount, pathsByCod[target]))
        if len(pcSets) > 1:
            ps = pathsByCod[target]
            ps.sort(key=len)
            ps.reverse()
            for p1 in ps:
                for p2 in ps:
                    if primCount(p1) != primCount(p2):
                        return (p1, p2)
            print(lmap(str, pcSets))
