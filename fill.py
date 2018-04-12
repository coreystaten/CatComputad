from ontology import *
from decomp import *
from primfamily import *

# Returns True if the generated prim (found from a source) is acceptable, False if not.
# Used to eliminate degenerate or expanding primitives.

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
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(x.source)), fEqAEMol2(fAEIdMol2(x.target)), verticalDecompEqAEMol2)

def horizontalSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(x.source.source)))), fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(x.target.target)))), horizontalDecompEqAEMol2)

def tensorSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(fMol(()))))), fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(fMol(()))))), horizontalDecompEqAEMol2)


# TODO: Implement finding paths under functors; cache results, and generate from all path results in cellsAway.
# keepCond and endCond should be True/False functions on paths.
# If keepCond is True, a record of the path will be returned.  If endCond is True, that path will no longer be explored.
def findPaths(start, keepCond, endCond, primFamilies, cellsAwayFunc, collateFunc, memory={}):
    # Maps 0-cells to the 1-cells from primFamilies having them as sources.
    cellGraph = {}

    paths = [{"cells": (), "seen": frozenset([start]), "current": start}]
    # If start == end, we don't include it in seen so that we don't preclude all paths.
    # TODO: Got rid of end parameter, does this formulation cause any problems?
    if endCond(paths[0]) and keepCond(paths[0]):
        paths[0]["seen"] = frozenset()

    finishedPaths = []
    while paths:
#        if len(paths) == 8:
#            for p in paths:
#                print(list(map(str, p["cells"])))
        #if len(paths[0]["cells"]) > 10:
            #print(str(paths[0]["cells"][2]))
        newPaths = []
        for path in paths:
            if path["current"] in cellGraph:
                cellsAway = cellGraph[path["current"]]
            else:
                cellsAway = cellsAwayFunc(path["current"], primFamilies, memory)
                cellGraph[path["current"]] = cellsAway
            for cell in cellsAway:
                if cell.target in path["seen"]:
                    continue
                p = dict(path)
                p["cells"] = p["cells"] + (cell,)
                p["seen"] = p["seen"].union(frozenset([cell.target]))
                p["current"] = cell.target
                if keepCond(p):
                    finishedPaths.append(p)
                if not(endCond(p)):
                    newPaths.append(p)
        paths = newPaths
    #return finishedPaths
    return list(set([collateFunc(tuple(p["cells"])) for p in finishedPaths]))

def cellsAway1(mol0, primFamilies, memory={}):
    # Path-cache is used to keep track of path graphs generated while looking under FunctorPrims
    if "pc" not in memory:
        memory["pc"] = {}

    cellsAway = []
    for site in tensorSitesMol0(mol0):
        sitePrims = []
        for fam in primFamilies:
            matches = fam.sourceAST.match(site[1], [ASTMatch()])
            matchPrims = [fam.buildPrim(m.varMatch) for m in matches if not(fam.isDegen(m.varMatch))]
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
        siteCells = [comp0s([site[0], p, site[2]]) for p in sitePrims]
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
                matchPrims = [fam.buildPrim(m.varMatch) for m in matches if not(fam.isDegen(m.varMatch))]
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
            siteCells = [comp1s([hSite[0], comp0s([tSite[0], p, tSite[2]]), hSite[2]]) for p in sitePrims]
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
                    matchPrims = [fam.buildPrim(m.varMatch) for m in matches if not(fam.isDegen(m.varMatch))]
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
            siteCells = [comp0s([site[0], p, site[2]]) for p in sitePrims]
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
        return []
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

# TODO: Split primFamilies into two variables, check all combinations of adjs in one.
# Takes a list of 1-cell paths.  Pairwise tries to find paths between them using given primitive families.  Prints results.
def searchForPathPairs2(paths, primFamilies):
    kk = len(paths)
    for ii in range(kk):
        for jj in range(kk):
            if ii == jj:
                continue
            pathsFound = findPaths2(paths[jj], paths[ii], primFamilies)
            if len(pathsFound) > 0:
                print("(%d,%d) -- %s -- %s" % (jj, ii, str(paths[jj]), str(paths[ii])))
                print("%d paths found." % (len(pathsFound),))
                for path in pathsFound:
                    print("\t" + str(path))
