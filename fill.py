from ontology import *
from decomp import *

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
    return sites(eqMol1, fEqMol1(fIdMol1(x.source)), fEqMol1(fIdMol1(x.target)), horizontalDecompEqMol1)

def tensorSitesEqMol1(eqMol1):
    return sites(eqMol1, fEqMol1(fIdMol1(fMol(()))), fEqMol1(fIdMol1(fMol(()))), tensorDecompEqMol1)

def verticalSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(x.source)), fEqAEMol2(fAEIdMol2(x.target)), verticalDecompEqAEMol2)

def horizontalSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(x.source.source)))), fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(x.target.target)))), horizontalDecompEqAEMol2)

def tensorSitesEqAEMol2(eqAEMol2):
    return sites(eqAEMol2, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(fMol(()))))), fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(fMol(()))))), horizontalDecompEqAEMol2)

def findPaths(start, end, primFamilies, cellsAwayFunc, collateFunc):
    # Maps 0-cells to the 1-cells from primFamilies having them as sources.
    cellGraph = {}
    # If start == end, we don't include it in seen so that we don't preclude all paths.
    if start == end:
        paths = [{"cells": (), "seen": frozenset(), "current": start}]
    else:
        paths = [{"cells": (), "seen": frozenset([start]), "current": start}]

    finishedPaths = []
    while paths:
        newPaths = []
        for path in paths:
            if path["current"] in cellGraph:
                cellsAway = cellGraph[path["current"]]
            else:
                cellsAway = cellsAwayFunc(path["current"], primFamilies)
                cellGraph[path["current"]] = cellsAway
            for cell in cellsAway:
                if cell.target in path["seen"]:
                    continue
                p = dict(path)
                p["cells"] = p["cells"] + (cell,)
                p["seen"] = p["seen"].union(frozenset([cell.target]))
                p["current"] = cell.target
                if cell.target == end:
                    finishedPaths.append(p)
                else:
                    newPaths.append(p)
        paths = newPaths
    #return finishedPaths
    return list(set([collateFunc(tuple(p["cells"])) for p in finishedPaths]))

def cellsAway1(mol0, primFamilies):
    cellsAway = []
    for site in tensorSitesMol0(mol0):
        sitePrims = []
        for fam in primFamilies:
            matches = fam.sourceAst.match(site[1], [{}])
            matchPrims = [fam.buildPrim(m) for m in matches]
            sitePrims.extend(matchPrims)
        siteCells = [comp0s([site[0], p, site[2]]) for p in sitePrims]
        cellsAway.extend(siteCells)
    return cellsAway

def collate1(cells):
    return fEqMol1(fNonIdMol1(cells))

# TODO: Add limited use prim families.
def findPaths1(start, end, primFamilies):
    return findPaths(start, end, primFamilies, cellsAway1, collate1)

def cellsAway2(eqMol1, primFamilies):
    cellsAway = []
    for hSite in horizontalSitesEqMol1(eqMol1):
        for tSite in tensorSiteEqMol1(hSite[1]):
            sitePrims = []
            for fam in primFamilies:
                matches = fam.sourceAst.match(tSite[1], [{}])
                matchPrims = [fam.buildPrim(m) for m in matches]
                sitePrims.extend(matchPrims)
            siteCells = [comp1s([hSite[0], comp0s([tSite[0], p, tSite[2]]), hSite[2]]) for p in sitePrims]
            cellsAway.extend(siteCells)
    return cellsAway

def collate2(cells):
    return fEqAEMol2(fAENonIdMol2(cells))



# TODO: In findPaths3, what's the justification for iterative decomposition in the order we do it?
# Should be something like: every AST describe a 2-cell from primitives can be reduced to a disjunctive normal form; check that we can move higher dimensional comps up the AST.

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
                matches = fam.sourceAst.match(site[1], [{}])
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
