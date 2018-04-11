from ontology import *
from decomp import *

class ASTNode(object):
    pass

class PrimFamilyNode(ASTNode):
    def __init__(self, family, paramNodes):
        if len(paramNodes) != len(family.signature):
            raise Exception("NodePrimFamily number of nodes does not match family signature.")
        self.family = family
        self.paramNodes = paramNodes
        self.size = 1 + sum([p.size for p in self.paramNodes])

    def match(self, target, matchList):
        if len(target) != 1:
            return []
        if dim(target) == 0:
            return []
        elif dim(target) == 1:
            atom1 = next(iter(target.mol1s)).atom1s[0]
            if (atom1.a0 == fMol0(())) and (atom1.b0 == fMol0(())) and isinstance(atom1.p1, FamilyPrim1):
                if atom1.p1.family == self.family:
                    newMatches = matchList
                    for ii in range(self.paramNodes):
                        newMatches = self.paramNodes[ii].match(atom1.p1.params[ii], newMatches)
                    return newMatches
            return []
        elif dim(target) == 2:
            atom2 = next(iter(target.aeMol2s)).eqAtom2s[0].rightmost
            if isinstance(atom2.l1, IdMol1) and isinstance(atom2.r1, IdMol1) and (atom2.a0 == fMol0(())) and (atom2.b0 == fMol0(())) and isinstance(atom2.p2, FamilyPrim2):
                if atom2.p2.family == self.family:
                    newMatches = matchList
                    for ii in range(self.paramNodes):
                        newMatches = self.paramNodes[ii].match(atom1.p1.params[ii], newMatches)
                    return newMatches
            return []

    def eval(self, assignment):
        params = [p.eval(assignment) for p in self.paramNodes]
        prim = self.family.buildPrim(params)
        if dim(prim) == 1:
            return fEqMol1(prim1ToMol1(prim))
        elif dim(prim) == 2:
            return fEqAEMol2(prim2ToAEMol2(prim))


class Comp0Node(ASTNode):
    def __init__(self, leftNode, rightNode):
        self.leftNode = leftNode
        self.rightNode = rightNode
        self.size = 1 + self.leftNode.size + self.rightNode.size

    def match(self, target, matchList):
        # No empty tensor decomps.
        if dim(target) == 0:
            decomps = tensorDecompMol0(target)
        elif dim(target) == 1:
            decomps = tensorDecompEqMol1(target)
        elif dim(target) == 2:
            decomps = tensorDecompEqAEMol2(target)

        newMatches = []
        for (left, right) in decomps:
            leftMatches = self.leftNode.match(left, matchList)
            rightMatches = self.rightNode.match(right, leftMatches)
            newMatches.extend(rightMatches)
        return newMatches

    def eval(self, assignment):
        return comp0(self.leftNode.eval(assignment), self.rightNode.eval(assignment))


class Comp1Node(ASTNode):
    def __init__(self, leftNode, rightNode):
        self.leftNode = leftNode
        self.rightNode = rightNode
        self.size = 1 + self.leftNode.size + self.rightNode.size

    def match(self, target, matchList):
        if dim(target) == 1:
            decomps = horizontalDecompEqMol1(target)
            decomps.add((fEqMol1(fIdMol1(target.source)), target))
            decomps.add((target, fEqMol1(fIdMol1(target.target))))
        elif dim(target) == 2:
            decomps = horizontalDecompEqAEMol2(target)
            decomps.add((fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(target.source.source)))), target))
            decomps.add((target, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(target.target.target))))))

        newMatches = []
        for (left, right) in decomps:
            leftMatches = self.leftNode.match(left, matchList)
            rightMatches = self.rightNode.match(right, leftMatches)
            newMatches.extend(rightMatches)
        return newMatches

    def eval(self, assignment):
        return comp1(self.leftNode.eval(assignment), self.rightNode.eval(assignment))


class Comp2Node(ASTNode):
    def __init__(self, leftNode, rightNode):
        self.leftNode = leftNode
        self.rightNode = rightNode
        self.size = 1 + self.leftNode.size + self.rightNode.size

    def match(self, target, matchList):
        # Target must be dimension 2
        decomps = verticalDecompEqAEMol2(target)
        decomps.add((fEqAEMol2(fAEIdMol2(target.source)), target))
        decomps.add((target, fEqAEMol2(fAEIdMol2(target.target))))

        newMatches = []
        for (left, right) in decomps:
            leftMatches = self.leftNode.match(left, matchList)
            rightMatches = self.rightNode.match(right, leftMatches)
            newMatches.extend(rightMatches)
        return newMatches

    def eval(self, assignment):
        return comp2(self.leftNode.eval(assignment), self.rightNode.eval(assignment))




class IdNode(ASTNode):
    def __init__(self, subNode):
        self.subNode = subNode
        self.size = 1 + self.subNode.size

    def match(self, target, matchList):
        if dim(target) == 0:
            return []
        elif dim(target) == 1:
            inst = next(iter(target.mol1s))
            if isinstance(inst, IdMol1):
                return self.subNode.match(inst.mol0, matchList)
            return []
        elif dim(target) == 2:
            inst = next(iter(target.aeMol2s))
            if isinstance(inst, AEIdMol2):
                return self.subNode.match(inst.eqMol1, matchList)
            return []
        return []

    def eval(self, assignment):
        base = self.subNode.eval(assignment)
        if dim(base) == 0:
            return fEqMol1(fIdMol1(base))
        if dim(base) == 1:
            return fEqAEMol2(fAEIdMol2(base))


class VarNode(ASTNode):
    # Name is an integer index when used in primitive family source/target ASTs
    def __init__(self, name):
        self.name = name
        self.size = 1

    def match(self, target, matchList):
        updatedMatches = []
        for match in matchList:
            if self.name in match:
                if target == match[self.name]:
                    updatedMatches.append(match)
                # Leave it out if it doesn't match.
            else:
                new = dict(match)
                new[self.name] = target
                updatedMatches.append(new)
        return updatedMatches

    def eval(self, assignment):
        return assignment[self.name]

class ConstNode(ASTNode):
    def __init__(self, prim):
        self.prim = prim
        self.dim = dim(prim)
        self.size = 1
        if self.dim == 0:
            self.matchTarget = prim0ToMol0(self.prim)
        elif self.dim == 1:
            self.matchTarget = fEqMol1(prim1ToMol1(self.prim))
        elif self.dim == 2:
            self.matchTarget = fEqAEMol2(prim2ToAEMol2(self.prim))

    def match(self, target, matchList):
        if target == self.matchTarget:
            return matchList
        return []

    def eval(self, assignment):
        return self.matchTarget


# Takes a list [x, y, z], and returns a dict {0:x, 1:y, 2:z}
def listToIndexDict(l):
    d = {}
    for ii in range(len(l)):
        d[ii] = l[ii]
    return d

def indexDictToList(d):
    l = []
    for ii in range(len(d)):
        l.append(d[ii])
    return l

class FamilyPrim1(Prim1):
    def __init__(self, family, params):
        self.family = family
        self.params = params
        self.source = self.family.sourceAst.eval(listToIndexDict(self.params))
        self.target = self.family.targetAst.eval(listToIndexDict(self.params))

    def __str__(self):
        return self.family.strParams(self.params)

class FamilyPrim2(Prim2):
    def __init__(self, family, params):
        self.family = family
        self.params = params
        self.source = self.family.sourceAst.eval(listToIndexDict(self.params))
        self.target = self.family.targetAst.eval(listToIndexDict(self.params))

    def __str__(self):
        return self.family.strParams(self.params)

class FamilyPrim3(Prim3):
    def __init__(self, family, params):
        self.family = family
        self.params = params
        self.source = self.family.sourceAst.eval(listToIndexDict(self.params))
        self.target = self.family.targetAst.eval(listToIndexDict(self.params))

    def __str__(self):
        return self.family.strParams(self.params)

class PrimitiveFamily(object):
    # Signature is a list of dimensions, one for each parameter slot.  E.g. [1,1,1] takes three arguments, all of dimension 1.
    # AST VarNodes are to be named after the index of each parameter, as an integer.
    def __init__(self, name, dim, signature, sourceAst, targetAst):
        self.name = name
        self.dim = dim
        self.signature = signature
        self.sourceAst = sourceAst
        self.targetAst = targetAst

    def buildPrim(self, params):
        # Need as an ordered list.
        if isinstance(params, dict):
            params = indexDictToList(params)

        if len(params) == len(self.signature):
            for ii in range(len(params)):
                if dim(params[ii]) != self.signature[ii]:
                    raise Exception("Signature error")

        if self.dim == 1:
            return FamilyPrim1(self, params)
        elif self.dim == 2:
            return FamilyPrim2(self, params)
        elif self.dim == 3:
            return FamilyPrim3(self, params)

    def strParams(self, params):
        return self.name + "(" + ", ".join(map(str, params)) + ")"

# All of the below take a molecule and attempt to make an AST from it, replacing paramPrims wherever they occur with VarNodes named by their index in the paramPrims list
def astFromMol0(mol0, paramPrims):
    if len(mol0) == 0:
        raise Exception("No ASTs from length 0 mol0's")

    if len(mol0) == 1:
        prim0 = mol0.prim0s[0]
        if prim0 in paramPrims:
            name = paramPrims.index(prim0)
            return VarNode(name)
        else:
            return ConstNode(fMol0((prim0,)))
    else:
        return Comp0Node(astFromMol0(takeMol0(mol0, 1)), astFromMol0(dropMol0(mol0, 1)))


def astFromMol1(mol1, paramPrims):
    if isinstance(mol1, IdMol1):
        return IdNode(astFromMol0(mol1.mol0))
    else:
        if len(mol1) == 1:
            atom1 = mol1.atom1s[0]
            if atom1.p1 in paramPrims:
                name = paramPrims.index(atom1.p1)
                node = VarNode(name)
            else:
                node = ConstNode(atom1.p1)
            if len(atom1.a0) > 0:
                node = Comp0Node(IdNode(astFromMol0(atom1.a0)), node)
            if len(atom1.b0) > 0:
                node = Comp0Node(node, IdNode(astFromMol0(atom1.b0)))
            return node
        else:
            return Comp1Node(astFromMol1(takeMol1(mol1, 1)), astFromMol1(dropMol1(mol1,1)))

def astFromAEMol2(aeMol2, paramPrims):
    if isinstance(aeMol2, AEIdMol2):
        return IdNode(astFromMol1(next(iter(aeMol2.eqMol1))))
    else:
        if len(aeMol2) == 1:
            atom2 = next(iter(aeMol2.eqAtom2s[0]))
            if atom2.p2 in paramPrims:
                name = paramPrims.index(atom2.p2)
                node = VarNode(name)
            else:
                node = ConstNode(atom2.p2)
            if len(atom1.a0) > 0:
                node = Comp0Node(IdNode(IdNode(astFromMol0(atom1.a0))), node)
            if len(atom1.b0) > 0:
                node = Comp0Node(node, IdNode(IdNode(astFromMol0(atom1.b0))))
            if len(atom1.l1) > 0:
                node = Comp1Node(IdNode(astFromMol1(atom1.l1)), node)
            if len(atom1.r1) > 0:
                node = Comp1Node(node, IdNode(astFromMol1(atom1.r1)))
            return node
        else:
            return Comp2Node(astFromAEMol2(takeAEMol2(aeMol2, 1)), astFromAEMol2(dropAEMol2(aeMol2,1)))

def minimalASTFromMol0(mol0, paramPrims):
    if len(mol0) == 0:
        raise Exception("No AST for empty mol0")
    elif len(mol0) == 1:
        prim0 = mol0.prim0s[0]
        if prim0 in paramPrims:
            name = paramPrims.index(prim0)
            return VarNode(name)
        else:
            return ConstNode(prim0)
    tensorDecomps = tensorDecompMol0(mol0)
    candidates = [Comp0Node(minimalASTFromMol0(l, paramPrims), minimalASTFromMol0(r, paramPrims)) for (l, r) in tensorDecomps]
    return min(candidates, key=lambda x: x.size)

def minimalASTFromEqMol1(eqMol1, paramPrims):
    if len(eqMol1) == 0:
        idMol1 = next(iter(eqMol1.mol1s))
        return IdNode(minimalASTFromMol0(idMol1.mol0, paramPrims))
    elif len(eqMol1) == 1:
        atom1 = next(iter(eqMol1.mol1s)).atom1s[0]
        if len(atom1.a0) == 0 and len(atom1.b0) == 0:
            if atom1.p1 in paramPrims:
                name = paramPrims.index(atom1.p1)
                return VarNode(name)
            else:
                return ConstNode(atom1.p1)
    tensorDecomps = tensorDecompEqMol1(eqMol1)
    horizontalDecomps = horizontalDecompEqMol1(eqMol1)

    candidates = []
    candidates.extend([Comp0Node(minimalASTFromEqMol1(l, paramPrims), minimalASTFromEqMol1(r, paramPrims)) for (l,r) in tensorDecomps])
    candidates.extend([Comp1Node(minimalASTFromEqMol1(l, paramPrims), minimalASTFromEqMol1(r, paramPrims)) for (l,r) in horizontalDecomps])
    return min(candidates, key=lambda x: x.size)


# Finds all possible ASTs for the AEMol2, then takes the minimal one.
def minimalASTFromEqAEMol2(eqAEMol2, paramPrims):
    if len(eqAEMol2) == 0:
        aeIdMol2 = next(iter(eqAEMol2.aeMol2s))
        return IdNode(minimalASTFromEqMol1(aeIdMol2.eqMol1))
    elif len(eqAEMol2) == 1:
        atom2 = next(iter(next(iter(eqAEMol2.aeMol2s)).eqAtom2s[0]))
        if len(atom2.a0) == 0 and len(atom2.b0) == 0 and len(atom2.l1) == 0 and len(atom2.r1) == 0:
            if atom2.p2 in paramPrims:
                name = paramPrims.index(atom2.p2)
                return VarNode(name)
            else:
                return ConstNode(atom2.p2)
    tensorDecomps = tensorDecompEqAEMol2(eqAEMol2)
    horizontalDecomps = horizontalDecompEqAEMol2(eqAEMol2)
    verticalDecomps = verticalDecompEqAEMol2(eqAEMol2)

    candidates = []
    candidates.extend([Comp0Node(minimalASTFromEqAEMol2(l, paramPrims), minimalASTFromEqAEMol2(r, paramPrims)) for (l,r) in tensorDecomps])
    candidates.extend([Comp1Node(minimalASTFromEqAEMol2(l, paramPrims), minimalASTFromEqAEMol2(r, paramPrims)) for (l,r) in horizontalDecomps])
    candidates.extend([Comp2Node(minimalASTFromEqAEMol2(l, paramPrims), minimalASTFromEqAEMol2(r, paramPrims)) for (l,r) in verticalDecomps])
    return min(candidates, key=lambda x: x.size)
