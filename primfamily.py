from ontology import *
from decomp import *

class ASTNode(object):
    pass

class ASTMatch(object):
    def __init__(self):
        # Dictionary name => cell
        self.varMatch = {}
        # Dictionary (name, "s" or "t", depth) => cell
        self.stInfo = {}
        self.functorMatches = {}

    def dup(self):
        m = ASTMatch()
        m.varMatch = dict(self.varMatch)
        m.stInfo = dict(self.stInfo)
        m.functorMatches = dict(self.functorMatches)
        return m

class PrimitiveFamilyNode(ASTNode):
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
                    for ii in range(len(self.paramNodes)):
                        newMatches = self.paramNodes[ii].match(atom1.p1.params[ii], newMatches)
                    return newMatches
            return []
        elif dim(target) == 2:
            atom2 = next(iter(target.aeMol2s)).eqAtom2s[0].righthand
            if isinstance(atom2.l1, IdMol1) and isinstance(atom2.r1, IdMol1) and (atom2.a0 == fMol0(())) and (atom2.b0 == fMol0(())) and isinstance(atom2.p2, FamilyPrim2):
                if atom2.p2.family == self.family:
                    newMatches = matchList
                    for ii in range(len(self.paramNodes)):
                        newMatches = self.paramNodes[ii].match(atom2.p2.params[ii], newMatches)
                    return newMatches
            return []

    def eval(self, assignment, functors=[]):
        params = [p.eval(assignment, functors) for p in self.paramNodes]
        prim = self.family.buildPrim(params, functors)
        return prim
        #if dim(prim) == 1:
        #    return fEqMol1(prim1ToMol1(prim))
        #elif dim(prim) == 2:
        #    return fEqAEMol2(prim2ToAEMol2(prim))

    # Func should take a node, and return either a new node or None.  If None, recurse.
    def alterAST(self, func):
        newParamNodes = []
        for p in paramNodes:
            res = func(p)
            if res is None:
                newParamNodes.append(p.alterAST(func))
            else:
                newParamNodes.append(res)
        return PrimitiveFamilyNode(self.family, newParamNodes)

    def __str__(self):
        return self.family.name + "(" + ", ".join(map(str, self.paramNodes)) + ")"


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

    def eval(self, assignment, functors):
        return comp0(self.leftNode.eval(assignment, functors), self.rightNode.eval(assignment, functors))

    def alterAST(self, func):
        leftRes = func(self.leftNode)
        rightRes = func(self.rightNode)
        if leftRes is None:
            newLeft = self.leftNode.alterAST(func)
        else:
            newLeft = leftRes
        if rightRes is None:
            newRight = self.rightNode.alterAST(func)
        else:
            newRight = rightRes
        return Comp0Node(newLeft, newRight)

    def __str__(self):
        return "(" + str(self.leftNode) + " @ " + str(self.rightNode) + ")"


class Comp1Node(ASTNode):
    def __init__(self, leftNode, rightNode):
        self.leftNode = leftNode
        self.rightNode = rightNode
        self.size = 1 + self.leftNode.size + self.rightNode.size

    def match(self, target, matchList):
        if dim(target) == 1:
            decomps = set(horizontalDecompEqMol1(target))
            decomps.add((fEqMol1(fIdMol1(target.source)), target))
            decomps.add((target, fEqMol1(fIdMol1(target.target))))
        elif dim(target) == 2:
            decomps = set(horizontalDecompEqAEMol2(target))
            decomps.add((fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(target.source.source)))), target))
            decomps.add((target, fEqAEMol2(fAEIdMol2(fEqMol1(fIdMol1(target.target.target))))))

        newMatches = []
        for (left, right) in decomps:
            leftMatches = self.leftNode.match(left, matchList)
            rightMatches = self.rightNode.match(right, leftMatches)
            newMatches.extend(rightMatches)
        return newMatches

    def eval(self, assignment, functors):
        return comp1(self.leftNode.eval(assignment, functors), self.rightNode.eval(assignment, functors))

    def alterAST(self, func):
        leftRes = func(self.leftNode)
        rightRes = func(self.rightNode)
        if leftRes is None:
            newLeft = self.leftNode.alterAST(func)
        else:
            newLeft = leftRes
        if rightRes is None:
            newRight = self.rightNode.alterAST(func)
        else:
            newRight = rightRes
        return Comp1Node(newLeft, newRight)

    def __str__(self):
        return "(" + str(self.leftNode) + " . " + str(self.rightNode) + ")"


class Comp2Node(ASTNode):
    def __init__(self, leftNode, rightNode):
        self.leftNode = leftNode
        self.rightNode = rightNode
        self.size = 1 + self.leftNode.size + self.rightNode.size

    def match(self, target, matchList):
        # Target must be dimension 2
        decomps = set(verticalDecompEqAEMol2(target))
        decomps.add((fEqAEMol2(fAEIdMol2(target.source)), target))
        decomps.add((target, fEqAEMol2(fAEIdMol2(target.target))))

        newMatches = []
        for (left, right) in decomps:
            leftMatches = self.leftNode.match(left, matchList)
            rightMatches = self.rightNode.match(right, leftMatches)
            newMatches.extend(rightMatches)
        return newMatches

    def eval(self, assignment, functors):
        return comp2(self.leftNode.eval(assignment, functors), self.rightNode.eval(assignment, functors))

    def alterAST(self, func):
        leftRes = func(self.leftNode)
        rightRes = func(self.rightNode)
        if leftRes is None:
            newLeft = self.leftNode.alterAST(func)
        else:
            newLeft = leftRes
        if rightRes is None:
            newRight = self.rightNode.alterAST(func)
        else:
            newRight = rightRes
        return Comp2Node(newLeft, newRight)

    def __str__(self):
        return "(" + str(self.leftNode) + " & " + str(self.rightNode) + ")"




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

    def eval(self, assignment, functors):
        base = self.subNode.eval(assignment, functors)
        if dim(base) == 0:
            return fIdMol1(base)
        if dim(base) == 1:
            return fAEIdMol2(ensureEqMol1(base))

    def alterAST(self, func):
        res = func(self.subNode)
        if res is None:
            return IdNode(self.subNode.alterAST(func))
        else:
            return IdNode(res)

    def __str__(self):
        return "1_{" + str(self.subNode) + "}"

def _repeat(f, x, n):
    for ii in range(n):
        x = f(x)
    return x

# Checks if it is ok to set the given depth's source or target of the given var to cell
def _isCompatibleSetSourceTarget(name, sOrT, depth, cell, match):
    if name in match.varMatch:
        if sOrT == "s":
            source = _repeat(lambda x: x.source, match.varMatch[name], depth)
            if source != cell:
                return False
        elif sOrT == "t":
            target = _repeat(lambda x: x.target, match.varMatch[name], depth)
            if target != cell:
                return False
    for key in match.stInfo:
        (name2, sOrT2, depth2) = key
        if name == name2:
            if depth == depth2 and sOrT == sOrT2:
                if cell != match.stInfo[key]:
                    return False
            elif depth < depth2:
                increment = depth2 - depth
                if sOrT2 == "s":
                    source = _repeat(lambda x: x.source, cell, increment)
                    if source != match.stInfo[key]:
                        return False
                elif sOrT2 == "t":
                    target = _repeat(lambda x: x.target, cell, increment)
                    if target != match.stInfo[key]:
                        return False
            elif depth > depth2:
                increment = depth - depth2
                if sOrT2 == "s":
                    source = _repeat(lambda x: x.source, match.stInfo[key], increment)
                    if source != cell:
                        return False
                elif sOrT2 == "t":
                    target = _repeat(lambda x: x.target, match.stInfo[key], increment)
                    if target != cell:
                        return False
    return True

# Checks if is ok set name=cell in the match.
def _isCompatibleSetCell(name, cell, match):
    if name in match.varMatch:
        if cell != match.varMatch[name]:
            return False
    for key in match.stInfo:
        (n, sOrT, depth) = key
        if n == name:
            if sOrT == "s":
                source = _repeat(lambda x: x.source, cell, depth)
                if source != match.stInfo[key]:
                    return False
            elif sOrT == "t":
                target = _repeat(lambda x: x.target, cell, depth)
                if target != match.stInfo[key]:
                    return False
    return True

class VarNode(ASTNode):
    # Name is an integer index when used in primitive family source/target ASTs
    def __init__(self, name):
        self.name = name
        self.size = 1

    def match(self, target, matchList):
        updatedMatches = []
        for match in matchList:
            if _isCompatibleSetCell(self.name, target, match):
                new = match.dup()
                new.varMatch[self.name] = target
                updatedMatches.append(new)
        return updatedMatches

    def eval(self, assignment, functors):
        return assignment[self.name]

    def alterAST(self, func):
        return self

    def __str__(self):
        return str(self.name)

# NOTE: For solvable families, source and target nodes can only be applied to VarNodes and other source/target nodes.
class SourceNode(ASTNode):
    def __init__(self, subNode):
        if isinstance(subNode, VarNode):
            self.depth = 1
            self.varName = subNode.name
        elif isinstance(subNode, SourceNode) or isinstance(subNode, TargetNode):
            if subNode.depth is not None:
                self.depth = 1 + subNode.depth
            self.varName = subNode.varName
        else:
            # Not solvable; used for eval-only trees.
            self.depth = None
            self.varName = None
        self.subNode = subNode
        self.size = 1 + subNode.size


    def match(self, target, matchList):
        if self.varName is None:
            raise Exception("Can not match on SourceNode; does not descend to VarNode.")
        updatedMatches = []
        for match in matchList:
            if _isCompatibleSetSourceTarget(self.varName, "s", self.depth, target, match):
                new = match.dup()
                new.stInfo[(self.varName, "s", self.depth)] = target
                updatedMatches.append(new)
        return updatedMatches

    def eval(self, assignment, functors):
        return self.subNode.eval(assignment, functors).source

    def alterAST(self, func):
        res = func(self.subNode)
        if res is None:
            return SourceNode(self.subNode.alterAST(func))
        else:
            return SourceNode(res)

    def __str__(self):
        return "#s(" + str(self.subNode) + ")"

class TargetNode(ASTNode):
    def __init__(self, subNode):
        if isinstance(subNode, VarNode):
            self.depth = 1
            self.varName = subNode.name
        elif isinstance(subNode, SourceNode) or isinstance(subNode, TargetNode):
            if subNode.depth is not None:
                self.depth = 1 + subNode.depth
            self.varName = subNode.varName
        else:
            # Not solvable; used for eval-only trees.
            self.depth = None
            self.varName = None
        self.subNode = subNode
        self.size = 1 + subNode.size


    def match(self, target, matchList):
        if self.varName is None:
            raise Exception("Can not match on TargetNode; does not descend to VarNode.")
        updatedMatches = []
        for match in matchList:
            if _isCompatibleSetSourceTarget(self.varName, "t", self.depth, target, match):
                new = match.dup()
                new.stInfo[(self.varName, "t", self.depth)] = target
                updatedMatches.append(new)
        return updatedMatches

    def eval(self, assignment, functors):
        return self.subNode.eval(assignment, functors).target

    def alterAST(self, func):
        res = func(self.subNode)
        if res is None:
            return TargetNode(self.subNode.alterAST(func))
        else:
            return TargetNode(res)

    def __str__(self):
        return "#t(" + str(self.subNode) + ")"

# Here functor should be an integer index starting from 0, indicating which of the functor parameters this is on.
class FunctorNode(ASTNode):
    def __init__(self, functor, subNode):
        self.functor = functor
        self.subNode = subNode
        self.size = 1 + subNode.size

    def match(self, target, matchList):
        if dim(target) == 0:
            prim = asPrim0(target)
        elif dim(target) == 1:
            prim = asPrim1(target)
        elif dim(target) == 2:
            prim = asPrim2(target)
        if isinstance(prim, FunctorPrim0):
            res = self.subNode.match(prim.mol0, matchList)
        elif isinstance(prim, FunctorPrim1):
            res = self.subNode.match(prim.eqMol1, matchList)
        elif isinstance(prim, FunctorPrim2):
            res = self.subNode.match(prim.eqAEMol2, matchList)
        else:
            return []

        f = prim.functor
        updatedMatches = []
        for match in res:
            if self.functor in match.functorMatches:
                if match.functorMatches[self.functor] == prim.functor:
                    updatedMatches.append(match)
            else:
                newMatch = match.dup()
                newMatch.functorMatches[self.functor] = f
                updatedMatches.append(newMatch)
        return updatedMatches

    def eval(self, assignment, functors):
        cell = self.subNode.eval(assignment, functors)
        if isinstance(functors, dict) and self.functor in functors:
            functor = functors[self.functor]
        elif isinstance(functors, list):
            functor = functors[self.functor]
        else:
            functor = self.functor
        if dim(cell) == 0:
            return fFunctorPrim0(functor, ensureMol0(cell))
        elif dim(cell) == 1:
            return fFunctorPrim1(functor, ensureEqMol1(cell))
        elif dim(cell) == 2:
            return fFunctorPrim2(functor, ensureEqAEMol2(cell))

    def alterAST(self, func):
        res = func(self.subNode)
        if res is None:
            return FunctorNode(self.functor, self.subNode.alterAST(func))
        else:
            return FunctorNode(self.functor, res)

    def __str__(self):
        return "F_{" + str(self.functor) + "}(" +  str(self.subNode) + ")"

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

    def eval(self, assignment, functors):
        return self.matchTarget

    def alterAST(self, func):
        return self

    def __str__(self):
        return str(self.prim)

def alterAST(func, ast):
    res = func(ast)
    if res is None:
        return ast.alterAST(func)
    else:
        return res

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

# class FamilyPrim0(Prim0):
#     def __init__(self, family, params, functors):
#         self.family = family
#         self.params = params
#         self.functors = functors
#
#     def __str__(self):
#         return self.family.strParams(self.params, self.functors)

class FamilyPrim1(Prim1):
    def __init__(self, family, params, functors):
        self.family = family
        self.params = params
        self.functors = functors
        self.source = family.sourceAST.eval(listToIndexDict(params), functors)
        self.target = family.targetAST.eval(listToIndexDict(params), functors)

    def __str__(self):
        return self.family.strParams(self.params, self.functors)

class FamilyPrim2(Prim2):
    def __init__(self, family, params, functors):
        self.family = family
        self.params = params
        self.functors = functors
        self.source = family.sourceAST.eval(listToIndexDict(params), functors)
        self.target = family.targetAST.eval(listToIndexDict(params), functors)
        # Atom2s won't have well-defined source/target if we don't select an instance.
        if isinstance(self.source, EqMol1):
            self.source = next(iter(self.source.mol1s))
        if isinstance(self.target, EqMol1):
            self.target = next(iter(self.target.mol1s))
        super(FamilyPrim2, self).__init__()


    def __str__(self):
        return self.family.strParams(self.params, self.functors)

    def __repr__(self):
        return "FamilyPrim2(" + self.family.name + ", " + repr(self.params) + ")"

class FamilyPrim3(Prim3):
    def __init__(self, family, params, functors):
        self.family = family
        self.params = params
        self.functors = functors
        self.source = self.family.sourceAST.eval(listToIndexDict(self.params), functors)
        self.target = self.family.targetAST.eval(listToIndexDict(self.params), functors)

    def __str__(self):
        return self.family.strParams(self.params, self.functors)

class AdjIdPrim2(Prim2):
    # NOTE: Don't build directly, use factory functions below.
    def __init__(self, left, right):
        if not(isinstance(self.left, FamilyPrim1) and isinstance(self.right, FamilyPrim1) and self.left.family == self.right.family.adj) or (self.left.params != self.right.params):
            raise Exception("FamilyPrim1s are not adjoint.")

        self.left = left
        self.right = right
        self.source = comp1(left, right)
        self.target = fIdMol1(left.source)
        super(AdjIdPrim2, self).__init__()

    def __str__(self):
        return "adjid(" + str(self.left) + ", " + str(self.right) + ")"

    def __repr__(self):
        return "AdjIdPrim2(" + repr(self.left) + ", " + repr(self.right) + ")"

# TODO: Mates

class AdjIdPrim3(Prim3):
    # NOTE: Don't build directly, use factory functions below.
    def __init__(self, left, right):
        if not(isinstance(self.left, FamilyPrim2) and isinstance(self.right, FamilyPrim2) and self.left.family == self.right.family.adj) or (self.left.params != self.right.params):
            raise Exception("FamilyPrim2s are not adjoint.")

        self.left = left
        self.right = right
        self.source = comp2(left, right)
        self.target = fAEIdMol2(left.source)
        super(AdjIdPrim3, self).__init__()

    def __str__(self):
        return "adjid(" + str(self.left) + ", " + str(self.right) + ")"

    def __repr__(self):
        return "AdjIdPrim3(" + repr(self.left) + ", " + repr(self.right) + ")"

# Keys are (params, family1, family2)
adjIdPrim2Repo = {}
adjIdPrim3Repo = {}

def fAdjIdPrim2(left, right):
    if left.params != right.params:
        raise Exception("Parameters do not match.")
    if left.family != right.family.adj:
        raise Exception("Prims are not adjoint.")
    famPair = (left.params, left.family, right.family)
    if famPair in adjIdPrim2Repo:
        return adjIdPrim2Repo[famPair]
    else:
        res = AdjIdPrim2(left, right)
        adjIdPrim2Repo[famPair] = res
        return res

def fAdjIdPrim3(left, right):
    if left.params != right.params:
        raise Exception("Parameters do not match.")
    if left.family != right.family.adj:
        raise Exception("Prims are not adjoint.")
    famPair = (left.params, left.family, right.family)
    if famPair in adjIdPrim3Repo:
        return adjIdPrim2Repo[famPair]
    else:
        res = AdjIdPrim3(left, right)
        adjIdPrim3Repo[famPair] = res
        return res

def adjointDim2AST(ast):
    if isinstance(ast, PrimitiveFamilyNode):
        if ast.family.dim == 2:
            return PrimitiveFamilyNode(ast.family.adj, ast.paramNodes)
        else:
            return ast
    elif isinstance(ast, FunctorNode):
        return FunctorNode(ast.functor, adjointDim2AST(ast.subNode))
    elif isinstance(ast, Comp2Node):
        return Comp2Node(adjointDim2AST(ast.rightNode), adjointDim2AST(ast.leftNode))
    elif isinstance(ast, Comp1Node):
        return Comp1Node(adjointDim2AST(ast.leftNode), adjointDim2AST(ast.rightNode))
    elif isinstance(ast, Comp0Node):
        return Comp0Node(adjointDim2AST(ast.leftNode), adjointDim2AST(ast.rightNode))
    elif isinstance(ast, ConstNode) or isinstance(ast, SourceNode) or isinstance(ast, TargetNode) or isinstance(ast, VarNode) or isinstance(ast, IdNode):
        # TODO: Not sure about VarNode here, double check this.
        return ast
    else:
        raise Exception("Adjoint of unknown node type %s" % (type(ast),))

class PrimitiveFamily(IdHashed):
    # Signature is a list of dimensions, one for each parameter slot.  E.g. [1,1,1] takes three arguments, all of dimension 1.
    # AST VarNodes are to be named after the index of each parameter, as an integer.
    def __init__(self, name, dim, functorCount, signature, sourceAST, targetAST, adj=None):
        self.name = name
        self.dim = dim
        self.signature = signature
        self.functorCount = functorCount
        self.sourceAST = sourceAST
        self.targetAST = targetAST
        if adj is None:
            self.adj = PrimitiveFamily(name + "_adj", dim, functorCount, signature, targetAST, sourceAST, self)
        else:
            self.adj = adj

        if self.adj is None:
            print("XXXXXX")

        # TODO: May not need adjId/mates/triangles in dimension 1?  Can we frame all dimension 1 diagrams with no adjoints?

#        if dim == 1 or dim == 2:
        if dim == 2:
            varNodes = [VarNode(ii) for ii in range(len(signature))]
            # Need higher dimensional cells witnessing the fact that composing with its adjoint gives identity.
            # NOTE: These are only source-solvable, not target solvable.  Adjoints do come up in mates.
            # adjId1 corresponds to counit, adjId2 corresponds to inverse/adjoint of unit
            adjId1SourceAST = Comp1Node(PrimitiveFamilyNode(self, varNodes), PrimitiveFamilyNode(self.adj, varNodes))
            adjId1TargetAST = IdNode(SourceNode(PrimitiveFamilyNode(self, varNodes)))
            adjId2SourceAST = Comp1Node(PrimitiveFamilyNode(self.adj, varNodes), PrimitiveFamilyNode(self, varNodes))
            adjId2TargetAST = IdNode(TargetNode(PrimitiveFamilyNode(self, varNodes)))
            self.adjId1 = PrimitiveFamily(name + "_adjid1", dim + 1, signature, functorCount, adjId1SourceAST, adjId1TargetAST)
            self.adjId2 = PrimitiveFamily(name + "_adjid2", dim + 1, signature, functorCount, adjId2SourceAST, adjId2TargetAST)

        # if dim == 1:
        #     # Also need higher dimensional cells witnessing triangle axioms.
        #     # NOTE: See notes; we can format the triangle axioms as solvable prims on single atoms.  These prims imply (as witnesses) the original axioms.
        #     tri1SourceAST = Comp1Node(PrimitiveFamilyNode(self.adjId1, varNodes), PrimitiveFamilyNode(self, varNodes))
        #     tri1TargetAST = Comp1Node(PrimitiveFamilyNode(self, varNodes), PrimitiveFamilyNode(self.adjId2, varNodes))
        #     tri2SourceAST = Comp1Node(PrimitiveFamilyNode(self.adj, varNodes), PrimitiveFamilyNode(self.adjId1, varNodes))
        #     tri2TargetAST = Comp1Node(PrimitiveFamilyNode(self.adjId2, varNodes), PrimitiveFamilyNode(self.adj, varNodes))
        #     self.tri1 = PrimitiveFamily(name + "_tri1", dim + 2, signature, tri1SourceAST, tri1TargetAST)
        #     self.tri2 = PrimitiveFamily(name + "_tri2", dim + 2, signature, tri2SourceAST, tri2TargetAST)

        # TODO: Mate ASTs for matching; unsure if this is needed.  If all dim1 diagrams are adjoint-free, then this is relevant for dim2 diagrams if they require matching under adjoints.
        # (instead of just cancellation of adjoints).
        # Instead of relating mate to original prims, since we're only doing this in dim3, we just make a new prim family.  Reduce externally when necessary.
        # NOTE: Had to move this elsewhere, since it caused recursion problems.
        # if dim == 3:
        #     self.mate = PrimitiveFamily(name + "_mate", dim, signature, functorCount, adjointDim2AST(self.sourceAST), adjointDim2AST(self.targetAST))

        self.primByParams = {}
        # Hook for filtering degenerate instances.  Used by filling, write a hook for families that have degenerate or expanding primitive instances.
        # TODO: Don't use hook if we're already limiting instances of the prim.
        self.isDegen = lambda params: False


    def getMate(self):
        if hasattr(self, "mate"):
            return self.mate
        else:
            if self.dim == 3:
                self.mate = PrimitiveFamily(self.name + "_mate", self.dim, self.signature, self.functorCount, adjointDim2AST(self.sourceAST), adjointDim2AST(self.targetAST))
                self.mate.isDegen = self.isDegen
                return self.mate
            else:
                raise Exception("Mate of Primitive Families not of dimension 3 not supported.")


    def fprim(self, *params):
        return self.buildPrim(params, [])

    def fprimf(self, functors, *params):
        return self.buildPrim(params, functors)

    def buildPrim(self, params, functors):
        # Need as an ordered list.
        if isinstance(params, dict):
            params = indexDictToList(params)
        if isinstance(functors, dict):
            functors = indexDictToList(functors)

        functors = take(functors, self.functorCount)

        if len(params) == len(self.signature):
            for ii in range(len(params)):
                if dim(params[ii]) != self.signature[ii]:
                    raise Exception("Signature error")

        # Make sure parameters are in equivalence form.
        adjParams = []
        for p in params:
            if dim(p) == 0:
                adjParams.append(ensureMol0(p))
            elif dim(p) == 1:
                adjParams.append(ensureEqMol1(p))
            elif dim(p) == 2:
                adjParams.append(ensureEqAEMol2(p))
        params = adjParams

        parTup = (tuple(params), tuple(functors))
        if parTup  in self.primByParams:
            x = self.primByParams[parTup]
            return self.primByParams[parTup]
        else:
            if self.dim == 0:
                prim = FamilyPrim0(self, params, functors)
            elif self.dim == 1:
                prim = FamilyPrim1(self, params, functors)
            elif self.dim == 2:
                prim = FamilyPrim2(self, params, functors)
            elif self.dim == 3:
                prim = FamilyPrim3(self, params, functors)
            self.primByParams[parTup] = prim
            return prim

    def strParams(self, params, functors):
        if functors is None or len(functors) == 0:
            return self.name + "(" + ", ".join(map(str, params)) + ")"
        else:
            return self.name + "{" + ", ".join(map(str, functors)) + "}(" + ", ".join(map(str, params)) + ")"

    def __str__(self):
        return self.name

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

def minimalASTFromMol0(mol0, paramPrims, functors = []):
    if len(mol0) == 0:
        raise Exception("No AST for empty mol0")
    elif len(mol0) == 1:
        prim0 = mol0.prim0s[0]
        if prim0 in paramPrims:
            name = paramPrims.index(prim0)
            return VarNode(name)
        elif isinstance(prim0, FunctorPrim0) and prim0.functor in functors:
            fi = functors.index(prim0.functor)
            return FunctorNode(fi, minimalASTFromMol0(prim0.mol0, paramPrims, functors))
        else:
            return ConstNode(prim0)
    tensorDecomps = tensorDecompMol0(mol0)
    candidates = [Comp0Node(minimalASTFromMol0(l, paramPrims, functors), minimalASTFromMol0(r, paramPrims, functors)) for (l, r) in tensorDecomps]
    return min(candidates, key=lambda x: x.size)

def minimalASTFromEqMol1(eqMol1, paramPrims, functors = []):
    if len(eqMol1) == 0:
        idMol1 = next(iter(eqMol1.mol1s))
        return IdNode(minimalASTFromMol0(idMol1.mol0, paramPrims, functors))
    elif len(eqMol1) == 1:
        atom1 = next(iter(eqMol1.mol1s)).atom1s[0]
        if len(atom1.a0) == 0 and len(atom1.b0) == 0:
            if atom1.p1 in paramPrims:
                name = paramPrims.index(atom1.p1)
                return VarNode(name)
            elif isinstance(atom1.p1, FunctorPrim1) and atom1.p1.functor in functors:
                fi = functors.index(atom1.p1.functor)
                return FunctorNode(fi, minimalASTFromEqMol1(atom1.p1.eqMol1, paramPrims, functors))
            elif isinstance(atom1.p1, FamilyPrim1):
                paramNodes = []
                for ii in range(len(atom1.p1.family.signature)):
                    dim = atom1.p1.family.signature[ii]
                    par = atom1.p1.params[ii]
                    if dim == 0:
                        paramNodes.append(minimalASTFromMol0(par, paramPrims, functors))
                    elif dim == 1:
                        paramNodes.append(minimalASTFromEqMol1(par, paramPrims, functors))
                    elif dim == 2:
                        paramNodes.append(minimalASTFromEqAEMol2(par, paramPrims, functors))
                return PrimitiveFamilyNode(atom1.p1.family, paramNodes)
            else:
                return ConstNode(atom1.p1)
    tensorDecomps = tensorDecompEqMol1(eqMol1)
    horizontalDecomps = horizontalDecompEqMol1(eqMol1)

    candidates = []
    candidates.extend([Comp0Node(minimalASTFromEqMol1(l, paramPrims, functors), minimalASTFromEqMol1(r, paramPrims, functors)) for (l,r) in tensorDecomps])
    candidates.extend([Comp1Node(minimalASTFromEqMol1(l, paramPrims, functors), minimalASTFromEqMol1(r, paramPrims, functors)) for (l,r) in horizontalDecomps])
    return min(candidates, key=lambda x: x.size)


# Finds all possible ASTs for the AEMol2, then takes the minimal one.
def minimalASTFromEqAEMol2(eqAEMol2, paramPrims, functors = []):
    if len(eqAEMol2) == 0:
        aeIdMol2 = next(iter(eqAEMol2.aeMol2s))
        return IdNode(minimalASTFromEqMol1(aeIdMol2.eqMol1, paramPrims, functors))
    elif len(eqAEMol2) == 1:
        atom2 = next(iter(next(iter(eqAEMol2.aeMol2s)).eqAtom2s[0].atom2s))
        if len(atom2.a0) == 0 and len(atom2.b0) == 0 and len(atom2.l1) == 0 and len(atom2.r1) == 0:
            if atom2.p2 in paramPrims:
                name = paramPrims.index(atom2.p2)
                return VarNode(name)
            elif isinstance(atom2.p2, FunctorPrim2) and atom2.p2.functor in functors:
                fi = functors.index(atom2.p2.functor)
                return FunctorNode(fi, minimalASTFromEqAEMol2(atom2.p2.eqAEMol2, paramPrims, functors))
            elif isinstance(atom2.p2, FamilyPrim2):
                paramNodes = []
                for ii in range(len(atom2.p2.family.signature)):
                    dim = atom2.p2.family.signature[ii]
                    par = atom2.p2.params[ii]
                    if dim == 0:
                        paramNodes.append(minimalASTFromMol0(par, paramPrims, functors))
                    elif dim == 1:
                        paramNodes.append(minimalASTFromEqMol1(par, paramPrims, functors))
                    elif dim == 2:
                        paramNodes.append(minimalASTFromEqAEMol2(par, paramPrims, functors))
                return PrimitiveFamilyNode(atom2.p2.family, paramNodes)
            else:
                return ConstNode(atom2.p2)
    tensorDecomps = tensorDecompEqAEMol2(eqAEMol2)
    horizontalDecomps = horizontalDecompEqAEMol2(eqAEMol2)
    verticalDecomps = verticalDecompEqAEMol2(eqAEMol2)

    candidates = []
    candidates.extend([Comp0Node(minimalASTFromEqAEMol2(l, paramPrims, functors), minimalASTFromEqAEMol2(r, paramPrims, functors)) for (l,r) in tensorDecomps])
    candidates.extend([Comp1Node(minimalASTFromEqAEMol2(l, paramPrims, functors), minimalASTFromEqAEMol2(r, paramPrims, functors)) for (l,r) in horizontalDecomps])
    candidates.extend([Comp2Node(minimalASTFromEqAEMol2(l, paramPrims, functors), minimalASTFromEqAEMol2(r, paramPrims, functors)) for (l,r) in verticalDecomps])
    return min(candidates, key=lambda x: x.size)


# Replaces the VarNodes by index with the nodes in newVarNodes.
def replaceVarNodes(newVarNodes, ast):
    ast.alterAST(lambda x: _replaceVarNodes(newVarNodes, x))


def _replaceVarNodes(newVarNodes, ast):
    if isinstance(ast, VarNode):
        return newVarNodes[ast.name]
    return None

# Returns an AST which is the source of a given AST.
# We do this instead of just sticking under SourceNode to preserve solvability.
# TODO: Make a note that we could solve something like SourceNode(ConstNode), but because we can do this we don't allow it for simplification.
def applySourceAST(ast):
    if isinstance(ast, VarNode) or isinstance(ast, SourceNode) or isinstance(ast, TargetNode):
        return SourceNode(ast)
    elif isinstance(ast, PrimitiveFamilyNode):
        return ast.family.sourceAST.alterAST(lambda x: _replaceVarNodes(ast.paramNodes, x))
    elif isinstance(ast, IdNode):
        return ast.subNode
    elif isinstance(ast, ConstNode):
        return ConstNode(ast.prim.source)
    elif isinstance(ast, FunctorNode):
        return FunctorNode(ast.functor, applySourceAST(ast.subNode))
    elif isinstance(ast, Comp0Node):
        return Comp0Node(applySourceAST(ast.leftNode), applySourceAST(ast.rightNode))
    elif isinstance(ast, Comp1Node):
        return Comp1Node(applySourceAST(ast.leftNode), applySourceAST(ast.rightNode))
    elif isinstance(ast, Comp2Node):
        return Comp2Node(applySourceAST(ast.leftNode), applySourceAST(ast.rightNode))
    raise Exception("Taking source of unknown node type: %s" % (type(ast),))

# Returns an AST which is the target of a given AST.
def applyTargetAST(ast):
    if isinstance(ast, VarNode) or isinstance(ast, SourceNode) or isinstance(ast, TargetNode):
        return TargetNode(ast)
    elif isinstance(ast, PrimitiveFamilyNode):
        return ast.family.targetAST.alterAST(lambda x: _replaceVarNodes(ast.paramNodes, x))
    elif isinstance(ast, IdNode):
        return ast.subNode
    elif isinstance(ast, ConstNode):
        return ConstNode(ast.prim.target)
    elif isinstance(ast, FunctorNode):
        return FunctorNode(ast.functor, applyTargetAST(ast.subNode))
    elif isinstance(ast, Comp0Node):
        return Comp0Node(applyTargetAST(ast.leftNode), applyTargetAST(ast.rightNode))
    elif isinstance(ast, Comp1Node):
        return Comp1Node(applyTargetAST(ast.leftNode), applyTargetAST(ast.rightNode))
    elif isinstance(ast, Comp2Node):
        return Comp2Node(applyTargetAST(ast.leftNode), applyTargetAST(ast.rightNode))
    raise Exception("Taking target of unknown node type: %s" % (type(ast),))
