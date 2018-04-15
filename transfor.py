from ontology import *
from primfamily import *

# Throughout we follow Gurski's "An algebraic theory of tricategories",
# incluiding choice of variables names where appropriate.



#####################################################
# Tritransformations between functors C^n -> C coming from ASTs; all such functors are strict.
# Input data is the 1-cell primitive family for the tritransformation.
def tritrans2CellSourceAST(famDim1):
    sig = famDim1.signature
    for ii in range(len(sig)):
        if sig[ii] != 0:
            raise Exception("Invalid sig")
    theta_b = PrimitiveFamilyNode(famDim1, [TargetNode(VarNode(ii)) for ii in range(len(sig))])
    # The original source AST also represents the source functor.
    return Comp1Node(famDim1.sourceAST, theta_b)

def tritrans2CellTargetAST(famDim1):
    sig = famDim1.signature
    for ii in range(len(sig)):
        if sig[ii] != 0:
            raise Exception("Invalid sig")
    theta_a = PrimitiveFamilyNode(famDim1, [SourceNode(VarNode(ii)) for ii in range(len(sig))])
    return Comp1Node(theta_a, famDim1.targetAST)

def tritrans2CellPrimFamily(famDim1):
    sigLen = len(famDim1.signature)
    return PrimitiveFamily(famDim1.name + "_dim2", 2, famDim1.functorCount, [1 for s in range(sigLen)], tritrans2CellSourceAST(famDim1), tritrans2CellTargetAST(famDim1))

# Note: when ascending 1-cell functors, don't ascend the parts within source/targets.  Thus will need to replace source nodes with double source nodes, etc.

def tritrans3CellSourceAST(famDim1, famDim2):
     sigLen = len(famDim1.signature)
     F_alpha = famDim1.sourceAST
     idTheta_b = IdNode(PrimitiveFamilyNode(famDim1, [TargetNode(TargetNode(VarNode(ii))) for ii in range(sigLen)]))
     theta_g = PrimitiveFamilyNode(famDim2, [TargetNode(VarNode(ii)) for ii in range(sigLen)])
     return Comp2Node(Comp1Node(F_alpha, idTheta_b), theta_g)

def tritrans3CellTargetAST(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    G_alpha = famDim1.targetAST
    idTheta_a = IdNode(PrimitiveFamilyNode(famDim1, [SourceNode(SourceNode(VarNode(ii))) for ii in range(sigLen)]))
    theta_f = PrimitiveFamilyNode(famDim2, [SourceNode(VarNode(ii)) for ii in range(sigLen)])
    return Comp2Node(theta_f, Comp1Node(idTheta_a, G_alpha))

def tritrans3CellPrimFamily(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamily(famDim1.name + "_dim3", 3, famDim1.functorCount, [2 for s in range(sigLen)], tritrans3CellSourceAST(famDim1, famDim2), tritrans3CellTargetAST(famDim1, famDim2))

# alterAST function which add inc to all the VarNode names.
def addVarIndex(x, inc):
    if isinstance(x, VarNode):
        return VarNode(x.name + inc)
    return None


def tritransPiSourceAST(famDim1, famDim2):
    # Pi takes as arguments composable (f_i) and (g_i), i.e. t(f_i) = s(g_i).
    # Use the first sigLen parameters for f_i, and the next sigLen parameters for g_i.  Note that (... + sigLen) in several locations for this.
    sigLen = len(famDim1.signature)
    idF_f = IdNode(famDim1.sourceAST)
    theta_g = PrimitiveFamilyNode(famDim2, [VarNode(ii + sigLen) for ii in range(sigLen)])
    theta_f = PrimitiveFamilyNode(famDim2, [VarNode(ii) for ii in range(sigLen)])
    idG_g = IdNode(famDim1.targetAST.alterAST(lambda x: addVarIndex(x, sigLen)))
    return Comp2Node(Comp1Node(idF_f, theta_g), Comp1Node(theta_f, idG_g))

def tritransPiTargetAST(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamilyNode(famDim2, [Comp1Node(VarNode(ii), VarNode(ii + sigLen)) for ii in range(sigLen)])

def tritransPiPrimFamily(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamily(famDim1.name + "_pi", 3, famDim1.functorCount, [1 for s in range(2 * sigLen)], tritransPiSourceAST(famDim1, famDim2), tritransPiTargetAST(famDim1, famDim2))

def tritransMSourceAST(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamilyNode(famDim2, [IdNode(VarNode(ii)) for ii in range(sigLen)])

def tritransMTargetAST(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return IdNode(PrimitiveFamilyNode(famDim1, [VarNode(ii) for ii in range(sigLen)]))

def tritransMPrimFamily(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamily(famDim1.name + "_m", 3, famDim1.functorCount, [0 for s in range(sigLen)], tritransMSourceAST(famDim1, famDim2), tritransMTargetAST(famDim1, famDim2))

#####################################################
# Trimodification between transformations on functors C^n -> C coming from ASTs; all such functors are strict.

# Takes an AST representing a tritrans dim1 cell, and returns the AST representing its dim2 cell.
# dim1ByDim2 is a dictionary from 1-cell families to 2-cell families.
def dim1TritransASTToDim2TritransAST(dim1ByDim2, node):
    return _dim1TritransASTToDim2TritransAST(dim1ByDim2, node)[1]

# Returns, at each stage, a tuple of (topAST, newNode, bottomAST)
# Where topAST and bottomAST are ASTs for the naturality 1-cells on each side of the naturality diagram, and newNode is an AST for the new 2-cell so far.
def _dim1TritransASTToDim2TritransAST(dim1ByDim2, node):
    if isinstance(node, PrimitiveFamilyNode):
        # For now, just fail if we encounter a primitive family not in the list.
        # Should be ok for primitive families contained in other 1-cell primitive families, since those won't be recursed on.
        dim2 = dim1ByDim2[node.family]
        return (PrimitiveFamilyNode(node.family, list(map(applySourceAST, node.paramNodes))),
                PrimitiveFamilyNode(dim2, node.paramNodes),
                PrimitiveFamilyNode(node.family, list(map(applyTargetAST, node.paramNodes))))
    elif isinstance(node, VarNode) or isinstance(node, ConstNode) or isinstance(node, IdNode) or isinstance(node, SourceNode) or isinstance(node, TargetNode):
        return (node, IdNode(node), node)
    elif isinstance(node, Comp1Node):   # No Comp2s, since these are 1-cells.
        (leftTop, leftPNode, leftBot) = _dim1TritransASTToDim2TritransAST(dim1ByDim2, node.leftNode)
        (rightTop, rightPNode, rightBot) = _dim1TritransASTToDim2TritransAST(dim1ByDim2, node.rightNode)
        return (
            Comp1Node(leftTop, rightTop),
            Comp2Node(
                Comp1Node(leftPNode, rightBot),
                Comp1Node(leftTop, rightPNode)
            ),
            Comp1Node(leftBot, rightBot))
    elif isinstance(node, Comp0Node):
        (leftTop, leftPNode, leftBot) = _dim1TritransASTToDim2TritransAST(dim1ByDim2, node.leftNode)
        (rightTop, rightPNode, rightBot) = _dim1TritransASTToDim2TritransAST(dim1ByDim2, node.rightNode)
        return (
            Comp0Node(leftTop, rightTop),
            Comp0Node(leftPNode, rightPNode),
            Comp0Node(leftBot, rightBot)
        )
    elif isinstance(node, FunctorNode):
        (top, pNode, bot) = _dim1TritransASTToDim2TritransAST(dim1ByDim2, node.subNode)
        return (FunctorNode(node.functor, top), FunctorNode(node.functor, pNode), FunctorNode(node.functor, bot))
    raise Exception("Unknown node type: %s" % (type(node),))

# Source AST for the 3-cell of a modification.
# famDim2 is the 2-cell family for the modification; dim1ByDim2 is a dictionary mapping all dim1 tritransformation families in the source/target of the modification 2-cell to their dim2 families.
def trimod3CellSourceAST(famDim2, dim1ByDim2):
    theta_f = dim1TritransASTToDim2TritransAST(dim1ByDim2, famDim2.sourceAST)
    idG_f = IdNode(applyTargetAST(famDim2.sourceAST))
    m_a = PrimitiveFamilyNode(famDim2, [SourceNode(VarNode(ii)) for ii in range(len(famDim2.signature))])
    return Comp2Node(theta_f, Comp1Node(idG_f, m_a))

def trimod3CellTargetAST(famDim2, dim1ByDim2):
    m_b = PrimitiveFamilyNode(famDim2, [TargetNode(VarNode(ii)) for ii in range(len(famDim2.signature))])
    idF_f = IdNode(applySourceAST(famDim2.targetAST))
    phi_f = dim1TritransASTToDim2TritransAST(dim1ByDim2, famDim2.targetAST)
    return Comp2Node(Comp1Node(m_b, idF_f), phi_f)

def trimod3CellPrimFamily(famDim2, dim1ByDim2):
    return PrimitiveFamily(famDim2.name + "_dim3", 3, famDim2.functorCount, [1 for ii in range(len(famDim2.signature))], trimod3CellSourceAST(famDim2, dim1ByDim2), trimod3CellTargetAST(famDim2, dim1ByDim2))
