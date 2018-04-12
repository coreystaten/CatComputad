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
    return PrimitiveFamily(famDim1.name + "_dim2", 2, [1 for s in range(sigLen)], tritrans2CellSourceAST(famDim1), tritrans2CellTargetAST(famDim1))

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
    return PrimitiveFamily(famDim1.name + "_dim3", 3, [2 for s in range(sigLen)], tritrans3CellSourceAST(famDim1, famDim2), tritrans3CellTargetAST(famDim1, famDim2))

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
    return PrimitiveFamilyNode(famDim2, [Comp2Node(VarNode(ii), VarNode(ii + sigLen)) for ii in range(sigLen)])

def tritransPiPrimFamily(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamily(famDim1.name + "_pi", 3, [1 for s in range(2 * sigLen)], tritransPiSourceAST(famDim1, famDim2), tritransPiTargetAST(famDim1, famDim2))

def tritransMSourceAST(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamilyNode(famDim2, [IdNode(VarNode(ii)) for ii in range(sigLen)])

def tritransMTargetAST(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return IdNode(PrimitiveFamilyNode(famDim1, [VarNode(ii) for ii in range(sigLen)]))

def tritransMPrimFamily(famDim1, famDim2):
    sigLen = len(famDim1.signature)
    return PrimitiveFamily(famDim1.name + "_m", 3, [0 for s in range(sigLen)], tritransMSourceAST(famDim1, famDim2), tritransMTargetAST(famDim1, famDim2))
