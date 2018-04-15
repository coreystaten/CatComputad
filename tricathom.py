from fill import *
from ontology import *
from primfamily import *
from transfor import *

_a = ConstPrim0("a")
_b = ConstPrim0("b")
_c = ConstPrim0("c")
_d = ConstPrim0("d")

_f = ConstPrim1("f", _a, _b)
_g = ConstPrim1("g", _b, _c)
_h = ConstPrim1("h", _c, _d)

_k = ConstPrim1("k", _a, _b)
_l = ConstPrim1("l", _a, _b)
_m = ConstPrim1("m", _a, _b)

_u = ConstPrim2("u", _k, _l)
_v = ConstPrim2("v", _l, _m)

_H = Functor("H")

#H(f) . H(g) -> H(f . g)
appComp1Source = minimalASTFromEqMol1(ensureEqMol1(
    comp1(app(_H, _f), app(_H, _g))), [_f, _g], [_H])
appComp1Target = minimalASTFromEqMol1(ensureEqMol1(
    app(_H, comp1(_f, _g))), [_f, _g], [_H])
appComp1 = PrimitiveFamily("appComp1", 2, 1, [1, 1], appComp1Source, appComp1Target)
appComp1.adj.isDegen = lambda params: isIdMol1(params[0]) or isIdMol1(params[1])

# TODO: Move this construction to transfor.py?
appComp1Dim3Source = Comp2Node(
    PrimitiveFamilyNode(appComp1, [SourceNode(VarNode(0)), SourceNode(VarNode(1))]),
    FunctorNode(0, Comp1Node(VarNode(0), VarNode(1)))
)
appComp1Dim3Target = Comp2Node(
    Comp1Node(FunctorNode(0, VarNode(0)), FunctorNode(0, VarNode(1))),
    PrimitiveFamilyNode(appComp1, [TargetNode(VarNode(0)), TargetNode(VarNode(1))])
)
appComp1Dim3 = PrimitiveFamily("appComp1Dim3", 3, 1, [2, 2], appComp1Dim3Source, appComp1Dim3Target)

#H(u) & H(g) -> H(u & g)
appComp2Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(app(_H, _u), app(_H, _v))), [_u, _v], [_H])
appComp2Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    app(_H, comp2(_u, _v))), [_u, _v], [_H])
appComp2 = PrimitiveFamily("appComp2", 3, 1, [2, 2], appComp2Source, appComp2Target)

# 1_{H(a)} -> H(1_a)
appIdentity1Source = minimalASTFromEqMol1(ensureEqMol1(
    IdMol1(ensureMol0(app(_H, _a)))), [_a], [_H])
appIdentity1Target = minimalASTFromEqMol1(ensureEqMol1(
    app(_H, IdMol1(ensureMol0(_a)))), [_a], [_H])
appIdentity1 = PrimitiveFamily("appIdentity1", 2, 1, [0], appIdentity1Source, appIdentity1Target)
# No dim3 cells, since no arrows in unit category.

#1_{H(f)} -> H(1_f)
appIdentity2Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    AEIdMol2(ensureEqMol1(app(_H, _f)))), [_f], [_H])
appIdentity2Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    app(_H, AEIdMol2(ensureEqMol1(_f)))), [_f], [_H])
appIdentity2 = PrimitiveFamily("appIdentity2", 3, 1, [1], appIdentity2Source, appIdentity2Target)


appComp1HexagonSource = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        comp1(appComp1.fprimf([_H], _f, _g), app(_H, _h)),
        appComp1.fprimf([_H], comp1(_f, _g), _h))), [_f, _g, _h], [_H])
appComp1HexagonTarget = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        comp1(app(_H, _f), appComp1.fprimf([_H], _g, _h)),
        appComp1.fprimf([_H], _f, comp1(_g, _h)))), [_f, _g, _h], [_H])
appComp1Hexagon = PrimitiveFamily("appComp1Hexagon", 3, 1, [1, 1, 1], appComp1HexagonSource, appComp1HexagonTarget)


trihomGammaSource = Comp2Node(
    Comp1Node(
        IdNode(FunctorNode(0, VarNode(0))),
        PrimitiveFamilyNode(appIdentity1, [TargetNode(VarNode(0))])
    ),
    PrimitiveFamilyNode(appComp1, [VarNode(0), IdNode(TargetNode(VarNode(0)))])
)
trihomGammaTarget = IdNode(FunctorNode(0, VarNode(0)))
trihomGamma = PrimitiveFamily("trihomGamma", 3, 1, [1], trihomGammaSource, trihomGammaTarget)

trihomDeltaSource = IdNode(FunctorNode(0, VarNode(0)))
trihomDeltaTarget = Comp2Node(
    Comp1Node(
        PrimitiveFamilyNode(appIdentity1, [SourceNode(VarNode(0))]),
        IdNode(FunctorNode(0, VarNode(0)))
    ),
    PrimitiveFamilyNode(appComp1, [IdNode(SourceNode(VarNode(0))), VarNode(0)])
)
trihomDelta = PrimitiveFamily("trihomDelta", 3, 1, [1], trihomDeltaSource, trihomDeltaTarget)
