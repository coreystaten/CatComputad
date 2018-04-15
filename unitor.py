from fill import *
from ontology import *
from primfamily import *
from transfor import *

unit = ConstPrim0("I")
_a = ConstPrim0("a")
_b = ConstPrim0("b")
_c = ConstPrim0("c")

unitor10Source = minimalASTFromMol0(comp0(unit, _a), [_a])
unitor10Target = minimalASTFromMol0(prim0ToMol0(_a), [_a])
unitor10 = PrimitiveFamily("unitor10", 1, 0, [0], unitor10Source, unitor10Target)
unitor10Dim2 = tritrans2CellPrimFamily(unitor10)
unitor10Dim2.adj.isDegen = lambda params: isIdMol1(params[0])
unitor10Dim3 = tritrans3CellPrimFamily(unitor10, unitor10Dim2)
unitor10Pi = tritransPiPrimFamily(unitor10, unitor10Dim2)
unitor10M = tritransMPrimFamily(unitor10, unitor10Dim2)

unitor11Source = minimalASTFromMol0(comp0(_a, unit), [_a])
unitor11Target = minimalASTFromMol0(prim0ToMol0(_a), [_a])
unitor11 = PrimitiveFamily("unitor11", 1, 0, [0], unitor11Source, unitor11Target)
unitor11Dim2 = tritrans2CellPrimFamily(unitor11)
unitor11Dim3 = tritrans3CellPrimFamily(unitor11, unitor11Dim2)
unitor11Pi = tritransPiPrimFamily(unitor11, unitor11Dim2)
unitor11M = tritransMPrimFamily(unitor11, unitor11Dim2)

unitorDim2ByDim1 = {unitor10: unitor10Dim2, unitor11: unitor11Dim2}

#unitor20Paths1 = findPaths1(comp0s(unit, _a, _b), comp0s(_a,_b), [unitor10, unitor11])
unitor20Source = minimalASTFromEqMol1(fEqMol1(comp0(prim1ToMol1(unitor10.fprim(_a)),_b)), [_a,_b])
unitor20Target = minimalASTFromEqMol1(fEqMol1(prim1ToMol1(unitor10.fprim(comp0(_a, _b)))), [_a,_b])
unitor20 = PrimitiveFamily("unitor20", 2, 0, [0, 0], unitor20Source, unitor20Target)
unitor20Dim3 = trimod3CellPrimFamily(unitor20, unitorDim2ByDim1)

#unitor21Paths1 = findPaths1(comp0s([_a, unit, _b]), comp0s([_a,_b]), [unitor10, unitor11])
unitor21Source = minimalASTFromEqMol1(fEqMol1(
    comp0(prim1ToMol1(unitor11.fprim(_a)), _b)), [_a, _b])
unitor21Target = minimalASTFromEqMol1(fEqMol1(
    comp0(_a, prim1ToMol1(unitor10.fprim(_b)))), [_a, _b])
unitor21 = PrimitiveFamily("unitor21", 2, 0, [0, 0], unitor21Source, unitor21Target)
unitor21Dim3 = trimod3CellPrimFamily(unitor21, unitorDim2ByDim1)

#unitor22Paths1 = findPaths1(comp0s([_a, _b, unit]), comp0s([_a, _b]), [unitor10, unitor11])
unitor22Source = minimalASTFromEqMol1(fEqMol1(
    prim1ToMol1(unitor11.fprim(comp0(_a, _b)))), [_a, _b])
unitor22Target = minimalASTFromEqMol1(fEqMol1(
    comp0(_a, prim1ToMol1(unitor11.fprim(_b)))), [_a, _b])
unitor22 = PrimitiveFamily("unitor22", 2, 0, [0, 0], unitor22Source, unitor22Target)
unitor22Dim3 = trimod3CellPrimFamily(unitor22, unitorDim2ByDim1)

#unitor30Paths1 = findPaths1(comp0s([unit, _a, _b, _c]), comp0s([_a, _b, _c]), [unitor10, unitor11])
# searchForPathPairs2AdjFam(unitor30Paths1, [unitor20, unitor21, unitor22], [unitor10Dim2, unitor11Dim2, unitor10Dim2.adj, unitor11Dim2.adj])
#unitor30Source1 = ensureEqMol1(comp0s([prim1ToMol1(unitor10.fprim(_a)), _b, _c]))
#unitor30Target1 = ensureEqMol1(unitor10.fprim(comp0s([_a, _b, _c])))
#unitor30Paths2 = findPaths2(unitor30Source1, unitor30Target1, [unitor20])
unitor30Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    unitor20.fprim(_a, comp0(_b, _c))), [_a, _b, _c])
unitor30Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        comp0(unitor20.fprim(_a, _b), _c),
        unitor20.fprim(comp0(_a, _b), _c))), [_a, _b, _c])
unitor30 = PrimitiveFamily("unitor30", 3, 0, [0, 0, 0], unitor30Source, unitor30Target)

#unitor31Paths1 = findPaths1(comp0s([_a, unit, _b, _c]), comp0s([_a, _b, _c]), [unitor10, unitor11])
#searchForPathPairs2AdjFam(unitor31Paths1, [unitor20, unitor21, unitor22], [unitor10Dim2, unitor11Dim2, unitor10Dim2.adj, unitor11Dim2.adj])
#unitor31Source1 = ensureEqMol1(comp0s([unitor11.fprim(_a), _b, _c]))
#unitor31Target1 = ensureEqMol1(comp0(_a, unitor10.fprim(comp0(_b,_c))))
#unitor31Paths2 = findPaths2(unitor31Source1, unitor31Target1, [unitor20, unitor21])
# [[unitor21(a, (b @ c))]] -- [[(unitor21(a, b) @ c)] & [(a @ unitor20(b, c))]]
unitor31Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    unitor21.fprim(_a, comp0(_b, _c))), [_a, _b, _c])
unitor31Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(comp0(unitor21.fprim(_a,_b), _c), comp0(_a, unitor20.fprim(_b, _c)))), [_a, _b, _c])
unitor31 = PrimitiveFamily("unitor31", 3, 0, [0, 0, 0], unitor31Source, unitor31Target)

#unitor32Paths1 = findPaths1(comp0s([_a, _b, unit, _c]), comp0s([_a, _b, _c]), [unitor10, unitor11])
#searchForPathPairs2AdjFam(unitor32Paths1, [unitor20, unitor21, unitor22], [unitor10Dim2, unitor11Dim2, unitor10Dim2.adj, unitor11Dim2.adj])
# [(unitor11((a @ b)) @ c)] -- [((a @ b) @ unitor10(c))]
#unitor32Source1 = ensureEqMol1(comp0(unitor11.fprim(comp0(_a,_b)), _c))
#unitor32Target1 = ensureEqMol1(comp0s([_a, _b, unitor10.fprim(_c)]))
#unitor32Paths2 = findPaths2(unitor32Source1, unitor32Target1, [unitor21, unitor22])
# [[unitor21((a @ b), c)]] -- [[(unitor22(a, b) @ c)] & [(a @ unitor21(b, c))]]
unitor32Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    unitor21.fprim(comp0(_a,_b), _c)), [_a, _b, _c])
unitor32Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        comp0(unitor22.fprim(_a, _b), _c),
        comp0(_a, unitor21.fprim(_b, _c))
    )), [_a, _b, _c])
unitor32 = PrimitiveFamily("unitor32", 3, 0, [0, 0, 0], unitor32Source, unitor32Target)

# unitor33Paths1 = findPaths1(comp0s([_a, _b, _c, unit]), comp0s([_a, _b, _c]), [unitor10, unitor11])
# searchForPathPairs2AdjFam(unitor33Paths1, [unitor20, unitor21, unitor22], [unitor10Dim2, unitor11Dim2, unitor10Dim2.adj, unitor11Dim2.adj])
# [(unitor11((a @ b @ c)))] -- [((a @ b) @ unitor11(c))]
#unitor33Source1 = ensureEqMol1(unitor11.fprim(comp0s([_a, _b, _c])))
#unitor33Target1 = ensureEqMol1(comp0s([_a, _b, unitor11.fprim(_c)]))
#unitor33Paths2 = findPaths2(unitor33Source1, unitor33Target1, [unitor22])
# [[unitor22((a @ b), c)]] -- [[unitor22(a, (b @ c))] & [(a @ unitor22(b, c))]]
unitor33Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    unitor22.fprim(comp0(_a, _b), _c)), [_a, _b, _c])
unitor33Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        unitor22.fprim(_a, comp0(_b, _c)),
        comp0(_a, unitor22.fprim(_b, _c)))), [_a, _b, _c])
unitor33 = PrimitiveFamily("unitor33", 3, 0, [0,0,0], unitor33Source, unitor33Target)
