from fill import *
from ontology import *
from primfamily import *

unit = ConstPrim0("I")
_a = ConstPrim0("a")
_b = ConstPrim0("b")

unitor10Source0 = minimalASTFromMol0(comp0(unit, _a), [_a])
unitor10Target0 = minimalASTFromMol0(prim0ToMol0(_a), [_a])
unitor10 = PrimitiveFamily("unitor10", 1, [0], unitor10Source0, unitor10Target0)

unitor11Source0 = minimalASTFromMol0(comp0(_a, unit), [_a])
unitor11Target0 = minimalASTFromMol0(prim0ToMol0(_a), [_a])
unitor11 = PrimitiveFamily("unitor11", 1, [0], unitor11Source0, unitor11Target0)



#unitor20Paths1 = findPaths1(comp0s([unit, _a, _b]), comp0s([_a,_b]), [unitor10, unitor11])
unitor20Source0 = minimalASTFromEqMol1(fEqMol1(comp0(prim1ToMol1(unitor10.buildPrim([_a])),_b)), [_a,_b])
unitor20Target0 = minimalASTFromEqMol1(fEqMol1(prim1ToMol1(unitor10.buildPrim(comp0(_a, _b)))), [_a,_b])
unitor20 = PrimitiveFamily("unitor20", 2, [0, 0], unitor20Source0, unitor20Target0)

unitor21Paths1 = findPaths1(comp0s([_a, unit, _b]), comp0s([_a,_b]), [unitor10, unitor11])
