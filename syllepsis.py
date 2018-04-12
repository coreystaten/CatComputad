from fill import *
from ontology import *
from primfamily import *
from transfor import *

from shuffle import *

_a = ConstPrim0("a")
_b = ConstPrim0("b")

syllepsis2Source = minimalASTFromEqMol1(ensureEqMol1(comp1(shuffle1.fprim(_a, _b), shuffle1.fprim(_b, _a))), [_a, _b])
syllepsis2Target = minimalASTFromEqMol1(ensureEqMol1(IdMol1(comp0(_a, _b))), [_a, _b])
syllepsis2 = PrimitiveFamily("syllepsis2", 2, [0, 0], syllepsis2Source, syllepsis2Target)
# TODO: Bimod cell

#syllepsis3Source1 = ensureEqMol1(comp1s(shuffle1.fprim(_a, _b), shuffle1.fprim(_b, _a), shuffle1.fprim(_a, _b)))
#syllepsis3Target1 = ensureEqMol1(shuffle1.fprim(_a, _b))
#syllepsis3Paths2 = findPaths2(syllepsis3Source1, syllepsis3Target1, [syllepsis2])
syllepsis3Source = minimalASTFromEqAEMol2(ensureEqAEMol2(comp1(syllepsis2.fprim(_a, _b), shuffle1.fprim(_a, _b))), [_a, _b])
syllepsis3Target = minimalASTFromEqAEMol2(ensureEqAEMol2(comp1(shuffle1.fprim(_a, _b), syllepsis2.fprim(_b, _a))))
syllepsis3 = PrimitiveFamily("syllepsis3", 3, [0, 0], syllepsis3Source, syllepsis3Target)
