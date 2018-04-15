import cProfile

from ontology import *
from decomp import *

from shuffle import *
from syllepsis import *
from unitor import *

from tricathom import *
from smt import *

a = prim0ToMol0(ConstPrim0("a"))
b = prim0ToMol0(ConstPrim0("b"))
c = prim0ToMol0(ConstPrim0("c"))
d = prim0ToMol0(ConstPrim0("d"))
e = prim0ToMol0(ConstPrim0("e"))
f = prim0ToMol0(ConstPrim0("f"))

x = prim1ToMol1(ConstPrim1("x",a,b))
y = prim1ToMol1(ConstPrim1("y",c,d))
z = prim1ToMol1(ConstPrim1("z",e,f))

x1 = prim1ToMol1(ConstPrim1("x1",a,b))
x2 = prim1ToMol1(ConstPrim1("x2",a,b))
y1 = prim1ToMol1(ConstPrim1("y1",c,d))
y2 = prim1ToMol1(ConstPrim1("y2",c,d))
z1 = prim1ToMol1(ConstPrim1("z1",c,d))
z2 = prim1ToMol1(ConstPrim1("z2",c,d))

h = prim2ToAEMol2(ConstPrim2("h", x1, x2))
k = prim2ToAEMol2(ConstPrim2("k",y1,y2))
l = prim2ToAEMol2(ConstPrim2("l",y1,y2))


aa1 = fEqMol1(comp0s(x,y,z))
#aa2 = tensorDecompEqMol1(aa1)

#cProfile.run("findPaths(comp0s(a, b, c, d), targetEndCond(comp0s(d, c, b, a)), isNotSingleBraid, [shuffle1], cellsAway1, collate1)")


test = comp0s(h,k,l)
p = test.eqAtom2s
p0 = p[0]
p1 = p[1]
test2 = fEqAEMol2(test)

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

#shuffle1_pi_adj([(shuffle1(b, c))], [1_{d}], [1_{(c @ b)}], [1_{d}])

#cProfile.run("searchForPathPairs3AdjFam(permutahedron4Paths2, [shuffle30, shuffle31, shuffle32.adj, permutahedron3.adj, shuffle1Dim3], [])")
