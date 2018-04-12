import cProfile

from ontology import *
from decomp import *

from shuffle import *
from unitor import *

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


aa1 = fEqMol1(comp0s([x,y,z]))
#aa2 = tensorDecompEqMol1(aa1)



test = comp0s([h,k,l])
p = test.eqAtom2s
p0 = p[0]
p1 = p[1]
test2 = fEqAEMol2(test)
