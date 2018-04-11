import cProfile

from ontology import *
from decomp import *
from fill import *
from unitor import *

a = prim0ToMol0(ConstPrim0("a"))
b = prim0ToMol0(ConstPrim0("b"))
c = prim0ToMol0(ConstPrim0("c"))
d = prim0ToMol0(ConstPrim0("d"))
e = prim0ToMol0(ConstPrim0("e"))
f = prim0ToMol0(ConstPrim0("f"))

x = prim1ToMol1(ConstPrim1("x",b,c))
y = prim1ToMol1(ConstPrim1("y",c,d))
z = prim1ToMol1(ConstPrim1("z",e,f))
w = prim1ToMol1(ConstPrim1("w",c,d))

x1 = prim1ToMol1(ConstPrim1("x1",a,b))
x2 = prim1ToMol1(ConstPrim1("x2",a,b))
x3 = prim1ToMol1(ConstPrim1("x3",a,b))
y1 = prim1ToMol1(ConstPrim1("y1",b,c))
y2 = prim1ToMol1(ConstPrim1("y2",b,c))
z1 = prim1ToMol1(ConstPrim1("z1",d,e))
z2 = prim1ToMol1(ConstPrim1("z2",d,e))

h = prim2ToAEMol2(ConstPrim2("h", x1, x2))
k = prim2ToAEMol2(ConstPrim2("k",y1,y2))
l = prim2ToAEMol2(ConstPrim2("l",z1,z2))

h1 = prim2ToAEMol2(ConstPrim2("h1", x1, x2))
h2 = prim2ToAEMol2(ConstPrim2("h2", x2, x3))

#aa1 = fEqMol1(comp0s([x,y,z]))
#aa2 = tensorDecompEqMol1(aa1)

#bb1 = comp1s([z,z,comp2s([h1,h2])])
#bb2 = idPrefixDecompAEMol2(bb1)
bb3 = fEqAEMol2(comp0s([w,h,w,l]))
bb4 = tensorDecompEqAEMol2(bb3)


cc1 = comp0s([unit, b, unit])
cc2 = b
cc3 = findPaths1(cc1, cc2, [unitor10, unitor11])
