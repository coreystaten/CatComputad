from decomp import *
from ontology import *
from fill import *
from primfamily import *

###############
# LISTING 1
a0 = ConstPrim0("a")
b0 = ConstPrim0("b")

f1 = ConstPrim1("f", a0, b0)
g1 = ConstPrim1("g", a0, b0)

u2 = ConstPrim2("u", f1, g1)
v2 = ConstPrim2("v", f1, g1)

p3 = ConstPrim3("p", u2, v2)

###############
# LISTING 2
y1 = prim0ToMol0(a0)

x1 = fMol0((a0, b0))
x2 = comp0(a0, b0)
# x1 == x2
# >>> True


###############
# LISTING 3
a, b, c, d = ConstPrim0("a"), ConstPrim0("b"), ConstPrim0("c"), ConstPrim0("d")
f, g = ConstPrim1("f", a, b), ConstPrim1("g", c, d)
emptyMol0 = fMol0(())

# Constructing an identity molecule.
idMol = fIdMol1(ensureMol0(a))
# Constructing the same equivalence class in two ways.
idMolEq = fEqMol1(idMol)
idMolEq2 = ensureEqMol1(a) # This automatically upcasts a to an identity molecule.
# idMolEq == idMolEq2
# >>> True

# Constructing a pair of atoms in two different ways.
firstAtom = fAtom1(emptyMol0, f, ensureMol0(c))
secondAtom = fAtom1(ensureMol0(b), g, emptyMol0)
# str(firstAtom)
# >>> 'f @ c'
# str(secondAtom)
# >>> 'b @ g'
firstAtom2 = comp0(f, c)
secondAtom2 = comp0(b, g)
# (firstAtom, secondAtom) == (firstAtom2, secondAtom2)
# >>> True

# Building a 1-molecule in two different ways.
mol = fNonIdMol1((firstAtom, secondAtom))
mol2 = comp1(firstAtom, secondAtom)
# mol == mol2
# >>> True

# Finding the equivalence class of a 1-molecule.
eqMol = fEqMol1(mol)
# print(lmap(str, eqMol.mol1s))
# >>> ['(f @ c) . (b @ g)', '(a @ g) . (f @ d)']

###############
# LISTING 3A
# f = ConstPrim1("f", fMol0(()), fMol0(()))


###############
# LISTING 4
a, b, c, d = ConstPrim0("a"), ConstPrim0("b"), ConstPrim0("c"), ConstPrim0("d")
f, g = ConstPrim1("f", a, b), ConstPrim1("g", a, b)
h, k = ConstPrim1("h", c, d), ConstPrim1("k", c, d)
u, v = ConstPrim2("u", f, g), ConstPrim2("v", h, k)

# Checking the equivalence class of a 2-atom.
atom = comp0(u, h)
eqAtom = fEqAtom2(atom)
# print(lmap(str, eqAtom.atom2s))
# >>> ['(u @ c) . (b @ h)', '(a @ h) . (u @ d)']

# Checking the equivalence class of a 2-molecule.
# Note that & is used to represent 2-composition.
mol = comp0(u, v)
eqMol = ensureEqAEMol2(mol)
# print(lmap(str, eqMol.aeMol2s))
# >>> ['[(u @ c) . (b @ h)] & [(a @ v) . (g @ d)]', '[(a @ v) . (f @ d)] & [(u @ c) . (b @ k)]']


###############
# LISTING 5
def displayDecomps(decomps):
    print(lmap(lambda x: (str(x[0]), str(x[1])), decomps))

td = tensorDecompEqAEMol2(ensureEqAEMol2(comp0(u,v)))
displayDecomps(td)

hd = horizontalDecompEqAEMol2(ensureEqAEMol2(comp0(u,v)))
displayDecomps(hd)

vd = verticalDecompEqAEMol2(ensureEqAEMol2(comp0(u,v)))
displayDecomps(vd)

###############
# LISTING 6
unit = ConstPrim0("unit")
a = ConstPrim0("a")
lSource = minimalASTFromMol0(comp0(unit, a), [a])
lTarget = minimalASTFromMol0(ensureMol0(a), [a])
rSource = minimalASTFromMol0(comp0(a, unit), [a])
rTarget = minimalASTFromMol0(ensureMol0(a), [a])

u1_L = PrimitiveFamily("unitor1L", 1, 0, [0], lSource, lTarget)
u1_R = PrimitiveFamily("unitor1R", 1, 0, [0], rSource, rTarget)

###############
# LISTING 7
b = ConstPrim0("b")
source0 = comp0s(unit, a, b)
target0 = comp0s(a, b)
paths = findPaths1(source0, target0, [u1_L, u1_R])
pathsXXX = paths
lmap(str, pathsXXX)
u2_0Source = minimalASTFromEqMol1(paths[0], [a, b])
u2_0Target = minimalASTFromEqMol1(paths[1], [a, b])
u2_0 = PrimitiveFamily("unitor2_0", 2, 0, [0, 0], u2_0Source, u2_0Target)

source1 = comp0s(a, unit, b)
target1 = comp0s(a, b)
paths = findPaths1(source1, target1, [u1_L, u1_R])
u2_1Source = minimalASTFromEqMol1(paths[0], [a, b])
u2_1Target = minimalASTFromEqMol1(paths[1], [a, b])
u2_1 = PrimitiveFamily("unitor2_2", 2, 0, [0, 0], u2_1Source, u2_1Target)

source2 = comp0s(a, b, unit)
target2 = comp0s(a, b)
paths = findPaths1(source2, target2, [u1_L, u1_R])
u2_2Source = minimalASTFromEqMol1(paths[0], [a, b])
u2_2Target = minimalASTFromEqMol1(paths[1], [a, b])
u2_2 = PrimitiveFamily("unitor2_2", 2, 0, [0, 0], u2_2Source, u2_2Target)


###############
# LISTING 8
c = ConstPrim0("c")
sourceDim0 = comp0s(unit, a, b, c)
targetDim0 = comp0s(a, b, c)
paths = findPaths1(sourceDim0, targetDim0, [u1_L, u1_R])

for ii in range(len(paths)):
    for jj in range(len(paths)):
            if ii != jj:
                print("(%d, %d) - %d paths" % (ii, jj, len(findPaths2(paths[ii], paths[jj], [u2_0, u2_1, u2_2]))))

###############
# LISTING 9

#pathsDim2 = findPaths2(paths[0], paths[1], [u2_0, u2_1, u2_2])
#str(pathsDim2[0])
#str(pathsDim2[1])
#u3_0Source = minimalASTFromEqAEMol2(pathsDim2[0], [a, b, c])
#u3_0Target = minimalASTFromEqAEMol2(pathsDim2[1], [a, b, c])
#u3_0 = PrimitiveFamily("unitor3_0", 3, 0, [0, 0, 0], u3_0Source, u3_0Target)


###############
# LISTING 10
print("----")
for ii in range(len(paths)):
    for jj in range(len(paths)):
            if ii != jj:
                print("(%d, %d) - %d paths" % (ii, jj, len(findPaths2(paths[ii], paths[jj], [u2_0, u2_1, u2_2, u2_0.adj, u2_1.adj, u2_2.adj]))))

###############
# LISTING 10
a, b, c, d, e = ConstPrim0("a"), ConstPrim0("b"), ConstPrim0("c"), ConstPrim0("d"), ConstPrim0("e")
H = Functor("H")

# We add the functor H in the third argument, indicating that it is a functor variable.
chiSource = minimalASTFromMol0(comp0(app(H, a), app(H, b)), [a, b], [H])
chiTarget = minimalASTFromMol0(ensureMol0(app(H, comp0(a, b))), [a, b], [H])
# We also set the third argument of chi to 1, indicating that its source and target contain a single functor variable.
chi = PrimitiveFamily("chi", 1, 1, [0,0], chiSource, chiTarget)

###############
# LISTING 11
Ha, Hb, Hc, Hd, He = app(H,a), app(H,b), app(H,c), app(H,d), app(H,e)

omega2Source0 = ensureMol0(comp0s(Ha, Hb, Hc))
omega2Target0 = ensureMol0(app(H, comp0s(a, b, c)))
omega2Paths1 = findPaths1(omega2Source0, omega2Target0, [chi])

lmap(str, omega2Paths1)
# >>> ['[(chi{H}(a, b) @ H(c)) . (chi{H}((a @ b), c))]', '[(H(a) @ chi{H}(b, c)) . (chi{H}(a, (b @ c)))]']

# The .fprimf method is used to get an instance of the primitive family with the given parameters.
omega2Source1AST = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(Ha, chi.fprimf([H], b, c)),
        chi.fprimf([H], a, comp0(b, c)))), [a, b, c], [H])
omega2Target1AST = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(chi.fprimf([H], a, b), Hc),
        chi.fprimf([H], comp0(a, b), c))), [a, b, c], [H])
omega2 = PrimitiveFamily("omega2", 2, 1, [0, 0, 0], omega2Source1AST, omega2Target1AST)


###############
# LISTING 12
omega3Source0 = ensureMol0(comp0s(Ha, Hb, Hc, Hd))
omega3Target0 = ensureMol0(app(H, comp0s(a, b, c, d)))
omega3Paths1 = findPaths1(omega3Source0, omega3Target0, [chi])
len(omega3Paths1)
# >>> 5

# Some filtering later...
# Source: [((H(a) @ H(b)) @ chi1{0}(c, d)) . (H(a) @ chi1{0}(b, (c @ d))) . (chi1{0}(a, (b @ c @ d)))]
# --
# Target: [(chi1{H}(a, b) @ (H(c) @ H(d))) . (chi1{H}((a @ b), c) @ H(d)) . (chi1{0}((a @ b @ c), d))]
omega3Source1 = ensureEqMol1(
    comp1s(
        comp0s(Ha, Hb, chi.fprimf([H], c, d)),
        comp0(Ha, chi.fprimf([H], b, comp0(c, d))),
        chi.fprimf([H], a, comp0s(b, c, d))))
omega3Target1 = ensureEqMol1(
    comp1s(
        comp0s(chi.fprimf([H], a, b), Hc, Hd),
        comp0(chi.fprimf([H], comp0(a, b), c), Hd),
        chi.fprimf([H], comp0s(a, b, c), d)))
omega3Paths2 = findPaths2(omega3Source1, omega3Target1, [omega2])
len(omega3Paths2)
# >>> 2
