from decomp import *
from ontology import *
from fill import *
from primfamily import *

a, b, c, d, e = ConstPrim0("a"), ConstPrim0("b"), ConstPrim0("c"), ConstPrim0("d"), ConstPrim0("e")
H = Functor("H")
Ha, Hb, Hc, Hd, He = app(H,a), app(H,b), app(H,c), app(H,d), app(H,e)

# We add the functor H in the third argument, indicating that it is a functor variable.
chiSource = minimalASTFromMol0(comp0(app(H, a), app(H,b)), [a, b], [H])
chiTarget = minimalASTFromMol0(ensureMol0(app(H, comp0(a, b))), [a, b], [H])
# We also set the third argument of chi to 1, indicating that its source and target contain a single functor variable.
chi = PrimitiveFamily("chi", 1, 1, [0,0], chiSource, chiTarget)

omega2Source0 = ensureMol0(comp0s(Ha, Hb, Hc))
omega2Target0 = ensureMol0(app(H, comp0s(a, b, c)))
omega2Paths1 = findPaths1(omega2Source0, omega2Target0, [chi])

# lmap(str, omega2Paths1)

omega2Source1AST = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(Ha, chi.fprimf([H], b, c)),
        chi.fprimf([H], a, comp0(b, c)))), [a, b, c], [H])
omega2Target1AST = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(chi.fprimf([H], a, b), Hc),
        chi.fprimf([H], comp0(a, b), c))), [a, b, c], [H])
omega2 = PrimitiveFamily("omega2", 2, 1, [0, 0, 0], omega2Source1AST, omega2Target1AST)



omega3Source0 = ensureMol0(comp0s(Ha, Hb, Hc, Hd))
omega3Target0 = ensureMol0(app(H, comp0s(a, b, c, d)))
omega3Paths1 = findPaths1(omega3Source0, omega3Target0, [chi])

# len(omega3Paths1)

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
